#!/usr/bin/env python3
"""Verify the assembled DAVE-ML model by evaluating it independently.

A generated file that nothing re-reads is the one artifact in a codebase
that can rot silently, so this parses the .dml from scratch -- it shares
no code with build-daveml.py deliberately, because an evaluator that
reused the writer's own helpers would agree with it by construction --
and then:

  1. checks structural invariants (every function's breakpoints exist and
     its table has exactly the product of their lengths, every varID is
     unique, every referenced varID is defined);
  2. evaluates each gridded table by multilinear interpolation;
  3. runs every checkData staticShot and compares against its stated
     tolerance.

No external DAVE-ML implementation is used. The MathML subset the file
emits is under our control precisely so a compact reader suffices, and
the alternative -- depending on a large third-party toolchain -- cuts
against the same instinct that keeps BEMT and VPM out of this repo.

Exit status is 0 only if every check passes, so this is usable directly
as a ctest.

Usage:  python verify-daveml.py [model.dml]
"""

import os
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT = os.path.join(HERE, "AetherionFlightModel.dml")

failures = []
checks = 0


def check(condition, message):
    global checks
    checks += 1
    if not condition:
        failures.append(message)
    return condition


def parse_numbers(text):
    """A DAVE-ML value list: comma and/or whitespace separated."""
    if text is None:
        return []
    return [float(tok) for tok in text.replace(",", " ").split()]


def interp(bp_values, table, point):
    """Multilinear interpolation on a rectilinear grid.

    bp_values : list of ascending breakpoint arrays, one per dimension
    table     : flat values, row-major (last dimension fastest)
    point     : one coordinate per dimension
    Out-of-range coordinates clamp to the edge, which is what a gridded
    table means -- extrapolation is the assembler's business, not the
    lookup's.
    """
    dims = len(bp_values)
    lo_idx, frac = [], []
    for d in range(dims):
        axis, x = bp_values[d], point[d]
        if len(axis) == 1:
            lo_idx.append(0)
            frac.append(0.0)
            continue
        if x <= axis[0]:
            lo_idx.append(0)
            frac.append(0.0)
        elif x >= axis[-1]:
            lo_idx.append(len(axis) - 2)
            frac.append(1.0)
        else:
            i = 0
            while i + 2 < len(axis) and axis[i + 1] < x:
                i += 1
            lo_idx.append(i)
            span = axis[i + 1] - axis[i]
            frac.append((x - axis[i]) / span if span > 0 else 0.0)

    # Strides for row-major order.
    strides = [1] * dims
    for d in range(dims - 2, -1, -1):
        strides[d] = strides[d + 1] * len(bp_values[d + 1])

    total = 0.0
    for corner in range(1 << dims):
        weight = 1.0
        offset = 0
        for d in range(dims):
            high = (corner >> d) & 1
            if len(bp_values[d]) == 1:
                if high:
                    weight = 0.0
                    break
                idx = 0
            else:
                idx = lo_idx[d] + high
                weight *= frac[d] if high else (1.0 - frac[d])
            offset += idx * strides[d]
        if weight == 0.0:
            continue
        if offset >= len(table):
            return None
        total += weight * table[offset]
    return total


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    if not os.path.exists(path):
        print(f"no model at {path} -- run build-daveml.py first", file=sys.stderr)
        return 1

    try:
        root = ET.parse(path).getroot()
    except ET.ParseError as err:
        print(f"FAIL: {path} is not well-formed XML: {err}", file=sys.stderr)
        return 1
    print(f"parsed {os.path.basename(path)}: well-formed")

    # DTD validation, where a validator exists. This is not redundant with
    # the structural checks below: the DTD encodes the STANDARD's grammar,
    # and it caught six classes of real deviation the first time it was
    # run -- element ordering at the top level, initialValue being an
    # attribute rather than an element, reference/documentRef ID discipline,
    # and the true shape of a gridded function. The in-repo checks then
    # verify what a DTD structurally cannot: that the numbers interpolate.
    #
    # lxml is deliberately not required; xmllint is used when present and
    # the step reports honestly when it is not, rather than passing quietly.
    dtd = os.path.join(os.path.dirname(os.path.abspath(path)), "DAVEfunc.dtd")
    if not os.path.exists(dtd):
        print("  DTD: models/DAVEfunc.dtd absent, skipping grammar validation")
    elif shutil.which("xmllint") is None:
        print("  DTD: no xmllint on PATH, skipping grammar validation "
              "(structural checks below still run)")
    else:
        proc = subprocess.run(["xmllint", "--noout", "--dtdvalid", dtd, path],
                              capture_output=True, text=True)
        if proc.returncode == 0:
            print("  DTD: validates against DAVEfunc.dtd 2.0.1")
        else:
            head = "\n".join(proc.stderr.strip().split("\n")[:12])
            check(False, f"DTD validation failed:\n{head}")

    # ---------------- structure ----------------
    var_ids = [v.get("varID") for v in root.iter("variableDef")]
    check(len(var_ids) == len(set(var_ids)),
          f"varIDs must be unique; {len(var_ids) - len(set(var_ids))} duplicate(s)")
    defined = set(var_ids)

    bps = {}
    for bp in root.iter("breakpointDef"):
        vals = parse_numbers(bp.findtext("bpVals"))
        check(len(vals) >= 1, f"breakpoint {bp.get('bpID')} is empty")
        check(all(b < a for b, a in zip(vals, vals[1:])),
              f"breakpoint {bp.get('bpID')} must be strictly ascending")
        bps[bp.get("bpID")] = vals

    tables = {}
    for gt in root.iter("griddedTableDef"):
        refs = [r.get("bpID") for r in gt.iter("bpRef")]
        vals = parse_numbers(gt.findtext("dataTable"))
        expected = 1
        ok = True
        for r in refs:
            if r not in bps:
                check(False, f"table {gt.get('name')} references unknown breakpoint {r}")
                ok = False
                break
            expected *= len(bps[r])
        if not ok:
            continue
        check(len(vals) == expected,
              f"table {gt.get('name')}: {len(vals)} values against "
              f"{expected} from its breakpoints {refs}")
        check(all(v == v for v in vals), f"table {gt.get('name')} contains NaN")
        tables[gt.get("gtID")] = (refs, vals)

    print(f"  {len(defined)} variables, {len(bps)} breakpoint sets, {len(tables)} tables")

    # ---------------- functions ----------------
    funcs = {}
    for fn in root.iter("function"):
        dep = fn.find("dependentVarRef")
        defn = fn.find("functionDefn")
        if dep is None or defn is None:
            continue
        out = dep.get("varID")
        # The DTD form: functionDefn CONTAINS a griddedTableRef rather than
        # carrying a gtID, and independentVarRef carries no bpID at all --
        # the breakpoint association comes from the table's own
        # breakpointRefs, matched POSITIONALLY against these.
        ref_el = defn.find("griddedTableRef")
        gtid = ref_el.get("gtID") if ref_el is not None else None
        check(out in defined, f"function {fn.get('name')} drives undefined varID {out}")
        check(gtid in tables, f"function {fn.get('name')} references unknown table {gtid}")
        in_vars = [iv.get("varID") for iv in fn.iter("independentVarRef")]
        for v in in_vars:
            check(v in defined,
                  f"function {fn.get('name')} reads undefined varID {v}")
        if gtid in tables:
            refs, _ = tables[gtid]
            check(len(in_vars) == len(refs),
                  f"function {fn.get('name')}: {len(in_vars)} independent variables against "
                  f"{len(refs)} table dimensions")
            if len(in_vars) == len(refs):
                funcs[out] = (list(zip(in_vars, refs)), gtid)

    # Every table must be evaluable at the centre of its own grid: a
    # cheap, total check that the encoding and the lookup agree.
    for out, (ins, gtid) in funcs.items():
        refs, vals = tables[gtid]
        mid = [bps[b][len(bps[b]) // 2] for _, b in ins]
        got = interp([bps[b] for _, b in ins], vals, mid)
        check(got is not None and got == got,
              f"{out} is not evaluable at the centre of its grid")

    # ---------------- checkData ----------------
    shots = list(root.iter("staticShot"))
    print(f"  {len(shots)} staticShot(s)")
    for shot in shots:
        name = shot.get("name")
        inputs = {}
        for sig in shot.findall("./checkInputs/signal"):
            inputs[sig.findtext("signalName")] = float(sig.findtext("signalValue"))
        for sig in shot.findall("./checkOutputs/signal"):
            out = sig.findtext("signalName")
            want = float(sig.findtext("signalValue"))
            tol = float(sig.findtext("tol") or 1e-9)
            if out not in funcs:
                check(False, f"shot {name}: no function drives {out}")
                continue
            ins, gtid = funcs[out]
            refs, vals = tables[gtid]
            missing = [v for v, _ in ins if v not in inputs]
            if missing:
                check(False, f"shot {name}: inputs {missing} not supplied for {out}")
                continue
            got = interp([bps[b] for _, b in ins], vals, [inputs[v] for v, _ in ins])
            ok = got is not None and abs(got - want) <= tol + 1e-12
            check(ok, f"shot {name}: {out} evaluated {got}, expected {want} (tol {tol})")
            if ok:
                print(f"    {name}: {out} = {got:.9g} OK")

    # ---------------- report ----------------
    print()
    if failures:
        for f in failures:
            print(f"FAIL: {f}", file=sys.stderr)
        print(f"{len(failures)} of {checks} checks failed", file=sys.stderr)
        return 1
    print(f"PASS: all {checks} checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
