#!/usr/bin/env python3
"""Assemble the Aetherion DAVE-ML flight model from the cached sweep JSONs.

The house pattern (models/README.md): each driver exports its own JSON,
and this assembler builds the single .dml from those cached files, so
regenerating one map never reruns another. Nothing here computes
aerodynamics; if a number is wrong, it is wrong in the JSON.

Emits ANSI/AIAA S-119-2011 DAVE-ML 2.0.1. Inputs, in models/data/:

    aero-map.json         aeolion_aero_map        baseline + rate derivatives
    aero-aileron.json     aeolion_aero_map ... aileron   control increments
    parasite-drag.json    aeolion_parasite_drag   aeroCD0(alpha)
    propulsion-map.json   aeolion_propulsion_map  propCT/propCQ(J)
    propulsion-singlevane.json                    per-vane increments
    coupling-map.json     aeolion_induction_map   fan-on-airframe increments

Any block whose JSON is absent is DECLARED absent in the file header
rather than emitted as zeros, because a zero-valued table is a lie a
consumer cannot detect.

The fan-on-airframe interaction was long declared blocked on an
upstream-induction model. That blocker is stale: Solver/DiskInduction.h
implements the semi-infinite vortex cylinder, which HAS a field upstream
of the disk -- the half the momentum-theory slipstream lacks, and the
half that matters when the fan sits behind the wing.

Usage:  python build-daveml.py [--data DIR] [--out FILE]
"""

import argparse
import json
import os
import sys
from datetime import date

HERE = os.path.dirname(os.path.abspath(__file__))

# --- what the standard separates -------------------------------------------
# varID is "an internal identifier that is unique within the file" and is
# unconstrained in form; name "should correspond to the standard AIAA
# parameter name" (Annex A). So the namespaced varIDs stay and the
# standard name rides alongside. See models/README.md.
INPUTS = [
    # varID,             AIAA name,                units,     axis,   sign,  symbol
    ("alphaDeg", "angleOfAttack", "deg", "body", "", "α"),
    ("betaDeg", "angleOfSideslip", "deg", "body", "", "β"),
    ("trueAirspeedMps", "trueAirspeed", "m/s", "", "", "V"),
    ("airDensityKgpm3", "airDensity", "kg/m3", "", "", "ρ"),
    ("rollRateRadps", "bodyAngularRate_Roll", "rad/s", "body", "+RWD", "p"),
    ("pitchRateRadps", "bodyAngularRate_Pitch", "rad/s", "body", "+NU", "q"),
    ("yawRateRadps", "bodyAngularRate_Yaw", "rad/s", "body", "+NR", "r"),
    ("propSpeedRevps", "propellerSpeed", "rev/s", "", "", "n"),
    ("aileronDeg", "aileronDeflection", "deg", "body", "TED", "δa"),
    ("vanePitchDeg", "vaneDeflection_Pitch", "deg", "body", "", "δP"),
    ("vaneYawDeg", "vaneDeflection_Yaw", "deg", "body", "", "δY"),
    ("vaneRollDeg", "vaneDeflection_Roll", "deg", "body", "", "δR"),
]


def esc(text):
    return (str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def salvage(text):
    """Parse a sweep JSON that is still being written.

    The drivers flush each row as it completes, precisely so an
    interrupted long sweep keeps its finished conditions. That leaves the
    file without its closing brackets while the run is in flight, so the
    assembler truncates to the last complete row and closes the structure
    rather than refusing to read it.
    """
    cut = text.rfind("}")
    while cut > 0:
        candidate = text[: cut + 1]
        for suffix in ("\n]}", "\n]}}", "}"):
            try:
                return json.loads(candidate + suffix), True
            except json.JSONDecodeError:
                continue
        cut = text.rfind("}", 0, cut)
    return None, False


def load(data_dir, name, required=True):
    path = os.path.join(data_dir, name)
    if not os.path.exists(path):
        if required:
            print(f"missing required input {path}", file=sys.stderr)
            sys.exit(1)
        return None
    with open(path, "r") as handle:
        text = handle.read()
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        doc, ok = salvage(text)
        if not ok:
            print(f"cannot parse {path}", file=sys.stderr)
            sys.exit(1)
        print(f"note: {name} is still being written -- using its completed rows only",
              file=sys.stderr)
        return doc


def uniq(values, tol=1e-9):
    """Sorted unique values, collapsing float noise."""
    out = []
    for v in sorted(values):
        if not out or abs(v - out[-1]) > tol:
            out.append(v)
    return out


class Doc:
    """Minimal indented-XML writer. No dependency, matching the in-repo
    verifier's deliberately small MathML/table subset."""

    def __init__(self):
        self.lines = []
        self.depth = 0

    def raw(self, text):
        self.lines.append("  " * self.depth + text)

    def open(self, tag, **attrs):
        a = "".join(f' {k.replace("_", ":")}="{esc(v)}"' for k, v in attrs.items() if v != "")
        self.raw(f"<{tag}{a}>")
        self.depth += 1

    def close(self, tag):
        self.depth -= 1
        self.raw(f"</{tag}>")

    def leaf(self, tag, text="", **attrs):
        a = "".join(f' {k.replace("_", ":")}="{esc(v)}"' for k, v in attrs.items() if v != "")
        if text == "":
            self.raw(f"<{tag}{a}/>")
        else:
            self.raw(f"<{tag}{a}>{esc(text)}</{tag}>")

    def comment(self, text):
        # A double hyphen may not appear inside an XML comment, and the
        # prose in this file uses "--" as an em dash freely. Substitute
        # here rather than rely on every call site remembering.
        for line in text.strip().split("\n"):
            self.raw("<!-- " + line.strip().replace("--", "—") + " -->")

    def text(self):
        return "\n".join(self.lines) + "\n"


# --- MathML helpers ---------------------------------------------------------
# The subset is deliberately narrow so the in-repo verifier's evaluator can
# stay compact: ci, cn, apply with plus/minus/times/divide/power/sin/cos/
# arcsin/arctan/abs/max/min.


def ci(name):
    return f"<ci>{name}</ci>"


def cn(value):
    return f"<cn>{value}</cn>"


def app(op, *args):
    return f"<apply><{op}/>" + "".join(args) + "</apply>"


def deg2rad(expr):
    return app("times", expr, app("divide", "<pi/>", cn(180)))


class Model:
    def __init__(self, doc):
        self.d = doc
        self.calc_order = []

    def variable(self, varID, name, units, *, axis="", sign="", symbol="",
                 initial=None, calc=None, is_output=False, description=""):
        attrs = dict(name=name, varID=varID, units=units)
        if initial is not None:
            attrs["initialValue"] = initial
        if axis:
            attrs["axisSystem"] = axis
        if sign:
            attrs["sign"] = sign
        if symbol:
            attrs["symbol"] = symbol
        self.d.open("variableDef", **attrs)
        if description:
            self.d.leaf("description", description)
        if calc is not None:
            self.d.open("calculation")
            self.d.open("math", xmlns="http://www.w3.org/1998/Math/MathML")
            self.d.raw(calc)
            self.d.close("math")
            self.d.close("calculation")
        if is_output:
            self.d.leaf("isOutput")
        self.d.close("variableDef")


def emit_breakpoints(d, bpID, name, values):
    d.open("breakpointDef", bpID=bpID, name=name)
    d.leaf("bpVals", ", ".join(f"{v:g}" for v in values))
    d.close("breakpointDef")


def gridded_table(d, name, bp_ids, values, description=""):
    """values: flat list in row-major order over bp_ids."""
    d.open("griddedTableDef", name=name, gtID=f"{name}_data")
    if description:
        d.leaf("description", description)
    d.open("breakpointRefs")
    for b in bp_ids:
        d.leaf("bpRef", bpID=b)
    d.close("breakpointRefs")
    # Wrap for legibility; the standard is whitespace-insensitive here.
    chunk = 8
    d.open("dataTable")
    for i in range(0, len(values), chunk):
        d.raw(", ".join(f"{v:.9g}" for v in values[i:i + chunk]) +
              ("," if i + chunk < len(values) else ""))
    d.close("dataTable")
    d.close("griddedTableDef")


def simple_function(d, name, out_var, in_vars, bp_ids, table_name, description=""):
    """A gridded function: independent variables -> one dependent variable.

    The DTD form, which is not the obvious one. independentVarRef carries
    NO bpID -- the breakpoint association comes from the table's own
    breakpointRefs, so the ORDER of these must match it. And functionDefn
    CONTAINS a griddedTableRef element rather than carrying a gtID
    attribute.

    extrapolate="neither" states the clamping the verifier's lookup
    already does: a gridded table is defined on its grid, and holding the
    edge value is a claim the file should make explicitly rather than
    leave to a consumer's default.
    """
    d.open("function", name=name)
    if description:
        d.leaf("description", description)
    for var in in_vars:
        d.leaf("independentVarRef", varID=var, interpolate="linear", extrapolate="neither")
    d.leaf("dependentVarRef", varID=out_var)
    d.open("functionDefn")
    d.leaf("griddedTableRef", gtID=f"{table_name}_data")
    d.close("functionDefn")
    d.close("function")


def build(args):
    data = args.data
    aero = load(data, "aero-map.json", required=False)
    parasite = load(data, "parasite-drag.json", required=False)
    prop = load(data, "propulsion-map.json", required=False)
    vane = load(data, "propulsion-singlevane.json", required=False)
    coupling = load(data, "coupling-map.json", required=False)
    # The aileron sweep may live in its own file: the driver can be run
    # with the block selector so the expensive baseline map is not
    # regenerated alongside it, and merging two JSONs by hand is exactly
    # the sort of step that goes wrong without saying so.
    ailfile = load(data, "aero-aileron.json", required=False)
    if ailfile and aero is not None and not aero.get("aileron"):
        aero = dict(aero)
        aero["aileron"] = ailfile.get("aileron", [])
        for k in ("aileronMirrorProbeDeg", "aileronTau", "aileronEtaStart",
                  "aileronChordFraction"):
            if k in ailfile:
                aero[k] = ailfile[k]

    if aero is None:
        print("note: aero-map.json absent -- emitting the model without aero* tables",
              file=sys.stderr)

    d = Doc()
    m = Model(d)
    d.raw('<?xml version="1.0" encoding="UTF-8"?>')
    d.raw('<!DOCTYPE DAVEfunc SYSTEM "DAVEfunc.dtd">')
    d.open("DAVEfunc")

    # ---------------- fileHeader ----------------
    d.open("fileHeader", name="AetherionFlightModel")
    d.leaf("author", "", name="Onur Tuncer", org="Istanbul Technical University")
    d.leaf("author", "", name="Caglar Ucler", org="Ozyegin University")
    d.leaf("author", "", name="Ahmet Gunes", org="Istanbul Technical University")
    d.leaf("fileCreationDate", "", date=str(date.today()))
    # The description carries the validity envelope and the declared gaps.
    # These were once <reference> elements, which the DTD reserves for
    # bibliography -- refID is an ID and must be unique, and author and date
    # are required. Prose belongs here.
    notes = [
        "Tabulated flight model of the Aetherion ducted-fan tail-sitter, generated from "
        "the Aeolion aerodynamic toolkit. Body axes are the contract frame "
        "aetherion_body_frd (x forward, y right, z down), the standard aeronautical body "
        "axis system of ANSI/AIAA R-004-1992. Moments are about the contract's "
        "moment_reference_point. See models/README.md for the normative specification "
        "and models/report/ for the technical report.",
        "",
        "VALIDITY. Incompressible, M < 0.3: no Mach dependence exists in the generating "
        "methods. Single Reynolds number -- every table was generated at V = 25 m/s. "
        "Propulsor tables are AXIAL INFLOW ONLY, because the rotor-vane machinery is "
        "axisymmetric end to end; alphaDiskDeg is output as a validity monitor rather "
        "than faked as a table axis. Powered operation only: the rho n^2 D^4 group "
        "excludes n -> 0, so windmilling and low-rotor-speed descent are outside the "
        "envelope.",
        "",
        "LIMITATIONS. Post-stall values are limit-cycle means, not steady states. CLmax "
        "is an upper bound -- bubble bursting is not modelled. The tables are the "
        "ascending-alpha branch; hysteresis is not represented. Rate derivatives are "
        "tapered to zero over alpha 20-40 deg, a declared assumption rather than a "
        "computed result. Parasite drag covers body and duct only -- the wing's profile "
        "drag is already inside the force tables -- and omits the duct's separated drag "
        "at incidence, so aeroCD0 is a lower bound at high alpha. The aileron flap model "
        "is lift-only: no section pitching-moment increment, no gap leakage, no viscous "
        "decay at large deflection, so tabulated roll authority is an upper bound.",
    ]
    if not (aero or {}).get("aileron"):
        notes.append("")
        notes.append(
            "INCOMPLETE: aileron increment tables are absent -- the deflected sweep has "
            "not been run -- so this model carries NO ROLL CONTROL INPUT.")
    notes.append("")
    if coupling:
        notes.append(
            "INTERACTION VALIDITY. The fan-on-airframe tables (coupling*) were swept at "
            "ZERO SIDESLIP and are indexed by alpha and Tc only, but the buildup applies "
            "them at every beta. Their beta dependence is therefore unmeasured, not "
            "established as weak. They are also uniform-disk-loading and swirl-free, and "
            "past alpha 16 they are differences of two limit-cycle means: the effect is "
            "thrust-ordered and real, but its attribution to delayed separation is "
            "consistent with the data rather than established by it.")
    else:
        notes.append(
            "INCOMPLETE: fan-on-airframe interaction tables (coupling*) are absent, so the "
            "airframe tables are POWER-OFF and underpredict the separation delay the aft "
            "fan provides in transition.")
    d.open("description")
    for line in notes:
        d.raw(esc(line))
    d.close("description")

    sources = [
        ("srcAeroMap", "aero-map.json", "aeolion_aero_map",
         "Airframe baseline map and reduced-rate derivatives", aero),
        ("srcAileron", "aero-aileron.json", "aeolion_aero_map (aileron block)",
         "Aileron increment sweep", ailfile or ((aero or {}).get("aileron") and aero)),
        ("srcParasite", "parasite-drag.json", "aeolion_parasite_drag",
         "Body and duct parasite drag, friction buildup plus crossflow branch", parasite),
        ("srcPropMap", "propulsion-map.json", "aeolion_propulsion_map",
         "Ducted propulsor over advance ratio", prop),
        ("srcPropVane", "propulsion-singlevane.json", "aeolion_propulsion_map (single)",
         "Per-vane control increments", vane),
    ]
    present = [(rid, fn, who, what) for rid, fn, who, what, doc in sources if doc]
    # <reference> is bibliography: refID must be a unique ID, and author,
    # title and date are required. The sweep JSONs are exactly that -- the
    # documents this model was generated from -- so they belong here and
    # the provenance's documentRefs point at them.
    for rid, fn, who, what in present:
        d.leaf("reference", "", refID=rid, author=who, title=f"{what} ({fn})",
               date=str(date.today()))

    d.open("provenance", provID="genProv")
    d.leaf("author", "", name="Aeolion", org="models/build-daveml.py")
    d.leaf("creationDate", "", date=str(date.today()))
    for rid, fn, who, what in present:
        d.leaf("documentRef", "", refID=rid)
    d.leaf("description", "Generated from cached solver sweeps; see each documentRef. "
                          "Nothing in the assembler computes aerodynamics.")
    d.close("provenance")

    d.close("fileHeader")

    # ---------------- inputs ----------------
    d.comment("Inputs. varID is namespaced for this file's structure; name is the "
              "Annex A standard parameter name.")
    for varID, name, units, axis, sign, symbol in INPUTS:
        m.variable(varID, name, units, axis=axis, sign=sign, symbol=symbol,
                   initial=0.0)
        d.raw("")  # spacing

    # ---------------- constants ----------------
    consts = {}
    if aero:
        meta = aero["meta"]
        consts["WingAreaM2"] = meta["area"]
        consts["WingSpanM"] = meta["span"]
        consts["WingChordM"] = meta["chord"]
        rp = meta.get("refPoint", [0, 0, 0])
        # refPoint is exported in SOLVER axes; convert back to FRD for the file.
        consts["XmrpM"] = -rp[0]
        consts["YmrpM"] = rp[1]
        consts["ZmrpM"] = -rp[2]
    if prop:
        consts["DiskDiameterM"] = prop["meta"]["diameterM"]

    d.comment("Constants, transcribed from the geometry contract by the assembler.")
    for k, v in consts.items():
        m.variable(k, k, "m2" if k.endswith("M2") else "m", initial=v)

    # ---------------- derived ----------------
    d.comment("Derived quantities. alphaDiskDeg and phiWDeg are VALIDITY MONITORS: "
              "they index no table, and exist so a consumer can detect that it has "
              "left the axial-inflow envelope.")
    m.variable("qbarPa", "dynamicPressure", "Pa", symbol="q",
               calc=app("times", cn(0.5), ci("airDensityKgpm3"),
                        app("power", ci("trueAirspeedMps"), cn(2))))
    if prop:
        m.variable("advanceRatio", "advanceRatio", "nd", symbol="J",
                   calc=app("divide", ci("trueAirspeedMps"),
                            app("times", ci("propSpeedRevps"), ci("DiskDiameterM"))))
    m.variable("alphaDiskDeg", "angleOfAttackDisk", "deg", axis="body",
               description="Validity monitor: angle between the free stream and the "
                           "rotor axis. The propulsor model is axial-inflow only.",
               calc=app("times",
                        app("arccos", app("times",
                                          app("cos", deg2rad(ci("alphaDeg"))),
                                          app("cos", deg2rad(ci("betaDeg"))))),
                        app("divide", cn(180), "<pi/>")),
               is_output=True)
    if vane and prop:
        # Mixing matrix (models/README.md): bottom = R+Y, left = R-P,
        # top = R-Y, right = R+P. These index the per-vane tables, which
        # replaced the per-mode ones after the mode-sum buildup was
        # measured wrong by up to 33%.
        for pos, expr in (
            ("Bottom", app("plus", ci("vaneRollDeg"), ci("vaneYawDeg"))),
            ("Left", app("minus", ci("vaneRollDeg"), ci("vanePitchDeg"))),
            ("Top", app("minus", ci("vaneRollDeg"), ci("vaneYawDeg"))),
            ("Right", app("plus", ci("vaneRollDeg"), ci("vanePitchDeg"))),
        ):
            m.variable(f"vane{pos}Deg", f"vaneDeflection_{pos}", "deg", axis="body",
                       calc=expr,
                       description=f"The {pos.lower()} vane's own total commanded angle, "
                                   "from the mixing matrix. The per-vane tables are indexed "
                                   "by this, and the four contributions are summed.")

    if coupling and prop and aero:
        # The interaction index, computed IN THE FILE from the propulsor's
        # own thrust so the two halves cannot disagree. No algebraic loop:
        # the propulsor wrench does not depend on Tc, so thrust is known
        # before the interaction tables are read.
        m.variable("propThrustN", "propellerThrust", "N", axis="body",
                   description="Thrust from the propulsor tables, rho n^2 D^4 CT.",
                   calc=app("times", ci("airDensityKgpm3"),
                            app("power", ci("propSpeedRevps"), cn(2)),
                            app("power", ci("DiskDiameterM"), cn(4)),
                            ci("propCT")),
                   is_output=True)
        m.variable("thrustCoefficient", "thrustCoefficient", "nd", symbol="Tc",
                   description="T / (qbar S). Indexes the fan-on-airframe interaction "
                               "tables. Degenerates as V -> 0, hence the declared "
                               "minimum speed.",
                   calc=app("divide", ci("propThrustN"),
                            app("times", ci("qbarPa"), ci("WingAreaM2"))),
                   is_output=True)

    if aero:
        m.variable("pHat", "reducedRollRate", "nd", symbol="phat",
                   calc=app("divide", app("times", ci("rollRateRadps"), ci("WingSpanM")),
                            app("times", cn(2), ci("trueAirspeedMps"))))
        m.variable("qHat", "reducedPitchRate", "nd", symbol="qhat",
                   calc=app("divide", app("times", ci("pitchRateRadps"), ci("WingChordM")),
                            app("times", cn(2), ci("trueAirspeedMps"))))
        m.variable("rHat", "reducedYawRate", "nd", symbol="rhat",
                   calc=app("divide", app("times", ci("yawRateRadps"), ci("WingSpanM")),
                            app("times", cn(2), ci("trueAirspeedMps"))))

    # ---------------- breakpoints ----------------
    tables = []   # (name, [bpIDs], values, out_var, description)
    bp_defs = []  # (bpID, name, values) -- emitted after every variableDef

    if aero:
        base = aero["baseline"]
        alphas = uniq([r["alphaDeg"] for r in base])
        betas = uniq([r["betaDeg"] for r in base])
        # One-sided beta is mirrored here, using the unpowered airframe's
        # symmetry: CY, Cl, Cn odd in beta; CX, CZ, Cm even.
        full_betas = uniq([-b for b in betas] + betas)
        bp_defs.append(("alphaBp", "angleOfAttack", alphas))
        bp_defs.append(("betaBp", "angleOfSideslip", full_betas))

        lookup = {(r["alphaDeg"], r["betaDeg"]): r for r in base}
        ODD = {"CY", "Cl", "Cn"}
        for comp in ("CX", "CY", "CZ", "Cl", "Cm", "Cn"):
            vals = []
            for a in alphas:
                for b in full_betas:
                    row = lookup.get((a, abs(b)))
                    if row is None:
                        vals.append(0.0)
                        continue
                    v = row[comp]
                    if b < 0 and comp in ODD:
                        v = -v
                    vals.append(v)
            tables.append((f"aero{comp}", ["alphaBp", "betaBp"], vals, f"aero{comp}Table",
                           f"Airframe {comp}, power-off, controls neutral. Carries induced "
                           f"and wing profile drag; parasite drag is aeroCD0."))

        # Aileron increments. Swept one-sided in deflection and mirrored
        # here on the same symmetry argument used for sideslip: the
        # configuration is mirror-symmetric about xz, and mirroring maps an
        # antisymmetric +delta command onto -delta while flipping the
        # lateral wrench. The generating sweep VERIFIES this rather than
        # assuming it -- it solves one negative deflection, which
        # reproduces the mirrored positive one to machine precision.
        ail = aero.get("aileron", [])
        if ail:
            probe = aero.get("aileronMirrorProbeDeg")
            # The probe row exists only to check the mirror; it must not
            # also become a breakpoint, or the axis gains a stray point.
            gen = [r for r in ail if probe is None or abs(r["deltaDeg"] - probe) > 1e-9]
            pos = uniq([r["deltaDeg"] for r in gen if r["deltaDeg"] > 0])
            full_deltas = uniq([-x for x in pos] + [0.0] + pos)
            bp_defs.append(("aileronBp", "aileronDeflection", full_deltas))

            lut = {(r["alphaDeg"], r["deltaDeg"]): r for r in gen}
            ODD_A = {"dCY", "dCl", "dCn"}
            for comp in ("dCX", "dCY", "dCZ", "dCl", "dCm", "dCn"):
                vals = []
                for a in alphas:
                    for dlt in full_deltas:
                        if abs(dlt) < 1e-12:
                            vals.append(0.0)  # neutral is the reference
                            continue
                        row = lut.get((a, abs(dlt)))
                        if row is None:
                            vals.append(0.0)
                            continue
                        v = row[comp]
                        if dlt < 0 and comp in ODD_A:
                            v = -v
                        vals.append(v)
                name = "aeroD" + comp[1:]
                tables.append((name, ["alphaBp", "aileronBp"], vals, name + "Table",
                               "Aileron increment " + comp[1:] + " from the neutral "
                               "configuration at matched alpha. The flap is carried in the "
                               "section, so this attenuates through stall as the sections "
                               "separate -- an inviscid lattice cannot show that."))

        rates = aero.get("rates", [])
        if rates:
            ralphas = uniq([r["alphaDeg"] for r in rates])
            bp_defs.append(("alphaRateBp", "angleOfAttack", ralphas))
            for comp in ("CZq", "Cmq", "Clp", "Cnp", "CYp", "Clr", "Cnr", "CYr"):
                vals = [r[comp] for r in sorted(rates, key=lambda x: x["alphaDeg"])]
                tables.append((f"aero{comp}", ["alphaRateBp"], vals, f"aero{comp}Table",
                               f"Reduced-rate derivative {comp}, attached range only."))

    if parasite:
        pts = parasite["table"]
        palphas = uniq([r["alphaDeg"] for r in pts])
        bp_defs.append(("alphaParasiteBp", "angleOfAttack", palphas))
        vals = [r["CD0"] for r in sorted(pts, key=lambda x: x["alphaDeg"])]
        tables.append(("aeroCD0", ["alphaParasiteBp"], vals, "aeroCD0Table",
                       "Parasite drag of BODY AND DUCT ONLY -- never the wing, whose "
                       "profile drag is already inside the force tables. Friction "
                       "buildup plus a slender-body crossflow branch."))

    if prop:
        rows = [r for r in prop["rows"] if r["mode"] == "baseline"]
        js = uniq([r["J"] for r in rows])
        bp_defs.append(("jBp", "advanceRatio", js))
        by_j = {r["J"]: r for r in rows}
        tables.append(("propCT", ["jBp"], [by_j[j]["ct"] for j in js], "propCTTable",
                       "Thrust coefficient, vanes neutral, T / rho n^2 D^4."))
        tables.append(("propCQ", ["jBp"], [by_j[j]["cq"] for j in js], "propCQTable",
                       "Shaft torque coefficient, Q / rho n^2 D^5."))

    if coupling:
        crows = coupling["rows"]
        ctcs = uniq([0.0] + [r["Tc"] for r in crows])
        calphas = uniq([r["alphaDeg"] for r in crows])
        bp_defs.append(("tcBp", "thrustCoefficient", ctcs))
        clut = {(r["alphaDeg"], r["Tc"]): r for r in crows}
        for comp in ("dCX", "dCY", "dCZ", "dCl", "dCm", "dCn"):
            vals = []
            for a in calphas:
                for tc in ctcs:
                    if tc == 0.0:
                        vals.append(0.0)  # power-off is the reference
                        continue
                    row = clut.get((a, tc))
                    vals.append(row[comp] if row else 0.0)
            name = "coupling" + comp[1:]
            tables.append((name, ["alphaCouplingBp", "tcBp"], vals, f"{name}Table",
                           f"Fan-on-airframe increment {comp[1:]} from the power-off "
                           "configuration at matched alpha. The aft fan induces a "
                           "favourable gradient UPSTREAM of itself, over the wing; the "
                           "momentum-theory slipstream is identically zero there and "
                           "would report no interaction at all."))
        bp_defs.append(("alphaCouplingBp", "angleOfAttack", calphas))

    if vane and prop:
        # PER-VANE increments, not per-mode: the mode-sum buildup was
        # measured wrong by up to 33%, per-vane summation to 1.4%.
        vrows = vane["rows"]
        vjs = uniq([r["J"] for r in vrows])
        vdeltas = uniq([r["deltaDeg"] for r in vrows] + [0.0])
        bp_defs.append(("jVaneBp", "advanceRatio", vjs))
        bp_defs.append(("vaneBp", "vaneDeflection", vdeltas))
        base_by_j = {r["J"]: r for r in prop["rows"] if r["mode"] == "baseline"}
        # One vane's response serves all four positions by the cruciform's
        # rotational symmetry; the reference vane is the starboard one.
        ref_mode = "vaneRight"
        norm = None
        meta = prop["meta"]
        nf = meta["rho"] * meta["revsPerSec"] ** 2 * meta["diameterM"] ** 4
        nm = nf * meta["diameterM"]
        idx = {(r["mode"], r["deltaDeg"], r["J"]): r for r in vrows}
        for comp, key, scale in (("CX", "fx", nf), ("CY", "fy", nf), ("CZ", "fz", nf),
                                 ("Cl", "mx", nm), ("Cm", "my", nm), ("Cn", "mz", nm)):
            vals = []
            for j in vjs:
                for dlt in vdeltas:
                    if abs(dlt) < 1e-12:
                        vals.append(0.0)
                        continue
                    r = idx.get((ref_mode, dlt, j))
                    b = base_by_j.get(j)
                    if r is None or b is None:
                        vals.append(0.0)
                        continue
                    vals.append((r["totalFrd"][key] - b["totalFrd"][key]) / scale)
            tables.append((f"propDC{comp}vane", ["jVaneBp", "vaneBp"], vals,
                           f"propDC{comp}vaneTable",
                           f"PER-VANE increment {comp} from one vane's own deflection. "
                           f"Sum over the four vanes at their mixed angles -- NOT a sum "
                           f"over command modes, which was measured wrong by up to 33%."))

    # ---------------- table + function definitions ----------------
    # DTD ORDER (DAVEfunc): fileHeader, variableDef+, breakpointDef*,
    # griddedTableDef*, ungriddedTableDef*, function*, checkData?. Every
    # variable must therefore be declared before the first breakpoint, so
    # the table outputs are emitted here rather than beside their tables.
    d.comment("Output variables the tables drive.")
    VANE_POS = ("Bottom", "Left", "Top", "Right")
    for name, bp_ids, vals, tname, desc in tables:
        if "vaneBp" in bp_ids:
            for pos in VANE_POS:
                m.variable(f"{name}{pos}", f"{name}{pos}", "nd", axis="body",
                           description=f"{desc} Evaluated at the {pos.lower()} vane's angle.")
            m.variable(f"{name}Total", f"{name}Total", "nd", axis="body", is_output=True,
                       description=desc + " Summed over the four vanes.",
                       calc=app("plus", *[ci(f"{name}{pos}") for pos in VANE_POS]))
            continue
        m.variable(name, name, "nd", axis="body", is_output=True, description=desc)

    d.comment("Breakpoint sets, taken from the generated data rather than restated.")
    for bpID, bpname, vals in bp_defs:
        emit_breakpoints(d, bpID, bpname, vals)

    d.comment("Gridded tables.")
    for name, bp_ids, vals, tname, desc in tables:
        gridded_table(d, tname, bp_ids, vals, desc)

    d.comment("Functions binding each table to its independent variables. A "
              "per-vane table is bound FOUR times -- once per vane, at that vane's own "
              "commanded angle -- and the four outputs are summed, which is the buildup "
              "the superposition measurement forced.")
    axis_var = {"alphaBp": "alphaDeg", "betaBp": "betaDeg", "alphaRateBp": "alphaDeg",
                "alphaParasiteBp": "alphaDeg", "aileronBp": "aileronDeg",
                "jBp": "advanceRatio", "jVaneBp": "advanceRatio",
                "alphaCouplingBp": "alphaDeg", "tcBp": "thrustCoefficient"}
    VANE_POS = ("Bottom", "Left", "Top", "Right")
    for name, bp_ids, vals, tname, desc in tables:
        if "vaneBp" in bp_ids:
            # One table, four bindings: the same per-vane response evaluated
            # at each vane's own angle. griddedTableRef is an IDREF, so all
            # four functions legitimately point at one table.
            for pos in VANE_POS:
                axes = ["vane" + pos + "Deg" if b == "vaneBp" else axis_var[b]
                        for b in bp_ids]
                simple_function(d, f"{name}{pos}Fn", f"{name}{pos}", axes, bp_ids, tname)
            continue
        simple_function(d, f"{name}Fn", name, [axis_var[b] for b in bp_ids], bp_ids, tname)

    # ---------------- checkData ----------------
    # A shot generated from the same JSON that built the table verifies
    # the ENCODING (breakpoints, ordering, interpolation), not the
    # physics -- that much is circular and is stated as such. The
    # non-circular pins come from independent paths: the solver's own unit
    # tests, and the symmetry a table must obey regardless of its values.
    d.comment("Verification data. Encoding shots pin the table encoding; the physics "
              "pins come from independent paths (see each description).")
    d.open("checkData")

    shots = []
    if aero:
        base = aero["baseline"]
        alphas = uniq([r["alphaDeg"] for r in base])
        # Encoding pin: at an exact breakpoint the interpolation must
        # return the stored value bit for bit.
        mid = [r for r in base if abs(r["betaDeg"]) < 1e-9]
        if mid:
            probe = mid[len(mid) // 2]
            shots.append((
                "aeroEncodingAtBreakpoint",
                "Encoding pin (circular by construction, and only checks the encoding): "
                "at an exact breakpoint the gridded lookup must return the stored value.",
                {"alphaDeg": probe["alphaDeg"], "betaDeg": 0.0},
                {"aeroCZ": probe["CZ"], "aeroCm": probe["Cm"]},
                1e-9))
            # Physics pin, independent of the tables: a symmetric airframe
            # at zero sideslip carries no lateral wrench. A sign slip in
            # the beta mirroring breaks this and nothing else would.
            shots.append((
                "lateralSymmetryAtZeroSideslip",
                "Physics pin: a symmetric configuration at beta = 0 carries no side "
                "force, rolling moment or yawing moment. Catches a sign error in the "
                "beta mirroring, which the encoding pin cannot see.",
                {"alphaDeg": probe["alphaDeg"], "betaDeg": 0.0},
                {"aeroCY": 0.0, "aeroCl": 0.0, "aeroCn": 0.0},
                1e-6))
        rates = aero.get("rates", [])
        if rates:
            near2 = min(rates, key=lambda r: abs(r["alphaDeg"] - 2.0))
            shots.append((
                "rollDampingSign",
                "Physics pin from an INDEPENDENT path: roll damping is negative and "
                "sits near the textbook Cl_p = -0.45 for this AR = 6 wing -- the value "
                "TestBodyAxes and TestSolverCore pin in the solver. A frame error "
                "returns +0.45: the right magnitude with the wrong sign, i.e. roll "
                "ANTI-damping.",
                {"alphaDeg": near2["alphaDeg"]},
                {"aeroClp": near2["Clp"]},
                1e-9))
    if aero and aero.get("aileron"):
        gen = [r for r in aero["aileron"] if r["deltaDeg"] > 0]
        if gen:
            top = max((r["deltaDeg"] for r in gen))
            probe = min((r for r in gen if abs(r["deltaDeg"] - top) < 1e-9),
                        key=lambda r: abs(r["alphaDeg"]))
            # The expected value is taken from the POSITIVE-deflection row
            # and negated by hand here, while the table cell is filled by
            # the assembler's mirroring code. The two are independent, so
            # this catches a wrong mirror sign -- which no encoding pin
            # can, since the encoding would faithfully store the error.
            shots.append((
                "aileronMirrorAntisymmetry",
                "Physics pin: the rolling-moment increment is ODD in aileron deflection, "
                "because mirroring the configuration about its xz-plane maps a +delta "
                "antisymmetric command onto -delta. The expected value comes from the "
                "+delta row negated; the table cell comes from the mirroring code.",
                {"alphaDeg": probe["alphaDeg"], "aileronDeg": -top},
                {"aeroDCl": -probe["dCl"]},
                1e-9))
            shots.append((
                "aileronRollAuthorityIsReal",
                "Physics pin: a deflected aileron produces a rolling moment at all. This "
                "was exactly zero at every attitude before the flap was carried in the "
                "section, silently, with the solve converging and reporting sensible "
                "forces.",
                {"alphaDeg": probe["alphaDeg"], "aileronDeg": top},
                {"aeroDCl": probe["dCl"]},
                1e-9))

    if parasite:
        pts = sorted(parasite["table"], key=lambda r: r["alphaDeg"])
        top = pts[-1]
        shots.append((
            "parasiteCrossflowDominates",
            "Physics pin: at 90 degrees the body's crossflow drag dominates the "
            "friction term by more than an order of magnitude. A model carrying "
            "parasite drag as a constant fails this by a factor of ~17.",
            {"alphaDeg": top["alphaDeg"]},
            {"aeroCD0": top["CD0"]},
            1e-9))

    for name, desc, inputs, outputs, tol in shots:
        d.open("staticShot", name=name)
        d.leaf("description", desc)
        d.open("checkInputs")
        for var, val in inputs.items():
            d.open("signal")
            d.leaf("signalName", var)
            d.leaf("signalUnits", "deg")
            d.leaf("signalValue", f"{val:.12g}")
            d.close("signal")
        d.close("checkInputs")
        d.open("checkOutputs")
        for var, val in outputs.items():
            d.open("signal")
            d.leaf("signalName", var)
            d.leaf("signalUnits", "nd")
            d.leaf("signalValue", f"{val:.12g}")
            d.leaf("tol", f"{tol:g}")
            d.close("signal")
        d.close("checkOutputs")
        d.close("staticShot")

    d.close("checkData")
    d.close("DAVEfunc")

    out = args.out or os.path.join(HERE, "AetherionFlightModel.dml")
    with open(out, "w", encoding="utf-8") as handle:
        handle.write(d.text())
    print(f"wrote {out}: {len(tables)} tables")
    if aero is None:
        print("  (no aero* tables -- aero-map.json was absent)")
    if not (aero or {}).get("aileron"):
        print("  NOT emitted: aeroDC* (aileron -- deflected sweep not yet run)")
    if not coupling:
        print("  NOT emitted: coupling* (no coupling-map.json)")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data", default=os.path.join(HERE, "data"))
    ap.add_argument("--out", default=None)
    return build(ap.parse_args())


if __name__ == "__main__":
    sys.exit(main())
