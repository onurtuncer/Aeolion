# render-coefficient-tables.py -- the Journal of Aircraft paper's
# coefficient and derivative tables, plus the separation-boundary figure.
#
# Input:  attachment-sweep.json, written by `aeolion_attachment_sweep`
#         (app/AttachmentSweepExport.cpp).
# Output: tables/coefficients.tex   -- \input-ed by paper.tex
#         tables/derivatives.tex
#         tables/separation.tex
#         separation-boundary.pdf/.png
#
# The tables are GENERATED rather than typed, for the same reason the
# figures are rendered from the solver's own export: a hand-copied table is
# a second source of truth that silently goes stale. Every number in the
# paper's Section VI tables traces to one solve.
#
# Attitudes at or beyond the separation boundary are reported but MARKED,
# never silently dropped: the reader is told which rows the potential-flow
# coefficients are trustworthy for, which is the whole point of running the
# separation march alongside the sweep.

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).parent
DATA = HERE / "attachment-sweep.json"
TABLES = HERE / "tables"

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 8.0,
    "pdf.fonttype": 42,
    "axes.linewidth": 0.6,
})

INK = "0.25"
GRID = dict(color="0.85", linewidth=0.5)
# Diverging: beta is a polarity about zero (see render-attachment-figures.py).
BETA_COLORS = {-10.0: "#b2182b", -5.0: "#d6604d", 0.0: "#737373",
               5.0: "#4393c3", 10.0: "#2166ac"}
BETA_MARKERS = {-10.0: "<", -5.0: "v", 0.0: "o", 5.0: "^", 10.0: ">"}

# A station counts as "separated" for the boundary when its upper-surface
# turbulent separation has moved forward of the trailing edge at all. The
# ENGINEERING boundary is stricter -- separation ahead of this chord
# station -- because a few percent of trailing-edge separation is normal
# and not a limit.
SEPARATION_LIMIT_PSI = 0.90


def load():
    return json.loads(DATA.read_text())


def by_attitude(data):
    return {(c["alphaDeg"], c["betaDeg"]): c for c in data["conditions"]}


def alphas(data):
    return sorted({c["alphaDeg"] for c in data["conditions"]})


def betas(data):
    return sorted({c["betaDeg"] for c in data["conditions"]})


def separation_extent(cond):
    """Fraction of resolved stations whose upper surface separates forward of
    SEPARATION_LIMIT_PSI, and the most forward separation station."""
    stations = cond["separation"]["stations"]
    if not stations:
        return 0.0, 1.0
    hits = [s for s in stations
            if s["mode"] == "turbulent" and s["psi"] < SEPARATION_LIMIT_PSI]
    forwardmost = min((s["psi"] for s in hits), default=1.0)
    return len(hits) / len(stations), forwardmost


def attached(cond):
    fraction, _ = separation_extent(cond)
    return fraction == 0.0


# ------------------------------------------------------------- tables ----
def fmt(value, places=4):
    return f"{value:.{places}f}"


def coefficient_table(data):
    """CL, CDi, Cm across the matrix; separated attitudes daggered."""
    rows = []
    for a in alphas(data):
        cells = []
        for b in betas(data):
            cond = by_attitude(data)[(a, b)]
            mark = "" if attached(cond) else r"$^\dagger$"
            cells.append(f"{fmt(cond['CL'], 3)}{mark}")
        rows.append(rf"  {a:5.0f} & " + " & ".join(cells) + r" \\")

    header = " & ".join(rf"$\beta={b:.0f}^\circ$" for b in betas(data))
    return "\n".join([
        r"\begin{tabular}{r" + "r" * len(betas(data)) + r"}",
        r"  \hline",
        rf"  $\alpha$ & {header} \\",
        r"  \hline",
        *rows,
        r"  \hline",
        r"\end{tabular}",
    ])


def full_coefficient_table(data):
    """Every coefficient at beta = 0 and beta = 10, attached rows only."""
    rows = []
    for a in alphas(data):
        for b in (0.0, 10.0):
            cond = by_attitude(data)[(a, b)]
            if not attached(cond):
                continue
            rows.append(
                rf"  {a:5.0f} & {b:5.0f} & {fmt(cond['CL'], 4)} & {fmt(cond['CDi'], 5)} & "
                rf"{fmt(cond['CY'], 5)} & {fmt(cond['Cm'], 4)} & {fmt(cond['Croll'], 5)} & "
                rf"{fmt(cond['Cn'], 5)} \\")
    return "\n".join([
        r"\begin{tabular}{rrrrrrrr}",
        r"  \hline",
        r"  $\alpha$ & $\beta$ & $C_L$ & $C_{D_i}$ & $C_Y$ & $C_m$ & $C_l$ & $C_n$ \\",
        r"  \hline",
        *rows,
        r"  \hline",
        r"\end{tabular}",
    ])


def derivative_table(data):
    """The derivative set at beta = 0, attached rows only."""
    rows = []
    for a in alphas(data):
        cond = by_attitude(data)[(a, 0.0)]
        if not attached(cond):
            continue
        d = cond["derivatives"]
        rows.append(
            rf"  {a:5.0f} & {fmt(d['CL_alpha'], 3)} & {fmt(d['Cm_alpha'], 3)} & "
            rf"{fmt(d['CY_beta'], 4)} & {fmt(d['Croll_beta'], 4)} & {fmt(d['Cn_beta'], 4)} & "
            rf"{fmt(d['Cm_q_nd'], 3)} & {fmt(d['Croll_p_nd'], 4)} & {fmt(d['Cn_r_nd'], 4)} \\")
    return "\n".join([
        r"\begin{tabular}{rrrrrrrrr}",
        r"  \hline",
        r"  $\alpha$ & $C_{L_\alpha}$ & $C_{m_\alpha}$ & $C_{Y_\beta}$ & $C_{l_\beta}$ & "
        r"$C_{n_\beta}$ & $C_{m_q}$ & $C_{l_p}$ & $C_{n_r}$ \\",
        r"  \hline",
        *rows,
        r"  \hline",
        r"\end{tabular}",
    ])


def separation_table(data):
    """Where the upper surface separates, across the matrix."""
    rows = []
    for a in alphas(data):
        cells = []
        for b in betas(data):
            fraction, forwardmost = separation_extent(by_attitude(data)[(a, b)])
            cells.append("---" if fraction == 0.0
                         else rf"{forwardmost:.2f} ({100 * fraction:.0f}\%)")
        rows.append(rf"  {a:5.0f} & " + " & ".join(cells) + r" \\")
    header = " & ".join(rf"$\beta={b:.0f}^\circ$" for b in betas(data))
    return "\n".join([
        r"\begin{tabular}{r" + "r" * len(betas(data)) + r"}",
        r"  \hline",
        rf"  $\alpha$ & {header} \\",
        r"  \hline",
        *rows,
        r"  \hline",
        r"\end{tabular}",
    ])


# ------------------------------------------------- separation boundary ----
# Sequential ramp: alpha is an ordered magnitude (see the palette note in
# render-attachment-figures.py). Validated with the dataviz validator.
ALPHA_RAMP = ["#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b"]


def fig_separation_boundary(data):
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(3.5, 4.4))

    # (a) WHERE separation sits across the span, as alpha grows. This is the
    # mechanism -- separation is most advanced just outboard of the wing
    # root and recovers toward the tips -- and it is what five collapsed
    # beta curves cannot show.
    shown = [a for a in alphas(data) if a >= 4.0]
    ramp = {a: ALPHA_RAMP[min(i, len(ALPHA_RAMP) - 1)] for i, a in enumerate(shown)}
    for a in shown:
        stations = sorted(by_attitude(data)[(a, 0.0)]["separation"]["stations"],
                          key=lambda s: s["eta"])
        eta = np.array([s["eta"] for s in stations])
        psi = np.array([s["psi"] if s["mode"] == "turbulent" else np.nan for s in stations])
        if np.all(np.isnan(psi)):
            continue
        k = int(np.searchsorted(eta, 0.0))  # break across the fuselage gap
        top.plot(np.insert(eta, k, np.nan), np.insert(psi, k, np.nan),
                 color=ramp[a], linewidth=1.2)
        # Label at the RIGHT tip, where the curves fan apart; the crowded
        # low-alpha end is left to the ramp and the two labelled extremes.
        finite = np.where(~np.isnan(psi))[0]
        if a in (shown[0], shown[len(shown) // 2], shown[-1]):
            top.annotate(rf"$\alpha={a:.0f}^\circ$", xy=(eta[finite[-1]], psi[finite[-1]]),
                         xytext=(4, 0), textcoords="offset points", fontsize=7,
                         color=ramp[a], va="center")

    top.set_xlabel(r"$2y/b$", color=INK)
    top.set_ylabel(r"separation $x/c_n$", color=INK)
    top.set_xlim(-1.05, 1.32)
    top.set_ylim(0.5, 1.02)

    # (b) The boundary itself. The five beta curves collapse -- that IS the
    # result, so they are drawn together and the caption says so.
    for b in betas(data):
        xs = list(alphas(data))
        extent = [100.0 * separation_extent(by_attitude(data)[(a, b)])[0] for a in xs]
        bottom.plot(xs, extent, color=BETA_COLORS[b], marker=BETA_MARKERS[b],
                    markersize=2.8, linewidth=1.0,
                    linestyle="--" if b == 0.0 else "-",
                    label=rf"$\beta = {b:+.0f}^\circ$")

    onset = next((a for a in alphas(data)
                  if separation_extent(by_attitude(data)[(a, 0.0)])[0] > 0.0), None)
    if onset is not None:
        bottom.axvline(onset, color="0.55", linewidth=0.7, linestyle=":")
        bottom.annotate(rf"onset $\alpha = {onset:.0f}^\circ$", xy=(onset, 60),
                        xytext=(4, 0), textcoords="offset points", fontsize=7, color="0.4")

    bottom.set_ylabel("span separated [\\%]", color=INK)
    bottom.set_xlabel(r"$\alpha$ [deg]", color=INK)
    bottom.set_ylim(-3, 103)
    bottom.legend(loc="upper left", frameon=False, fontsize=6.8, handletextpad=0.4,
                  labelcolor=INK)

    for ax in (top, bottom):
        ax.grid(True, **GRID)
        ax.set_axisbelow(True)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        ax.tick_params(colors=INK, labelsize=7.5)
        for spine in ax.spines.values():
            spine.set_color(INK)

    fig.subplots_adjust(left=0.17, right=0.97, bottom=0.09, top=0.97, hspace=0.38)
    return fig


def main():
    data = load()
    TABLES.mkdir(exist_ok=True)

    (TABLES / "coefficients.tex").write_text(full_coefficient_table(data) + "\n")
    (TABLES / "cl-matrix.tex").write_text(coefficient_table(data) + "\n")
    (TABLES / "derivatives.tex").write_text(derivative_table(data) + "\n")
    (TABLES / "separation.tex").write_text(separation_table(data) + "\n")
    print(f"wrote {TABLES.name}/: coefficients, cl-matrix, derivatives, separation")

    fig = fig_separation_boundary(data)
    fig.savefig(HERE / "separation-boundary.pdf")
    fig.savefig(HERE / "separation-boundary.png", dpi=220)
    plt.close(fig)
    print("wrote separation-boundary.pdf / .png")

    # A short console summary, so a rerun says what moved without opening
    # the paper.
    print("\nseparation onset (first alpha with any station separated forward of "
          f"{SEPARATION_LIMIT_PSI}):")
    for b in betas(data):
        onset = next((a for a in alphas(data)
                      if separation_extent(by_attitude(data)[(a, b)])[0] > 0.0), None)
        print(f"  beta = {b:+5.1f} deg: "
              + (f"alpha = {onset:+.0f} deg" if onset is not None else "attached throughout"))


if __name__ == "__main__":
    main()
