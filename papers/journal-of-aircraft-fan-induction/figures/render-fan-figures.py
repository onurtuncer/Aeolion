# render-fan-figures.py -- figures and tables for the aft-fan induction
# paper.
#
# Input:  fan-induction.json  -- written by `aeolion_fan_induction`
#         fine-off.json / fine-on.json -- 0.5-degree alpha sweeps written by
#         `aeolion_attachment_sweep` with thrust 0 and 25 N at V = 12 m/s
#         (see README.md for the exact commands).
# Output: axis-induction.pdf   -- the closed-form check, and the upstream half
#         chordwise-gradient.pdf -- the mechanism, across the span
#         separation-shift.pdf -- the deliverable
#         tables/*.tex         -- \input-ed by paper.tex
#
# Encoding follows the data's job. The operating points are an ordered
# progression through transition, so they take a single-hue sequential
# ramp; fan on/off is a two-state comparison and takes two clearly
# separated hues rather than a ramp.

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).parent
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
# Sequential: the operating points are an ordered walk through transition.
RAMP = ["#9ecae1", "#4292c6", "#2171b5", "#08306b"]
# Two-state comparison, clearly separated.
OFF, ON = "#737373", "#2166ac"


def style(ax):
    ax.grid(True, **GRID)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(colors=INK, labelsize=7.5)
    for s in ax.spines.values():
        s.set_color(INK)


def load():
    d = json.loads((HERE / "fan-induction.json").read_text())
    off = json.loads((HERE / "fine-off.json").read_text())
    on = json.loads((HERE / "fine-on.json").read_text())
    return d, off, on



# ------------------------------- fig 0: the configuration being solved -----
def fig_configuration(d):
    """Planform and meridional views of the coupled wing-body-duct solve,
    with the induction field the fan puts on it. This figure exists to SHOW
    the paper's two geometric claims rather than assert them: the duct sits
    aft of the wing, and its radius covers the span station that separates
    first."""
    cfg = d["configuration"]
    m = d["meta"]
    fig, (top, bot) = plt.subplots(2, 1, figsize=(5.0, 4.4), sharex=True)

    # --- (a) planform: the lattice actually solved --------------------------
    wing = np.array(cfg["wing"])          # y, xLE, xTE per strip
    body = np.array(cfg["body"])          # x, r
    duct = cfg["duct"]
    semi = m["semiSpan"]

    for y, xle, xte in wing:
        top.plot([xle, xte], [y, y], color="#4292c6", linewidth=0.7, solid_capstyle="butt")
    top.plot(body[:, 0], body[:, 1], color=INK, linewidth=0.9)
    top.plot(body[:, 0], -body[:, 1], color=INK, linewidth=0.9)

    dw = duct["x1"] - duct["x0"]
    for sgn in (+1, -1):
        top.add_patch(plt.Rectangle((duct["x0"], sgn * duct["rInner"]), dw,
                                    sgn * (duct["rOuter"] - duct["rInner"]),
                                    facecolor="0.7", edgecolor=INK, linewidth=0.7,
                                    zorder=4))
    # the fan disk, edge on
    top.plot([m["diskX"], m["diskX"]], [-m["radius"], -m["hubRadius"]], color="#b2182b",
             linewidth=1.6)
    top.plot([m["diskX"], m["diskX"]], [m["hubRadius"], m["radius"]], color="#b2182b",
             linewidth=1.6, label="fan disk")

    # The alignment claim, drawn.
    for sgn in (+1, -1):
        top.axhline(sgn * m["radius"], color="#b2182b", linewidth=0.5, linestyle=":")
    top.axhspan(0.15 * semi, 0.23 * semi, color="#fddbc7", zorder=0)
    top.axhspan(-0.23 * semi, -0.15 * semi, color="#fddbc7", zorder=0)
    label = "first to separate" + chr(10) + r"($|2y/b| = 0.15$--$0.23$)"
    top.annotate(label, xy=(0.615, 0.19 * semi), fontsize=6.4, color="0.45",
                 va="center", ha="left")
    top.annotate("wing", xy=(0.29, 0.42), fontsize=7, color="#2171b5", ha="center")
    top.annotate("duct", xy=(duct["x0"] - 0.005, 0.155), fontsize=7, color=INK, ha="right")

    top.set_ylabel(r"$y$ [m]", color=INK)
    top.set_aspect("equal")
    top.set_xlim(0.0, 0.60)
    top.set_ylim(-0.58, 0.58)
    top.legend(loc="lower right", frameon=False, fontsize=6.8, labelcolor=INK)
    style(top)

    # --- (b) meridional: the induction field on the same geometry ----------
    f = np.array(d["field"]["points"])
    xs = np.unique(f[:, 0]); ys = np.unique(f[:, 1])
    U = f[:, 2].reshape(len(xs), len(ys)).T / d["field"]["vi"]
    mesh = bot.pcolormesh(xs, ys, U, cmap="RdBu_r",
                          norm=matplotlib.colors.TwoSlopeNorm(vmin=-0.15, vcenter=0.0, vmax=0.6),
                          shading="gouraud", zorder=0)
    bot.plot(body[:, 0], body[:, 1], color=INK, linewidth=0.9)
    bot.plot(body[:, 0], -body[:, 1], color=INK, linewidth=0.9)
    for y, xle, xte in wing[::3]:
        bot.plot([xle, xte], [y, y], color="0.15", linewidth=0.5)
    bot.plot([m["diskX"]] * 2, [-m["radius"], -m["hubRadius"]], color="#111111", linewidth=1.6)
    bot.plot([m["diskX"]] * 2, [m["hubRadius"], m["radius"]], color="#111111", linewidth=1.6)

    bot.set_xlabel(r"$x$ [m] (aft)", color=INK)
    bot.set_ylabel(r"$y$ [m]", color=INK)
    bot.set_aspect("equal")
    bot.set_ylim(-0.58, 0.58)
    style(bot)
    bar = fig.colorbar(mesh, ax=(top, bot), fraction=0.025, pad=0.015)
    bar.set_label(r"$u_x / v_i$", color=INK)
    bar.ax.tick_params(colors=INK, labelsize=6.5)
    bar.outline.set_visible(False)

    fig.subplots_adjust(left=0.10, right=0.86, bottom=0.09, top=0.98, hspace=0.02)
    return fig


# ------------------------------------------------ fig 1: the axis check ----
def fig_axis(d):
    a = np.array(d["axis"])
    fig, ax = plt.subplots(figsize=(3.5, 2.5))

    ax.plot(a[:, 0], a[:, 1], color=ON, linewidth=1.4, label="solid disk, closed form")
    ax.plot(a[::8, 0], a[::8, 2], color=ON, linestyle="none", marker="o", markersize=2.4,
            markerfacecolor="white", markeredgewidth=0.7, label="discretized sheet")
    ax.plot(a[:, 0], a[:, 3], color=OFF, linewidth=1.2, linestyle="--",
            label="annulus (with bore)")

    ax.axvline(0.0, color="0.55", linewidth=0.7, linestyle=":")
    ax.annotate("disk plane", xy=(0.0, 1.85), xytext=(3, 0), textcoords="offset points",
                fontsize=6.8, color="0.4")
    for y, t in ((1.0, r"$v_i$"), (2.0, r"$2v_i$")):
        ax.axhline(y, color="0.85", linewidth=0.6)
        ax.annotate(t, xy=(-5.9, y), xytext=(0, 2), textcoords="offset points",
                    fontsize=6.8, color="0.4")

    # The number the interaction argument turns on.
    ax.plot([-1.0], [1.0 - 1.0 / np.sqrt(2.0)], marker="o", markersize=3.5, color=INK,
            zorder=5)
    ax.annotate(r"$0.293\,v_i$ at $s=-R$", xy=(-1.0, 1.0 - 1.0 / np.sqrt(2.0)),
                xytext=(-46, 26), textcoords="offset points", fontsize=7, color=INK,
                arrowprops=dict(arrowstyle="->", color="0.5", linewidth=0.6))

    ax.set_xlabel(r"$s/R$   (negative is upstream)", color=INK)
    ax.set_ylabel(r"$u_x / v_i$", color=INK)
    ax.set_xlim(-6, 10)
    ax.set_ylim(-0.35, 2.15)
    ax.legend(loc="center right", frameon=False, fontsize=7, labelcolor=INK)
    style(ax)
    fig.subplots_adjust(left=0.14, right=0.97, bottom=0.18, top=0.96)
    return fig


# ------------------------------- fig 2: the chordwise gradient, spanwise ----
def fig_gradient(d):
    fig, ax = plt.subplots(figsize=(3.5, 2.6))
    for i, op in enumerate(d["operating"]):
        s = np.array(op["span"])
        ax.plot(s[:, 0], 100.0 * s[:, 3], color=RAMP[i], linewidth=1.3,
                label=rf"$\mu={op['mu']:.2f}$  ($V_\infty={op['Vinf']:.0f}$)")
    ax.axhline(0.0, color="0.5", linewidth=0.8)

    # The sign change is the non-obvious result; mark where it happens.
    s = np.array(d["operating"][1]["span"])
    k = int(np.argmin(np.abs(s[:, 3])))
    ax.annotate("fan penalizes\nthe outer wing", xy=(s[k, 0], 0.0), xytext=(6, -26),
                textcoords="offset points", fontsize=7, color=INK,
                arrowprops=dict(arrowstyle="->", color="0.5", linewidth=0.6))

    ax.set_xlabel(r"$2y/b$", color=INK)
    ax.set_ylabel(r"chordwise rise $\Delta u_x / V_\infty$ [\%]", color=INK)
    ax.set_xlim(0.0, 1.0)
    ax.legend(loc="upper right", frameon=False, fontsize=7, labelcolor=INK)
    style(ax)
    fig.subplots_adjust(left=0.16, right=0.97, bottom=0.17, top=0.96)
    return fig


# --------------------------------------- fig 3: the separation-curve shift ---
def curve(sweep):
    return {c["alphaDeg"]: min((s["psi"] for s in c["separation"]["stations"]
                                if s["mode"] == "turbulent"), default=1.0)
            for c in sweep["conditions"]}


def alpha_at(c, level):
    a = sorted(c)
    for i in range(len(a) - 1):
        y0, y1 = c[a[i]], c[a[i + 1]]
        if (y0 - level) * (y1 - level) <= 0 and y0 != y1:
            return a[i] + (a[i + 1] - a[i]) * (y0 - level) / (y0 - y1)
    return None


LEVELS = (0.90, 0.85, 0.80, 0.75, 0.70, 0.60, 0.55)


def fig_shift(off, on):
    U, P = curve(off), curve(on)
    fig, ax = plt.subplots(figsize=(3.5, 2.6))
    for C, col, lab in ((U, OFF, "fan off"), (P, ON, "fan on, $T=25$ N")):
        a = np.array(sorted(C))
        ax.plot(a, [C[x] for x in a], color=col, linewidth=1.4, label=lab)

    # The shift itself, drawn at one criterion so the reader sees what is
    # being measured rather than only the table of numbers.
    lvl = 0.80
    ao, an = alpha_at(U, lvl), alpha_at(P, lvl)
    ax.plot([ao, an], [lvl, lvl], color=INK, linewidth=0.9,
            marker="|", markersize=6)
    ax.annotate(rf"$\Delta\alpha = {an - ao:+.2f}^\circ$",
                xy=(0.5 * (ao + an), lvl), xytext=(0, 8), textcoords="offset points",
                ha="center", fontsize=7.5, color=INK)
    ax.axhline(lvl, color="0.85", linewidth=0.6, zorder=0)

    ax.set_xlabel(r"$\alpha$ [deg]", color=INK)
    ax.set_ylabel(r"separation station $x_{\mathrm{sep}}/c_n$", color=INK)
    ax.set_ylim(0.45, 1.02)
    ax.legend(loc="lower left", frameon=False, fontsize=7, labelcolor=INK)
    style(ax)
    fig.subplots_adjust(left=0.17, right=0.97, bottom=0.17, top=0.96)
    return fig


# ----------------------------------------------------------------- tables ---
def table_shift(off, on):
    U, P = curve(off), curve(on)
    rows, ds = [], []
    for lvl in LEVELS:
        ao, an = alpha_at(U, lvl), alpha_at(P, lvl)
        if ao is None or an is None:
            continue
        ds.append(an - ao)
        rows.append(rf"  {lvl:.2f} & {ao:.2f} & {an:.2f} & {an - ao:+.2f} \\")
    mean = np.mean(ds[1:])  # exclude the near-onset level, reported separately
    return "\n".join([
        r"\begin{tabular}{rrrr}", r"  \hline",
        r"  $x_{\mathrm{sep}}/c_n$ & $\alpha$ off [deg] & $\alpha$ on [deg]"
        r" & $\Delta\alpha$ [deg] \\",
        r"  \hline", *rows, r"  \hline",
        rf"  \multicolumn{{3}}{{r}}{{mean, $0.85$ to $0.55$}} & {mean:+.2f} \\",
        r"  \hline", r"\end{tabular}",
    ])


def table_operating(d):
    rows = []
    for op in d["operating"]:
        s = np.array(op["span"])
        def at(eta):
            return 100.0 * s[int(np.argmin(np.abs(s[:, 0] - eta))), 3]
        rows.append(rf"  {op['Vinf']:.0f} & {op['thrust']:.0f} & {op['vi']:.1f} & "
                    rf"{op['mu']:.2f} & {at(0.10):.1f} & {at(0.20):.1f} & {at(0.30):.1f} & "
                    rf"{at(0.60):.1f} \\")
    return "\n".join([
        r"\begin{tabular}{rrrrrrrr}", r"  \hline",
        r"  $V_\infty$ & $T$ & $v_i$ & $\mu$ & \multicolumn{4}{c}{"
        r"$\Delta u_x/V_\infty$ [\%] at $2y/b =$} \\",
        # A row starting with "[" right after \ is read as \[<length>];
        # the empty group stops LaTeX looking for an optional argument.
        r"  {}[m/s] & [N] & [m/s] & & 0.10 & 0.20 & 0.30 & 0.60 \\",
        r"  \hline", *rows, r"  \hline", r"\end{tabular}",
    ])


def main():
    d, off, on = load()
    TABLES.mkdir(exist_ok=True)
    for name, fig in (("configuration", fig_configuration(d)),
                      ("axis-induction", fig_axis(d)),
                      ("chordwise-gradient", fig_gradient(d)),
                      ("separation-shift", fig_shift(off, on))):
        fig.savefig(HERE / f"{name}.pdf")
        fig.savefig(HERE / f"{name}.png", dpi=220)
        plt.close(fig)
        print(f"wrote {name}.pdf / .png")
    (TABLES / "shift.tex").write_text(table_shift(off, on) + "\n")
    (TABLES / "operating.tex").write_text(table_operating(d) + "\n")
    print("wrote tables/shift.tex, tables/operating.tex")


if __name__ == "__main__":
    main()
