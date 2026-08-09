# render-attachment-figures.py -- publication figures for the Journal of
# Aircraft paper's application section: the alpha/beta sweep of the coupled
# wing--body solve on the real geometry handoff.
#
# Input:  attachment-sweep.json, written by `aeolion_attachment_sweep`
#         (app/AttachmentSweepExport.cpp) -- rebuild and rerun that tool
#         whenever the solver or the handoff fixture changes, then rerun
#         this script.
# Output: attachment-offset.pdf    -- stagnation offset vs. span, alpha sweep
#         sideslip-asymmetry.pdf   -- effective sweep + Rbar vs. span, beta sweep
#         body-topology.pdf        -- nose attachment node migration
#         body-streamlines.pdf     -- fuselage surface streamlines colored by Cp
#         (+ .png previews of each)
#
# Encoding follows the data's job: the alpha levels are ordered magnitude,
# so they take a single-hue sequential ramp; beta is a polarity about zero,
# so it takes a diverging warm/cool pair with a neutral dashed zero; Cp is
# a signed field about the freestream reference, so it takes the diverging
# colormap centered at zero. Every curve carries its own marker and a
# direct label, so identity survives grayscale print.

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection

HERE = Path(__file__).parent
DATA = HERE / "attachment-sweep.json"

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 8.0,
    "pdf.fonttype": 42,
    "axes.linewidth": 0.6,
})

INK = "0.25"  # annotation ink: recessive gray, never pure black
GRID = dict(color="0.85", linewidth=0.5)
ALPHA_COLORS = {0.0: "#6baed6", 4.0: "#2171b5", 8.0: "#08306b"}  # sequential: alpha is a magnitude
ALPHA_MARKERS = {0.0: "o", 4.0: "s", 8.0: "^"}
BETA_COLORS = {-10.0: "#b2182b", 0.0: "#737373", 10.0: "#2166ac"}  # diverging: beta is a polarity
BETA_MARKERS = {-10.0: "<", 0.0: "o", 10.0: ">"}
CMAP = plt.get_cmap("RdBu_r")


def load():
    return json.loads(DATA.read_text())


def condition(data, alpha, beta):
    for c in data["conditions"]:
        if c["alphaDeg"] == alpha and c["betaDeg"] == beta:
            return c
    raise KeyError(f"no condition alpha={alpha} beta={beta}")


def line_arrays(cond, key="attachmentLine"):
    """Found, non-carry stations as (eta, station-dict) sorted by eta."""
    stations = [s for s in cond[key] if s["found"] and not s["carryStrip"]]
    stations.sort(key=lambda s: s["eta"])
    return stations


def gapped(stations, key):
    """(eta, value) arrays with a NaN break across the fuselage gap, so the
    two semi-spans are not drawn joined through a wing that is not there."""
    eta = np.array([s["eta"] for s in stations])
    val = np.array([float(s[key]) for s in stations])
    k = int(np.searchsorted(eta, 0.0))
    return np.insert(eta, k, np.nan), np.insert(val, k, np.nan)


def style_axes(ax):
    ax.grid(True, **GRID)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(colors=INK, labelsize=7.5)
    for spine in ax.spines.values():
        spine.set_color(INK)


# --------------------------------------------------- fig 1: offset vs span
def fig_attachment_offset(data):
    fig, ax = plt.subplots(figsize=(3.5, 2.4))
    for alpha in (0.0, 4.0, 8.0):
        cond = condition(data, alpha, 0.0)
        stations = line_arrays(cond)
        eta, off = gapped(stations, "stagnationOffset")
        ax.plot(eta, off, color=ALPHA_COLORS[alpha], marker=ALPHA_MARKERS[alpha],
                markersize=2.6, linewidth=1.1, markevery=3, clip_on=False)
        ax.annotate(rf"$\alpha = {alpha:.0f}^\circ$", xy=(eta[-1], off[-1]),
                    xytext=(4, 0), textcoords="offset points", va="center",
                    fontsize=7.5, color=ALPHA_COLORS[alpha])
    ax.set_xlabel(r"$2y/b$", color=INK)
    ax.set_ylabel(r"$s_{\mathrm{stag}}/c_n$", color=INK)
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(0.0, 0.062)
    style_axes(ax)
    fig.subplots_adjust(left=0.15, right=0.87, bottom=0.18, top=0.96)
    return fig


# ------------------------------------------- fig 2: sideslip asymmetry
def fig_sideslip(data):
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(3.5, 4.2), sharex=True)
    labels = {-10.0: r"$\beta = -10^\circ$", 0.0: r"$\beta = 0$", 10.0: r"$\beta = +10^\circ$"}
    label_at = {-10.0: 0.36, 0.0: 0.62, 10.0: -0.62}  # eta where each direct label sits

    for beta in (-10.0, 0.0, 10.0):
        cond = condition(data, 8.0, beta)
        stations = line_arrays(cond)
        eta, sweep = gapped(stations, "effectiveSweepDeg")
        _, rbar = gapped(stations, "Rbar")
        style = dict(color=BETA_COLORS[beta], marker=BETA_MARKERS[beta], markersize=2.6,
                     linewidth=1.1, markevery=3,
                     linestyle="--" if beta == 0.0 else "-")
        top.plot(eta, sweep, **style)
        bottom.plot(eta, rbar, **style)
        k = int(np.nanargmin(np.abs(eta - label_at[beta])))
        top.annotate(labels[beta], xy=(eta[k], sweep[k]), xytext=(0, 5),
                     textcoords="offset points", ha="center", fontsize=7.5,
                     color=BETA_COLORS[beta])

    top.set_ylabel(r"$\Lambda_{\mathrm{eff}}$ [deg]", color=INK)
    top.set_ylim(0, 34)

    bottom.set_yscale("log")
    bottom.set_ylim(0.3, 900)
    for value, name in ((245.0, "contamination"), (583.0, "transition")):
        bottom.axhline(value, color="0.55", linewidth=0.7, linestyle=":")
        bottom.annotate(rf"$\bar{{R}} = {value:.0f}$ ({name})", xy=(-1.0, value),
                        xytext=(0, 2), textcoords="offset points", fontsize=6.8, color="0.4")
    bottom.set_ylabel(r"$\bar{R}$", color=INK)
    bottom.set_xlabel(r"$2y/b$", color=INK)
    bottom.set_xlim(-1.05, 1.05)

    for ax in (top, bottom):
        style_axes(ax)
    fig.subplots_adjust(left=0.15, right=0.97, bottom=0.10, top=0.97, hspace=0.12)
    return fig


# --------------------------------------------- fig 3: body node migration
def fig_body_topology(data):
    # Solver frame: x aft (nose at 0), y right, z up. Meridian angle measured
    # from the windward-bottom meridian, positive toward +y (right).
    meta = data["meta"]
    surface = condition(data, meta["figureAlphaDeg"], meta["figureBetaDeg"])["surface"]
    body_length = max(p[0] for p in surface["points"])
    fig, ax = plt.subplots(figsize=(3.5, 2.6))
    for beta in (-10.0, 0.0, 10.0):
        xs, phis, alphas = [], [], []
        for alpha in (0.0, 4.0, 8.0):
            cond = condition(data, alpha, beta)
            nodes = [p for p in cond["body"]["criticalPoints"] if p["type"] == "attachment_node"]
            if cond["body"]["attachesUpstream"]:
                xs.append(0.0)
                phis.append(0.0)
                alphas.append(alpha)
                continue
            if not nodes:
                continue
            x, y, z = nodes[0]["point"]
            xs.append(x / body_length)
            phis.append(np.degrees(np.arctan2(y, -z)))
            alphas.append(alpha)
        ax.plot(xs, phis, color=BETA_COLORS[beta], linewidth=0.8, zorder=1,
                linestyle="--" if beta == 0.0 else "-")
        for x, phi, alpha in zip(xs, phis, alphas):
            ax.plot([x], [phi], marker=BETA_MARKERS[beta], color=ALPHA_COLORS[alpha],
                    markersize=5, zorder=2, linestyle="none")

    # The alpha = beta = 0 case attaches upstream of the first ring: the
    # honest apex answer, drawn open rather than pretending to a location.
    ax.plot([0.0], [0.0], marker="o", markerfacecolor="white", markeredgecolor=INK,
            markersize=5, zorder=3, linestyle="none")
    ax.annotate("apex\n(attaches upstream)", xy=(0.0, 0.0), xytext=(6, -18),
                textcoords="offset points", fontsize=6.8, color=INK)

    # Positive beta blows toward +y, so its stagnation point sits on the
    # LEFT (negative-phi) flank; the two families mirror exactly.
    ax.annotate(r"$\beta = -10^\circ$", xy=(0.0053, 36), fontsize=7.5, color=BETA_COLORS[-10.0])
    ax.annotate(r"$\beta = 0$", xy=(0.0032, -12), fontsize=7.5, color=BETA_COLORS[0.0])
    ax.annotate(r"$\beta = +10^\circ$", xy=(0.0058, -70), fontsize=7.5, color=BETA_COLORS[10.0])
    for alpha in (0.0, 4.0, 8.0):
        ax.plot([], [], marker="s", color=ALPHA_COLORS[alpha], linestyle="none",
                markersize=4, label=rf"$\alpha = {alpha:.0f}^\circ$")
    ax.legend(loc="upper right", frameon=False, fontsize=7, handletextpad=0.2)

    ax.set_xlabel(r"$x/L$ (from nose)", color=INK)
    ax.set_ylabel(r"meridian angle from windward $\phi$ [deg]", color=INK)
    ax.set_xlim(-0.0004, 0.008)
    ax.set_ylim(-100, 100)
    style_axes(ax)
    fig.subplots_adjust(left=0.16, right=0.97, bottom=0.17, top=0.96)
    return fig


# ------------------------------------------ fig 4: streamlines colored by Cp
def fig_body_streamlines(data):
    meta = data["meta"]
    cond = condition(data, meta["figureAlphaDeg"], meta["figureBetaDeg"])
    surface = cond["surface"]
    stations, sectors = surface["stations"], surface["sectors"]
    pts = np.array(surface["points"]).reshape(stations, sectors, 4)
    length = pts[:, :, 0].max()

    fig, (top, bottom) = plt.subplots(
        2, 1, figsize=(5.0, 5.4), height_ratios=[1.35, 1.0])

    # --- (a) orthographic three-quarter view from the windward-left side ----
    view = np.array([-0.45, -0.72, -0.53])  # camera direction: forward, left, below
    view /= np.linalg.norm(view)
    up_seed = np.array([0.0, 0.0, 1.0])
    right = np.cross(up_seed, view)
    right /= np.linalg.norm(right)
    up = np.cross(view, right)

    def project(p):
        return np.dot(p, right), np.dot(p, up), np.dot(p, view)

    quads, depth, face = [], [], []
    cp_min, cp_max = -0.6, 1.0
    norm = matplotlib.colors.TwoSlopeNorm(vmin=cp_min, vcenter=0.0, vmax=cp_max)
    xyz = pts[:, :, :3]
    for i in range(stations - 1):
        for j in range(sectors):
            j1 = (j + 1) % sectors
            corners = [xyz[i, j], xyz[i + 1, j], xyz[i + 1, j1], xyz[i, j1]]
            projected = [project(c) for c in corners]
            quads.append([(u, v) for u, v, _ in projected])
            depth.append(np.mean([d for _, _, d in projected]))
            cp_mean = np.mean([pts[i, j, 3], pts[i + 1, j, 3],
                               pts[i + 1, j1, 3], pts[i, j1, 3]])
            face.append(CMAP(norm(np.clip(cp_mean, cp_min, cp_max))))

    # The wing lattice, colored by its solved bound circulation -- the same
    # convention the companion SciTech figure uses -- painter-sorted into
    # one scene with the body.
    wing = cond["wingPanels"]
    gamma = np.array([p["gamma"] for p in wing])
    wing_cmap = plt.get_cmap("Greens")
    # Min-to-max stretch: a rectangular wing's loading varies little over
    # mid-span, and anchoring the ramp at zero would paint the whole
    # surface one green. The colorbar states the true range.
    gamma_norm = matplotlib.colors.Normalize(vmin=float(gamma.min()), vmax=float(gamma.max()))
    for panel in wing:
        corners = np.array(panel["corners"])
        projected = [project(c) for c in corners]
        quads.append([(u, v) for u, v, _ in projected])
        depth.append(np.mean([d for _, _, d in projected]))
        face.append(wing_cmap(gamma_norm(panel["gamma"])))

    order = np.argsort(depth)  # painter: far first
    collection = PolyCollection([quads[k] for k in order],
                                facecolors=[face[k] for k in order],
                                edgecolors="none", zorder=1)
    top.add_collection(collection)

    # The wing's attachment line, in the 3-D positions the analysis reports.
    # The camera looks up from the windward side, so the line -- which sits
    # on the lower surface -- faces it.
    line_stations = line_arrays(cond)
    for side in (lambda e: e < 0, lambda e: e > 0):
        pts3 = np.array([s["point"] for s in line_stations if side(s["eta"])])
        if len(pts3) < 2:
            continue
        projected = np.array([project(p) for p in pts3])
        # White casing under the ink line, so it reads on any green.
        top.plot(projected[:, 0], projected[:, 1], color="white", linewidth=1.6, zorder=3)
        top.plot(projected[:, 0], projected[:, 1], color="0.1", linewidth=0.7, zorder=4)
    tip = np.array(project(np.array(line_stations[0]["point"])))
    top.annotate("attachment line", xy=(tip[0], tip[1]), xytext=(-46, -12),
                 textcoords="offset points", fontsize=7, color=INK)

    # A streamline point on a body of revolution faces the camera when its
    # outward radial direction points against the view direction.
    for line in cond["streamlines"]:
        arr = np.array(line["points"])
        radial = arr[:, :3].copy()
        radial[:, 0] = 0.0
        norms = np.linalg.norm(radial, axis=1)
        norms[norms < 1e-12] = 1.0
        radial /= norms[:, None]
        visible = radial @ view < -0.05
        projected = np.array([project(p) for p in arr[:, :3]])
        u = np.where(visible, projected[:, 0], np.nan)
        v = np.where(visible, projected[:, 1], np.nan)
        top.plot(u, v, color="0.1", linewidth=0.45, zorder=2)

    top.set_aspect("equal")
    top.axis("off")

    # --- (b) unwrapped meridian map -----------------------------------------
    def meridian_deg(y, z):
        """Angle from the windward-bottom meridian, positive toward +y."""
        return np.degrees(np.arctan2(y, -z))

    x_over_l = pts[:, 0, 0] / length
    phi_columns = meridian_deg(xyz[0, :, 1], xyz[0, :, 2])
    order = np.argsort(phi_columns)
    phi_sorted = phi_columns[order]
    cp_sorted = pts[:, order, 3]
    # Pad one wrapped column on each side so the +-180 seam has no gap.
    phi_padded = np.concatenate(([phi_sorted[-1] - 360.0], phi_sorted, [phi_sorted[0] + 360.0]))
    cp_padded = np.concatenate((cp_sorted[:, -1:], cp_sorted, cp_sorted[:, :1]), axis=1)
    mesh = bottom.pcolormesh(
        np.repeat(x_over_l[:, None], cp_padded.shape[1], axis=1),
        np.repeat(phi_padded[None, :], stations, axis=0),
        cp_padded, cmap=CMAP, norm=norm, shading="gouraud", zorder=1)

    for line in cond["streamlines"]:
        arr = np.array(line["points"])
        x = arr[:, 0] / length
        phi = meridian_deg(arr[:, 1], arr[:, 2])
        # Break the curve where it wraps across +-180.
        phi[1:][np.abs(np.diff(phi)) > 180.0] = np.nan
        bottom.plot(x, phi, color="0.1", linewidth=0.45, zorder=2)

    nodes = [p for p in cond["body"]["criticalPoints"] if p["type"] == "attachment_node"]
    if nodes:
        x, y, z = nodes[0]["point"]
        bottom.plot([x / length], [np.degrees(np.arctan2(y, -z))], marker="o",
                    markerfacecolor="white", markeredgecolor="0.1", markersize=5, zorder=3)
        bottom.annotate("attachment node", xy=(x / length, np.degrees(np.arctan2(y, -z))),
                        xytext=(10, 10), textcoords="offset points", fontsize=7, color=INK)

    bottom.set_xlabel(r"$x/L$ (from nose)", color=INK)
    bottom.set_ylabel(r"meridian angle from windward $\phi$ [deg]", color=INK)
    bottom.set_xlim(0.0, 1.0)
    bottom.set_ylim(-180, 180)
    bottom.set_yticks([-180, -90, 0, 90, 180])
    style_axes(bottom)

    gamma_bar = fig.colorbar(
        matplotlib.cm.ScalarMappable(norm=gamma_norm, cmap=wing_cmap),
        ax=top, fraction=0.04, pad=0.02)
    gamma_bar.set_label(r"$\Gamma$ [m$^2$/s]", color=INK)
    gamma_bar.ax.tick_params(colors=INK, labelsize=7)
    gamma_bar.outline.set_visible(False)

    bar = fig.colorbar(mesh, ax=bottom, fraction=0.04, pad=0.02)
    bar.set_label(r"$C_p$", color=INK)
    bar.ax.tick_params(colors=INK, labelsize=7)
    bar.outline.set_visible(False)

    fig.subplots_adjust(left=0.11, right=0.86, bottom=0.10, top=0.98, hspace=0.15)
    return fig


def main():
    data = load()
    for name, build in (
        ("attachment-offset", fig_attachment_offset),
        ("sideslip-asymmetry", fig_sideslip),
        ("body-topology", fig_body_topology),
        ("body-streamlines", fig_body_streamlines),
    ):
        fig = build(data)
        fig.savefig(HERE / f"{name}.pdf")
        fig.savefig(HERE / f"{name}.png", dpi=220)
        plt.close(fig)
        print(f"wrote {name}.pdf / .png")


if __name__ == "__main__":
    main()
