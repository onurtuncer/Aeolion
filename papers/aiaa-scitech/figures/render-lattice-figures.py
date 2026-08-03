# render-lattice-figures.py -- publication figure for the AIAA SciTech
# paper: the reference ducted-fan configuration (two-blade rotor, snug
# shroud, four-vane cruciform) drawn from the solver's own exported
# lattice and solved circulation, as a vector PDF.
#
# Input:  lattice-solution.json, written by `aeolion_lattice_export`
#         (app/LatticeFigureExport.cpp) -- rebuild and rerun that tool
#         whenever the reference case changes, then rerun this script.
# Output: lattice-views.pdf (for paper.tex) + lattice-views.png (preview).
#
# Panel (a): aft-quarter orthographic view, hand-projected and
# painter-sorted (global per-quad depth sort; the geometry has no
# interpenetrating panels, so the painter's algorithm suffices). A true
# aft view is useless here -- the undeflected vanes are edge-on to it --
# so the camera sits a quarter off the axis, downstream.
#
# Panel (b): meridional section in the x-z plane, the classic
# presentation for an axisymmetric machine: hatched duct wall cuts, the
# z-axis vane pair at true profile (those plates LIE in the section
# plane), the rotor disk indicated edge-on, dash-dot centerline.
#
# Blades and vanes are colored by solved bound circulation on one
# diverging scale centered at zero -- the two populations carry opposite
# signs, the vanes recovering the rotor's swirl -- and the shroud is a
# neutral Lambert-shaded gray, carrying no field.

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.patches import FancyArrowPatch, Rectangle

HERE = Path(__file__).parent
DATA = HERE / "lattice-solution.json"

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman", "STIXGeneral", "DejaVu Serif"],
    "mathtext.fontset": "stix",
    "font.size": 8.0,
    "pdf.fonttype": 42,
    "hatch.linewidth": 0.4,
})

INK = "0.25"          # annotation ink: recessive gray, never pure black
DUCT_BASE = np.array([0.82, 0.82, 0.835])
DUCT_EDGE = "0.55"
LIGHT = np.array([0.35, -0.45, 0.82])
LIGHT = LIGHT / np.linalg.norm(LIGHT)
CMAP = plt.get_cmap("RdBu_r")


def load():
    d = json.loads(DATA.read_text())
    lifting = d["bladePanels"] + d["vanePanels"]
    gmax = max(abs(p["gamma"]) for p in lifting)
    return d["meta"], lifting, d["ductPanels"], gmax


def gamma_color(gamma, glim):
    return CMAP(0.5 * (1.0 + gamma / glim))


def duct_facecolor(normal):
    shade = 0.62 + 0.38 * abs(float(np.dot(normal, LIGHT)))
    return np.clip(DUCT_BASE * shade, 0.0, 1.0)


def view_basis(toward_camera):
    """Orthographic screen basis for a camera looking along -toward_camera,
    with world +z as the up reference."""
    toward = np.asarray(toward_camera, float)
    toward = toward / np.linalg.norm(toward)
    right = np.cross([0.0, 0.0, 1.0], toward)
    right = right / np.linalg.norm(right)
    up = np.cross(toward, right)
    return right, up, toward


def project(points, right, up):
    pts = np.asarray(points, float)
    return np.column_stack((pts @ right, pts @ up))


# --- panel (a): aft-quarter view -------------------------------------------

def draw_quarter_view(ax, lifting, duct, glim, meta):
    right, up, toward = view_basis((0.82, 0.47, 0.33))
    quads = []  # (depth, verts2d, facecolor, edgecolor, linewidth)

    for p in duct:
        c = np.asarray(p["corners"], float)
        quads.append((float(c.mean(axis=0) @ toward), project(c, right, up),
                      duct_facecolor(np.asarray(p["normal"], float)), DUCT_EDGE, 0.25))
    for p in lifting:
        c = np.asarray(p["corners"], float)
        quads.append((float(c.mean(axis=0) @ toward), project(c, right, up),
                      gamma_color(p["gamma"], glim), "0.78", 0.3))

    quads.sort(key=lambda q: q[0])
    ax.add_collection(PolyCollection(
        [q[1] for q in quads], facecolors=[q[2] for q in quads],
        edgecolors=[q[3] for q in quads], linewidths=[q[4] for q in quads],
        joinstyle="round"))

    # rotation-direction arc, drawn in the rotor plane and projected like
    # the geometry (+Omega about +x carries y into z)
    r_arc = 1.13 * meta["shroudOuter"]
    theta = np.linspace(np.deg2rad(96), np.deg2rad(148), 32)
    arc3d = np.column_stack((np.zeros_like(theta),
                             r_arc * np.cos(theta), r_arc * np.sin(theta)))
    arc = project(arc3d, right, up)
    ax.plot(arc[:, 0], arc[:, 1], color=INK, lw=0.7, solid_capstyle="round")
    ax.add_patch(FancyArrowPatch(arc[-2], arc[-1], arrowstyle="-|>",
                                 mutation_scale=6, color=INK, lw=0.7))
    label = project([[0.0, 1.16 * r_arc * np.cos(np.deg2rad(122)),
                      1.16 * r_arc * np.sin(np.deg2rad(122))]], right, up)[0]
    ax.text(label[0], label[1], r"$\Omega$", ha="center", va="center",
            fontsize=9, color=INK)

    ax.set_aspect("equal")
    ax.autoscale_view()
    ax.set_axis_off()


# --- panel (b): meridional section -----------------------------------------

def draw_section(ax, lifting, meta, glim):
    le, te = -0.5 * meta["shroudChord"], 0.5 * meta["shroudChord"]
    bore, outer = meta["shroudInner"], meta["shroudOuter"]
    tip = meta["radius"]

    # dash-dot centerline, the axisymmetry convention
    ax.plot([le - 0.030, 0.185], [0.0, 0.0], color="0.55", lw=0.6,
            dashes=(6, 1.5, 1, 1.5), zorder=0)

    # duct wall cross-sections, hatched as cut material
    for sign in (+1, -1):
        ax.add_patch(Rectangle((le, sign * bore), te - le, sign * (outer - bore),
                               facecolor="0.90", edgecolor="0.45", lw=0.6,
                               hatch="/////"))

    # the rotor disk, edge-on: indicated, not resolved (the blades are
    # perpendicular to the section plane; panel (a) carries their field)
    ax.add_patch(Rectangle((-0.012, -tip), 0.024, 2 * tip,
                           facecolor="0.80", edgecolor="0.45", lw=0.6))

    # vane-hinge / duct-exit plane
    ax.plot([te, te], [-bore, bore], color="0.55", lw=0.5, ls=(0, (1.5, 2)),
            zorder=0)

    # the z-axis vane pair lies in the section plane: true-size profiles
    quads, colors = [], []
    for p in lifting:
        c = np.asarray(p["corners"], float)
        if not p["surface"].startswith("vane") or abs(c.mean(axis=0)[1]) > 0.01:
            continue
        quads.append(c[:, [0, 2]])
        colors.append(gamma_color(p["gamma"], glim))
    ax.add_collection(PolyCollection(quads, facecolors=colors,
                                     edgecolors="0.78", linewidths=0.3))

    # annotations
    ann = dict(fontsize=7.5, color=INK)
    ax.annotate("duct", xy=(0.2 * te, bore + 0.55 * (outer - bore)),
                xytext=(-0.055, outer + 0.028), ha="center", **ann,
                arrowprops=dict(arrowstyle="-", lw=0.5, color=INK,
                                shrinkA=2, shrinkB=2))
    ax.annotate("rotor", xy=(-0.006, -0.62 * tip), xytext=(-0.062, -0.66 * tip),
                ha="right", va="center", **ann,
                arrowprops=dict(arrowstyle="-", lw=0.5, color=INK,
                                shrinkA=1, shrinkB=2))
    ax.annotate("vane", xy=(0.72 * te + 0.03, 0.80 * tip), ha="center",
                xytext=(0.115, outer + 0.028), **ann,
                arrowprops=dict(arrowstyle="-", lw=0.5, color=INK,
                                shrinkA=2, shrinkB=2))
    ax.text(te - 0.005, -0.042, "exit plane", ha="right", va="top",
            rotation=90, fontsize=6.5, color="0.45")

    # jet-direction arrow on the centerline
    ax.add_patch(FancyArrowPatch((0.135, 0.0), (0.175, 0.0), arrowstyle="-|>",
                                 mutation_scale=7, color=INK, lw=0.8))
    ax.text(0.155, 0.010, "jet", ha="center", va="bottom", fontsize=8,
            color=INK)

    # scale bar
    sx, sz = -0.0375, -outer - 0.052
    ax.plot([sx, sx + 0.05], [sz, sz], color=INK, lw=0.9, solid_capstyle="butt")
    for xt in (sx, sx + 0.05):
        ax.plot([xt, xt], [sz - 0.0045, sz + 0.0045], color=INK, lw=0.9)
    ax.text(sx + 0.025, sz - 0.011, "0.05 m", ha="center", va="top",
            fontsize=7.5, color=INK)

    ax.set_aspect("equal")
    ax.set_xlim(-0.115, 0.195)
    ax.set_ylim(-0.275, 0.235)
    ax.set_axis_off()


def main():
    meta, lifting, duct, gmax = load()
    glim = float(np.ceil(gmax * 20.0) / 20.0)  # symmetric limit, rounded up to 0.05

    fig = plt.figure(figsize=(6.5, 3.45))
    ax_a = fig.add_axes([0.005, 0.20, 0.44, 0.79])
    ax_b = fig.add_axes([0.46, 0.20, 0.535, 0.79])
    ax_cbar = fig.add_axes([0.355, 0.125, 0.29, 0.025])

    draw_quarter_view(ax_a, lifting, duct, glim, meta)
    draw_section(ax_b, lifting, meta, glim)

    sm = plt.cm.ScalarMappable(cmap=CMAP, norm=plt.Normalize(-glim, glim))
    cbar = fig.colorbar(sm, cax=ax_cbar, orientation="horizontal")
    cbar.set_ticks([-0.4, -0.2, 0.0, 0.2, 0.4])
    cbar.ax.tick_params(labelsize=7, length=2, width=0.5, color=INK,
                        labelcolor=INK)
    cbar.outline.set_linewidth(0.5)
    cbar.outline.set_edgecolor(INK)
    cbar.set_label(r"bound circulation $\Gamma$ (m$^2$/s)", fontsize=8,
                   color=INK, labelpad=2)

    fig.text(0.225, 0.04, "(a)", ha="center", fontsize=8.5, color=INK)
    fig.text(0.73, 0.04, "(b)", ha="center", fontsize=8.5, color=INK)

    fig.savefig(HERE / "lattice-views.pdf")
    fig.savefig(HERE / "lattice-views.png", dpi=300)
    print("wrote", HERE / "lattice-views.pdf")


if __name__ == "__main__":
    main()
