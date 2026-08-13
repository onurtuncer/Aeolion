"""Method sketches for the unsteady cross-check (Part II).

Left: the 2-D LESP-modulated discrete-vortex section -- Glauert bound
sheet on the chord, trailing-edge vortices every step, leading-edge
vortices while |A0| exceeds the critical LESP, loads from the impulse of
the tracked vorticity.

Right: the 3-D particle-wake solve -- one ring per strip, the one-step
buffer ring whose near segment cancels the ring's closing segment, the
merged conversion to particles, and the counter-rotating leading-edge
flux pair on separated strips.

Hand-drawn with matplotlib primitives; palette matches the other figures.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, Polygon

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
SECONDARY = "#52514e"
MUTED = "#898781"
BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "DejaVu Sans"],
    "text.color": INK, "figure.facecolor": SURFACE, "font.size": 8.5,
})

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.6, 2.7), dpi=300)
for ax in (ax1, ax2):
    ax.set_facecolor(SURFACE)
    ax.set_aspect("equal")
    ax.axis("off")

# ------------------------------------------------------------- 2-D panel
alpha = np.radians(25.0)
chord = np.array([np.cos(alpha), -np.sin(alpha)])
le = np.array([0.0, 0.0])
te = le + chord

ax1.plot([le[0], te[0]], [le[1], te[1]], color=INK, lw=2.2, solid_capstyle="round")

# Freestream arrow.
ax1.add_patch(FancyArrowPatch((-0.55, 0.12), (-0.15, 0.12), arrowstyle="-|>",
                              mutation_scale=9, color=SECONDARY, lw=1.2))
ax1.text(-0.55, 0.17, r"$U_\infty$", color=SECONDARY)

# Bound sheet: Glauert-weighted circulation ticks along the chord.
for s in np.linspace(0.08, 0.92, 9):
    p = le + chord * s
    strength = 0.09 * np.sqrt((1 - s) / max(s, 0.02))
    ax1.add_patch(plt.Circle(p, min(strength, 0.075), fill=False, color=BLUE, lw=1.0))
ax1.text(0.16, 0.14, r"bound sheet $\gamma(\theta,t)$", color=BLUE)

# TEV chain: convecting downstream of the TE.
rng = np.random.default_rng(4)
for k in range(8):
    p = te + np.array([0.10 + 0.13 * k, 0.015 * k + 0.02 * rng.standard_normal()])
    ax1.add_patch(plt.Circle(p, 0.028, color=ORANGE, alpha=0.85))
ax1.text(te[0] + 0.28, te[1] - 0.20, "TEV / step\n(Kelvin closes it)", color=ORANGE,
         ha="center")

# LEV chain: rolling up above the plate.
for k in range(6):
    p = le + np.array([0.12 + 0.11 * k, 0.16 + 0.05 * k + 0.02 * rng.standard_normal()])
    ax1.add_patch(plt.Circle(p, 0.028, fill=False, color=ORANGE, lw=1.3))
ax1.text(0.42, 0.58, r"LEV while $|A_0| > \mathrm{LESP}_{crit}$", color=ORANGE)

ax1.set_xlim(-0.65, 1.75)
ax1.set_ylim(-0.55, 0.75)
ax1.set_title("2-D: LESP-modulated discrete vortices", fontsize=9, color=INK, loc="left")

# ------------------------------------------------------------- 3-D panel
# Planform in shear projection: span to the right, chord down-and-right.
def proj(x, y):
    return np.array([0.60 * y + 0.34 * x, -0.62 * x])

nstrips = 5
for i in range(nstrips):
    quad = [proj(0.25, i), proj(0.25, i + 1), proj(1.0, i + 1), proj(1.0, i)]
    ax2.add_patch(Polygon(quad, closed=True, fill=False, color=BLUE, lw=1.1))
b0, b1 = proj(0.25, 0), proj(0.25, nstrips)
ax2.plot([b0[0], b1[0]], [b0[1], b1[1]], color=BLUE, lw=2.2)
lbl = proj(0.55, -0.15)
ax2.text(lbl[0] - 1.15, lbl[1], "rings $\\Gamma_i(t)$", color=BLUE)

# Buffer row, dashed, one convection step deep.
for i in range(nstrips):
    quad = [proj(1.0, i), proj(1.0, i + 1), proj(1.3, i + 1), proj(1.3, i)]
    ax2.add_patch(Polygon(quad, closed=True, fill=False, color=AQUA, lw=1.0, ls=(0, (3, 2))))
blbl = proj(1.15, nstrips + 0.15)
ax2.text(blbl[0] + 0.05, blbl[1] + 0.05, "buffer ring $\\Gamma_i(t{-}1)$\n(cancels the TE closer)",
         color=AQUA, fontsize=7.5, va="top")

# Trailing-edge particle stream: fans down-right from the buffer's far edge,
# aging particles fading and coarsening.
chord_dir = proj(1.0, 0.0) - proj(0.0, 0.0)
chord_dir = chord_dir / np.linalg.norm(chord_dir)
for k in range(34):
    age = k // 8
    y = rng.uniform(0.1, nstrips - 0.1)
    start = proj(1.3, y)
    p = start + chord_dir * (0.14 + 0.34 * age + 0.10 * (k % 8) / 8.0)
    p = p + (0.03 + 0.03 * age) * rng.standard_normal(2)
    ax2.add_patch(plt.Circle(p, 0.024 + 0.006 * age, color=ORANGE, alpha=0.85 - 0.16 * age))
# Leading-edge stream on the separated strips: lifts over the top, then bends
# downstream -- the counter-rotating layer (open circles).
for k in range(16):
    age = k // 5
    y = rng.uniform(1.2, nstrips - 0.5)
    start = proj(0.0, y)
    p = start + np.array([0.16, 0.24]) * (0.4 + 0.5 * age) + chord_dir * (0.30 * age)
    p = p + (0.03 + 0.02 * age) * rng.standard_normal(2)
    ax2.add_patch(plt.Circle(p, 0.026 + 0.005 * age, fill=False, color=ORANGE,
                             lw=1.2, alpha=0.9 - 0.15 * age))
le_lbl = proj(0.0, 2.6)
ax2.text(le_lbl[0], le_lbl[1] + 0.55, "LE flux pair, separated strips",
         color=ORANGE, fontsize=7.5)
te_lbl = proj(2.15, 1.1)
ax2.text(te_lbl[0], te_lbl[1] - 0.12, "particles: merged conversion,\nstretching, $N^2$ induction",
         color=ORANGE, fontsize=7.5, va="top")

# Freestream, aligned with the projected chordwise direction.
f0 = proj(-0.85, -0.55)
ax2.add_patch(FancyArrowPatch(f0, f0 + chord_dir * 0.42, arrowstyle="-|>", mutation_scale=9,
                              color=SECONDARY, lw=1.2))
ax2.text(f0[0] - 0.30, f0[1] + 0.02, r"$V_\infty$", color=SECONDARY)

ax2.set_xlim(-1.5, 4.6)
ax2.set_ylim(-2.15, 0.95)
ax2.set_title("3-D: ring lattice + particle wake", fontsize=9, color=INK, loc="left")

fig.tight_layout()
fig.savefig("crosscheck-sketch.pdf", bbox_inches="tight", facecolor=SURFACE)
fig.savefig("crosscheck-sketch.png", bbox_inches="tight", facecolor=SURFACE)
print("wrote crosscheck-sketch.pdf / .png")
