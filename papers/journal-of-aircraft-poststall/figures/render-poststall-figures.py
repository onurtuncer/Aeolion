# Part II figure renderer: the post-stall map figures from the anchored
# solve (poststall-map-anchored.json) with the Phase-0 baseline
# (poststall-map-baseline.json) alongside where the comparison IS the
# result. Regenerate both JSONs with app/PostStallSweepExport.cpp:
#
#   aeolion_poststall_sweep <handoff> figures/poststall-map-anchored.json \
#       25 0.05 0 1000 all anchored
#   aeolion_poststall_sweep <handoff> figures/poststall-map-baseline.json \
#       25 0.05 0 1000 all analytic
#
# Rerun after any change to ViscousCoupling.h, PostStallSection.h, or the
# separation march.
import json
import math
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

INK = "#0b0b0b"
SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
BLUE = "#2a78d6"    # anchored (Part II) solve
ORANGE = "#eb6834"  # Phase-0 baseline / references
RAMP = ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]  # beta = 0/10/20/30
DRAWN_BETAS = [0.0, 10.0, 20.0, 30.0]

plt.rcParams.update({
    "font.family": "serif", "font.size": 9,
    "text.color": INK, "axes.edgecolor": BASELINE, "axes.labelcolor": INK,
    "xtick.color": SECONDARY, "ytick.color": SECONDARY, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.5, "axes.axisbelow": True,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False,
})


def load(path):
    with open(path) as f:
        data = json.load(f)
    by_beta = defaultdict(list)
    for c in data["conditions"]:
        by_beta[c["betaDeg"]].append(c)
    for b in by_beta:
        by_beta[b].sort(key=lambda c: c["alphaDeg"])
    return data["meta"], by_beta


meta, anchored = load("poststall-map-anchored.json")
_, baseline = load("poststall-map-baseline.json")

# --- figure 1: the sigma collapse and the suction loss ----------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.8, 2.9), dpi=300)
for i, b in enumerate(DRAWN_BETAS):
    sig, val = [], []
    for c in anchored[b]:
        s = math.sin(math.radians(c["sigmaDeg"]))
        if c["alphaDeg"] >= 6.0 and abs(s) > 1e-6 and c["alphaDeg"] <= 75.0:
            sig.append(c["sigmaDeg"])
            val.append(c["CN"] / s)
    ax1.plot(sig, val, color=RAMP[i], lw=1.4, solid_capstyle="round",
             label=rf"$\beta={b:.0f}^\circ$")
ax1.axhline(meta["viternaCDmax"], color=MUTED, lw=0.9, ls=(0, (4, 3)))
ax1.text(3, meta["viternaCDmax"] + 0.1, rf"$C_{{d,\max}}(A\!R{{=}}6)={meta['viternaCDmax']:.2f}$",
         color=SECONDARY, fontsize=7.5)
ax1.set_xlabel(r"total inclination $\sigma$ [deg]")
ax1.set_ylabel(r"$C_N/\sin\sigma$")
ax1.legend(fontsize=7, loc="upper right")

for i, b in enumerate(DRAWN_BETAS):
    a = [c["alphaDeg"] for c in anchored[b] if c["alphaDeg"] >= 10 and c["attachedAngleDeg"] > 1e-6]
    s = [c["forceAngleDeg"] / c["attachedAngleDeg"] for c in anchored[b]
         if c["alphaDeg"] >= 10 and c["attachedAngleDeg"] > 1e-6]
    ax2.plot(a, s, color=RAMP[i], lw=1.4, solid_capstyle="round")
ax2.axhline(0.10, color=MUTED, lw=0.9, ls=(0, (4, 3)))
ax2.text(12, 0.12, "10% suction left", color=SECONDARY, fontsize=7.5)
ax2.set_xlabel(r"angle of attack $\alpha$ [deg]")
ax2.set_ylabel("leading-edge suction fraction")
ax2.set_ylim(0, 1.0)
fig.tight_layout()
fig.savefig("poststall-collapse.pdf", bbox_inches="tight")

# --- figure 2: beta=0 polar + xcp walk vs baseline and anchors --------------
fig2, (bx1, bx2) = plt.subplots(1, 2, figsize=(6.8, 2.9), dpi=300)
a_b = [c["alphaDeg"] for c in baseline[0.0]]
a_a = [c["alphaDeg"] for c in anchored[0.0]]

bx1.plot(a_b, [c["CN"] for c in baseline[0.0]], color=ORANGE, lw=1.2, ls=(0, (4, 3)),
         label="hand-tuned blend (baseline)")
bx1.plot(a_a, [c["CN"] for c in anchored[0.0]], color=BLUE, lw=1.4, label="anchored model")
bx1.axhline(meta["viternaCDmax"], color=MUTED, lw=0.9, ls=(0, (1, 2)))
bx1.text(1, meta["viternaCDmax"] - 0.16, rf"Viterna $C_{{d,\max}}={meta['viternaCDmax']:.2f}$",
         color=SECONDARY, fontsize=7.5)
bx1.set_xlabel(r"angle of attack $\alpha$ [deg]")
bx1.set_ylabel(r"$C_N$")
bx1.legend(fontsize=7, loc="lower right")

chord, le = meta["rootChord"], meta["wingLEx"]
bx2.plot(a_b, [(c["xcp"] - le) / chord if c["xcpValid"] else float("nan") for c in baseline[0.0]],
         color=ORANGE, lw=1.2, ls=(0, (4, 3)))
bx2.plot(a_a, [(c["xcp"] - le) / chord if c["xcpValid"] else float("nan") for c in anchored[0.0]],
         color=BLUE, lw=1.4)
alphas = [c["alphaDeg"] for c in anchored[0.0] if c["alphaDeg"] >= 0]
bx2.plot(alphas, [0.5 - 0.75 * math.cos(math.radians(a)) / (4 + math.pi * math.sin(math.radians(a)))
                  for a in alphas], color=MUTED, lw=0.9, ls=(0, (4, 3)))
bx2.text(35, 0.47, "Rayleigh plate", color=SECONDARY, fontsize=7.5)
bx2.set_xlabel(r"angle of attack $\alpha$ [deg]")
bx2.set_ylabel(r"$x_{cp}/c$")
bx2.set_xlim(0, 92)
bx2.set_ylim(0.2, 0.56)
fig2.tight_layout()
fig2.savefig("poststall-beta0.pdf", bbox_inches="tight")
print("wrote poststall-collapse.pdf, poststall-beta0.pdf")
