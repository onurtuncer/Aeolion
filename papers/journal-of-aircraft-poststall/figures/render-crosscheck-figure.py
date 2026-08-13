"""The 3-D cross-check against the quasi-steady map: wing-only normal force at
the declared attitudes, map cycle means vs particle-wake mean +- RMS.

Inputs (this directory):
  particle-crosscheck.json        -- aeolion_particle_crosscheck
  poststall-map-anchored.json  -- the map (wing-only forces per condition)
"""
import json
import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
BLUE = "#2a78d6"    # quasi-steady map (cycle means)
AQUA = "#1baf7a"    # 3-D particle cross-check (mean +- RMS)

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "DejaVu Sans"],
    "text.color": INK, "axes.edgecolor": BASELINE, "axes.labelcolor": SECONDARY,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "font.size": 9,
})

ref = json.load(open("particle-crosscheck.json"))
rows = np.array(ref["crosscheck"])  # a, b, CL, CD, CY, CN, rmsCL, rmsCN, circCL, particles

mapd = json.load(open("poststall-map-anchored.json"))
qS = 0.5 * mapd["meta"]["rho"] * mapd["meta"]["Vinf"] ** 2 * mapd["meta"]["area"]
map_cn = {}
for c in mapd["conditions"]:
    map_cn[(c["alphaDeg"], c["betaDeg"])] = c["forceWing"][2] / qS

fig, ax = plt.subplots(figsize=(4.6, 3.0), dpi=300)
markers = {0.0: "o", 15.0: "^"}
off = {0.0: -1.0, 15.0: 1.0}
for b in (0.0, 15.0):
    sel = rows[rows[:, 1] == b]
    ax.errorbar(sel[:, 0] + off[b], sel[:, 5], yerr=sel[:, 7], fmt=markers[b], ms=4.0,
                color=AQUA, elinewidth=1.0, capsize=2.0, mfc=AQUA if b == 0.0 else "none")
    a_map = [a for (a, bb) in map_cn if bb == b and a in sel[:, 0]]
    a_map.sort()
    ax.plot([a + off[b] for a in a_map], [map_cn[(a, b)] for a in a_map], markers[b],
            ms=4.5, color=BLUE, mfc=BLUE if b == 0.0 else "none")

ax.set_xlabel(r"$\alpha$ [deg]")
ax.set_ylabel(r"wing-only $C_N$")
ax.set_xticks([30, 60, 90])
ax.text(31, 0.6, "map (cycle mean)", color=BLUE, fontsize=8)
ax.text(31, 0.4, "3-D cross-check (mean $\\pm$ RMS)", color=AQUA, fontsize=8)
ax.text(31, 0.2, "filled $\\beta=0$, open $\\beta=15^\\circ$", color=SECONDARY, fontsize=8)

fig.tight_layout()
fig.savefig("particle-crosscheck.pdf", bbox_inches="tight", facecolor=SURFACE)
fig.savefig("particle-crosscheck.png", bbox_inches="tight", facecolor=SURFACE)

print("attitude  beta   map CN    cross-check CN (rms)")
for a, b, cl, cd, cy, cn, rcl, rcn, ccl, npart in rows:
    m = map_cn.get((a, b), float("nan"))
    print(f"  {a:4.0f}  {b:5.0f}  {m:7.3f}   {cn:7.3f} ({rcn:.3f})   [{int(npart)} particles]")
print("wrote particle-crosscheck.pdf / .png")
