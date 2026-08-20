"""Section-level validation figure: the anchored post-stall model against
the Sheldahl-Klimas NACA 0015 table at Re = 3.6e5.

Inputs (this directory):
  section-validation.json        -- aeolion_section_validation's model polar
  NACA_0015_SheldahlKlimas.dat   -- Sandia's distributed 360-deg table
                                    (github.com/sandialabs/CACTUS,
                                    test/Airfoil_Section_Data/NACA_0015.dat;
                                    the data of SAND80-2114)

Output: section-validation.pdf/.png and the headline numbers for the prose
(printed to stdout).
"""
import json
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RE_BLOCK = "3.6e5"

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASELINE = "#c3c2b7"
BLUE = "#2a78d6"    # anchored model
ORANGE = "#eb6834"  # Sheldahl-Klimas table

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Segoe UI", "DejaVu Sans"],
    "text.color": INK, "axes.edgecolor": BASELINE, "axes.labelcolor": SECONDARY,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True,
    "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE,
    "axes.spines.top": False, "axes.spines.right": False,
    "legend.frameon": False, "font.size": 9,
})

model = json.load(open("section-validation.json"))
polar = np.array(model["polar"])  # alpha, f, cl, cd, cm
meta = model["meta"]

txt = open("NACA_0015_SheldahlKlimas.dat").read()
block = re.split(r"Reynolds Number: ", txt)
data = None
for b in block[1:]:
    if b.split()[0] == RE_BLOCK:
        rows = [l.split() for l in b.splitlines() if l and l[0] in "-0123456789"]
        rows = [r for r in rows if len(r) == 4]
        data = np.array([[float(v) for v in r] for r in rows])  # aoa cl cd cm25
        break
assert data is not None, f"Re block {RE_BLOCK} not found"
mask = (data[:, 0] >= 0.0) & (data[:, 0] <= 90.0)
sk = data[mask]

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(6.6, 2.9), dpi=300)

ax1.plot(polar[:, 0], polar[:, 2], color=BLUE, lw=1.6)
ax1.plot(sk[:, 0], sk[:, 1], "o", ms=3.2, mfc="none", mec=ORANGE, mew=1.0)
ax1.set_xlabel(r"$\alpha$ [deg]")
ax1.set_ylabel(r"$c_l$")
ax1.set_xlim(0, 90)
ax1.text(2, 1.42, "anchored model", color=BLUE, fontsize=8)
ax1.text(38, 0.55, "Sheldahl--Klimas", color=ORANGE, fontsize=8)

ax2.plot(polar[:, 0], polar[:, 3], color=BLUE, lw=1.6)
ax2.plot(sk[:, 0], sk[:, 2], "o", ms=3.2, mfc="none", mec=ORANGE, mew=1.0)
ax2.set_xlabel(r"$\alpha$ [deg]")
ax2.set_ylabel(r"$c_d$")
ax2.set_xlim(0, 90)

fig.tight_layout()
fig.savefig("section-validation.pdf", bbox_inches="tight", facecolor=SURFACE)
fig.savefig("section-validation.png", bbox_inches="tight", facecolor=SURFACE)

# --- headline numbers for the prose ---------------------------------------
mi = np.argmax(polar[:, 2])
di = np.argmax(sk[:, 1])
print(f"model:  clmax {polar[mi, 2]:.3f} at {polar[mi, 0]:.0f} deg; stall (meta) {meta['stallDeg']:.1f}")
print(f"data:   clmax {sk[di, 1]:.3f} at {sk[di, 0]:.0f} deg")
cd90_m = polar[polar[:, 0] == 90.0, 3][0]
cd90_d = sk[sk[:, 0] == 90.0, 2][0]
print(f"cd(90): model {cd90_m:.3f}  data {cd90_d:.3f}  ({100 * (cd90_m / cd90_d - 1):+.1f}%)")

# CN comparison beyond 40 deg (the declared acceptance metric)
a = polar[:, 0]
cn_m = polar[:, 2] * np.cos(np.radians(a)) + polar[:, 3] * np.sin(np.radians(a))
cn_d = sk[:, 1] * np.cos(np.radians(sk[:, 0])) + sk[:, 2] * np.sin(np.radians(sk[:, 0]))
cn_mi = np.interp(sk[:, 0], a, cn_m)
sel = sk[:, 0] >= 40.0
err = (cn_mi[sel] - cn_d[sel]) / cn_d[sel]
print(f"CN, 40..90 deg: mean {100 * np.mean(err):+.1f}%  max |{100 * np.max(np.abs(err)):.1f}|%")
selmid = (sk[:, 0] >= 25.0) & (sk[:, 0] < 40.0)
errmid = (cn_mi[selmid] - cn_d[selmid]) / cn_d[selmid]
print(f"CN, 25..40 deg: mean {100 * np.mean(errmid):+.1f}%  max |{100 * np.max(np.abs(errmid)):.1f}|%")
print("wrote section-validation.pdf / .png")
