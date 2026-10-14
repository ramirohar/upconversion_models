"""Digitize Suta et al. 2025 (sources/suta2025.pdf) Fig. 1(c): LIR R21 = I(2H11/2) / I(4S3/2) vs 10^3/T.

fig1.png is the full Fig. 1 image extracted from the PDF (page 4); panel (c) is cropped from it.
Axes are calibrated from tick pixel positions (log y: 120 px per decade). Each point is located from its error
bar: the midpoint between the two caps (bars are symmetric in log), or, where a cap is lost under the red
fit line or a neighbour, from the dot itself. Below 10^3/T ~ 2.4 (T > 420 K) the bars overlap and are skipped.
Sample: microcrystalline β-NaYF4:2% Er, 18% Yb; 980 nm CW, ~5.5 mW/cm²; bands 510-533 / 533-555 nm.
Suta's fit (eqn 1 with k_nr(0) offset): ΔE21 = 652 ± 8 cm-1, C = 1.84 ± 0.07, k_nr(0) = 2.29 ± 0.54 µs-1.

Usage: pixi run python sources/suta2025_digitized/digitize.py sources/suta2025_digitized sources/suta2025_digitized
"""
import sys, numpy as np, pandas as pd
from PIL import Image
from scipy import ndimage
figs, out = sys.argv[1], sys.argv[2]
a = np.asarray(Image.open(f"{figs}/fig1.png").convert("RGB").crop((90, 740, 780, 1320))).astype(int)
black = (a.max(-1) < 110) & (a.max(-1) - a.min(-1) < 40)
m = np.zeros_like(black)
m[36:541, 27:664] = black[36:541, 27:664]
m[215:541, 75:135] = False  # rotated "optimum working range" label

fx = lambda p: 2 + (np.asarray(p) - 66) / 57.2  # 10^3/T: ticks 2..12 at 66..638 px
fy = lambda p: 10 ** (-1 - (np.asarray(p) - 228) / 120)  # 1e-1 at 228 px, 1e-3 at 468 px

points = {}

# 1) caps: horizontal segments 11-20 px wide, paired top/bottom in the same column
lab, _ = ndimage.label(ndimage.binary_opening(m, np.ones((1, 11))))
caps = sorted(
    ((s[1].start + s[1].stop - 1) / 2, (s[0].start + s[0].stop - 1) / 2)
    for s in ndimage.find_objects(lab)
    if 11 <= s[1].stop - s[1].start <= 20
)
used = set()
for i, (x, _) in enumerate(caps):
    ys = sorted((caps[j][1], j) for j in range(len(caps)) if abs(caps[j][0] - x) <= 3 and j not in used)
    for (y1, j1), (y2, j2) in zip(ys, ys[1:]):
        if 6 < y2 - y1 < 45 and j1 not in used and j2 not in used:
            xc = (caps[j1][0] + caps[j2][0]) / 2
            points[round(float(fx(xc)), 2)] = (fy((y1 + y2) / 2), fy(y2), fy(y1), "caps")
            used |= {j1, j2}

# 2) dots: thick part of a vertical bar, for bars whose caps were not paired
vert = ndimage.binary_opening(m, np.ones((12, 1)))
cols = np.where(vert.any(0))[0]
for g in np.split(cols, np.where(np.diff(cols) > 2)[0] + 1):
    xc = g.mean()
    if fx(xc) < 2.4 or any(abs(k - fx(xc)) < 0.1 for k in points):
        continue
    ys = np.where(vert[:, g].any(1))[0]
    top, bot = ys[0], ys[-1]
    w = m[top : bot + 1, int(xc) - 5 : int(xc) + 6].sum(1)
    core = np.where(w >= 5)[0]
    core = core[(core > 3) & (core < len(w) - 4)]
    if len(core) >= 3:
        points[round(float(fx(xc)), 2)] = (fy(top + core.mean()), np.nan, np.nan, "dot")

df = pd.DataFrame(points, index=["LIR", "LIR_lo", "LIR_hi", "method"]).T.sort_index()
df.index.name = "inv_T_1e3"
df.insert(0, "T", (1e3 / df.index.to_numpy()).round(1))
df.to_csv(f"{out}/fig1c_LIR.csv", float_format="%.4g")
print(df)
