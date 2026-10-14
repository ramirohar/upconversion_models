"""Digitize Yu et al. 2014 (sources/yu2014.pdf) Fig. 4, bulk β-NaYF4 (black squares) only.

fig4.png is the figure image extracted from the PDF (page 5). Axes are calibrated from tick pixel positions.
Markers are found as connected components of black pixels after a morphological opening (removes the
connecting lines). Squares partly covered by other markers are located from their visible edge
(marker side ~15 px); squares fully hidden are NaN.
Excitation: 980 nm, 2.9 W/mm² = 290 W/cm², in vacuum.

Usage: pixi run python sources/yu2014_digitized/digitize.py sources/yu2014_digitized sources/yu2014_digitized
"""
import sys, numpy as np, pandas as pd
from PIL import Image
from scipy import ndimage
figs, out = sys.argv[1], sys.argv[2]
a = np.asarray(Image.open(f"{figs}/fig4.png").convert("RGB")).astype(int)
black = (a.max(-1) < 110) & (a.max(-1) - a.min(-1) < 50)

Ts = np.array([10, 30, 50, 100, 150, 200, 250, 300, 350, 400])
half = 7.5  # half marker side [px]

# name: frame rows, (px, value) of top and zero y ticks, (px at 0 K, px at 400 K), partial markers
panels = {
    # 350 K: bottom hidden by the 6 nm triangle -> center from top edge
    "blue": ((71, 442), (94.5, 7e3), 419, (221.5, 622.5), {350: "top"}),
    # 350 K: top hidden by the 25 nm triangle -> center from bottom edge
    "S": ((619, 969), (619, 6e4), 943.5, (221.5, 622.5), {350: "bottom"}),
    # 300 K: bottom hidden; 400 K: top hidden (25 nm triangles)
    "red": ((1179, 1561), (1179.5, 2.5e4), 1525, (223, 639), {300: "top", 400: "bottom"}),
}

df = pd.DataFrame(index=pd.Index(Ts, name="T"))
for name, ((y0, y1), (ptop, vtop), pzero, (x0, x4), partial) in panels.items():
    m = np.zeros_like(black)
    m[y0 + 3 : y1 - 3, 176:668] = black[y0 + 3 : y1 - 3, 176:668]
    m[y0 : y0 + 120, 470:668] = False  # legend
    m = ndimage.binary_opening(m, np.ones((4, 4)))
    lab, _ = ndimage.label(m)
    col = {}
    for sl in ndimage.find_objects(lab):
        w, h = sl[1].stop - sl[1].start, sl[0].stop - sl[0].start
        T = Ts[np.argmin(abs(Ts - ((sl[1].start + sl[1].stop - 1) / 2 - x0) * 400 / (x4 - x0)))]
        if T in partial and w >= 12:
            yc = sl[0].start + half if partial[T] == "top" else sl[0].stop - 1 - half
        elif w >= 13 and h >= 13:
            yc = (sl[0].start + sl[0].stop - 1) / 2
        else:
            continue
        col[T] = (pzero - yc) / (pzero - ptop) * vtop
    df[name] = pd.Series(col)

df.to_csv(f"{out}/fig4_bulk.csv", float_format="%.4g")
print(df)
