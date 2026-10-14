"""Digitize Anderson et al. 2013 (sources/anderson2013.pdf) Figs. 2-4 from the embedded raster images.

figN.png are the figure images extracted from the PDF. Axes are calibrated from axis/gridline/tick pixel
positions; curves are traced per column by color mask (median of matching pixels), markers by
connected components. Fig. 2 y-axis is unlabeled: the x-axis line is taken as zero, curves normalized to max.

Usage: pixi run python sources/anderson2013_digitized/digitize.py sources/anderson2013_digitized sources/anderson2013_digitized
"""
import sys, numpy as np, pandas as pd
from PIL import Image
from scipy import ndimage
figs, out = sys.argv[1], sys.argv[2]
load = lambda n: np.asarray(Image.open(f"{figs}/{n}.png").convert("RGB")).astype(int)

def lin(p0, v0, p1, v1):
    return lambda p: v0 + (np.asarray(p) - p0) * (v1 - v0) / (p1 - p0)

def mask_color(a, ref, tol):
    return np.sqrt(((a - np.array(ref)) ** 2).sum(-1)) < tol

def trace(a, ref, tol, box, excl=(), fx=None, fy=None, step=1):
    x0, x1, y0, y1 = box
    m = mask_color(a, ref, tol)
    for ex0, ex1, ey0, ey1 in excl:
        m[ey0:ey1, ex0:ex1] = False
    xs, ys = [], []
    for x in range(x0, x1, step):
        col = np.where(m[y0:y1, x:x+step].any(1))[0]
        if len(col) >= 2:
            xs.append(x + (step-1)/2); ys.append(np.median(col) + y0)
    return fx(np.array(xs)), fy(np.array(ys))

def to_df(curves, xname):
    grid = np.unique(np.concatenate([c[0] for c in curves.values()]))
    df = pd.DataFrame({xname: grid})
    for k, (x, y) in curves.items():
        df[k] = np.interp(grid, x, y, left=np.nan, right=np.nan)
    return df

# ---- Fig 4 (both panels): QE in % --------------------------------------------
a = load("fig4")
cols = {"total": (112, 48, 160), "red": (240, 0, 0), "green": (0, 176, 80), "blue": (68, 114, 196)}
panels = {
    "pulsed": dict(box=(112, 962, 5, 625), fx=lin(110.5, 0, 961.5, 175), fy=lin(622, 0, 11.5, 8),
                   excl=[(640, 960, 120, 330)], xname="pulse_energy_density_mJ_cm2"),
    "cw": dict(box=(1161, 1944, 5, 622), fx=lin(1159.5, 0, 1943.5, 100), fy=lin(619.5, 0, 13.5, 8),
               excl=[(1640, 1968, 115, 320)], xname="power_density_W_cm2"),
}
for name, p in panels.items():
    curves = {f"{k}_QE_pct": trace(a, c, 60, p["box"], p["excl"], p["fx"], p["fy"]) for k, c in cols.items()}
    to_df(curves, p["xname"]).round(4).to_csv(f"{out}/fig4_{name}_QE.csv", index=False)

# ---- Fig 2 top-left: UC decays after 943 nm pulse (normalized to max) -------
a = load("fig2")
fx, fy = lin(56.5, 0, 992, 2.0), lin(504.5, 0, 0, 1)  # y: axis line = 0, arbitrary scale
box = (58, 993, 0, 503)
legend = (690, 1000, 0, 290)
cols = {"1um": ((0, 0, 0), 60), "650nm": ((240, 0, 0), 70), "540nm": ((0, 170, 72), 50), "410nm": ((0, 90, 240), 60)}
curves = {}
for k, (c, tol) in cols.items():
    x, y = trace(a, c, tol, box, [legend], fx, fy, step=4)
    curves[k] = (x, y / y.max())
to_df(curves, "time_ms").round(4).to_csv(f"{out}/fig2_decay_943nm.csv", index=False)

# ---- Fig 3 bottom-left: green/red ratio markers -------------------------------
a = load("fig3")
fx, fy = lin(98, 0, 653, 200), lin(1006, 0, 613, 16)
rows = []
for kind, c, tol in [("experimental", (0, 112, 192), 50), ("calculated_new_model", (240, 0, 0), 80)]:
    m = mask_color(a, c, tol)
    m[:, :99] = m[:, 660:] = False; m[:560] = m[1006:] = False
    m[560:730, 340:660] = False  # legend
    lab, n = ndimage.label(m)
    for i, (sy, sx) in enumerate(ndimage.find_objects(lab)):
        if (sy.stop - sy.start) < 8: continue
        cx, cy = (sx.start + sx.stop - 1) / 2, (sy.start + sy.stop - 1) / 2
        rows.append((kind, float(fx(cx)), float(fy(cy))))
# Manual fixes for occluded diamonds (marker height = 19 px, so center = top vertex + 9):
# - the 2.5 mJ/cm2 diamond is split in two by the square drawn over it
# - at 36 mJ/cm2 only the diamond's top vertex (y=942) sticks out of the square
rows = [r for r in rows if not (r[0] == "experimental" and (r[1] < 4 or 30 < r[1] < 40))]
rows += [("experimental", float(fx(104.5)), float(fy(584))), ("experimental", float(fx(197)), float(fy(951)))]
pd.DataFrame(rows, columns=["series", "pulse_energy_density_mJ_cm2", "green_red_ratio"]).sort_values(
    ["series", "pulse_energy_density_mJ_cm2"]).round(3).to_csv(f"{out}/fig3_green_red_ratio.csv", index=False)
print(pd.read_csv(f"{out}/fig3_green_red_ratio.csv").to_string())
