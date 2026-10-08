"""Iteración 5: decaimiento de la UC tras un pulso de 980 nm contra T (Xu 2024, Fig. 5),
con el protocolo pre-registrado en CUADERNO_NOCHE.md (commit 40ce23e).

    pixi run python noche/iter5_decays.py
"""

import itertools
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))

import numpy as np  # noqa: E402
from dynamics import fit_single_exp, pulse_curves, rise_decay  # noqa: E402
from pint import get_application_registry  # noqa: E402

from upconversion_models.anderson_thermal import AndersonThermal as M  # noqa: E402

u = get_application_registry()
OUT = pathlib.Path(__file__).parent / "data" / "iter5_decays.json"
T0 = 300.0  # K; the block parameter is stored as a plain magnitude

CONFIGS = {
    "final (iter. 3)": {},
    "iter. 2": {M.et_thermal: 0},
    "iter. 1": {M.et_thermal: 0, M.cr_thermal: 0},
    "ablación thermal=0": {M.thermal: 0},
    "sin T en norad3": {M.norad3.T: T0},
    "sin T en norad4": {M.norad4.T: T0},
    "sin T en norad8": {M.norad8.T: T0},
    "sin T en norad9": {M.norad9.T: T0},
}
HW = [300, 359, 450]
X0 = [0.003, 0.01, 0.1]
TEMPS = {"red": (55.0, 230.0), "541": (353.0, 453.0), "520": (353.0, 453.0)}


def taus(t, I):
    tau_r, tau_d, tp = rise_decay(t, I)
    return {"slope": tau_d, "xu_fit": fit_single_exp(t, I)}


if __name__ == "__main__":
    res = {}
    for (cname, extra), hw, x0 in itertools.product(CONFIGS.items(), HW, X0):
        if cname.startswith("sin T") and (hw != 359 or x0 != 0.01):
            continue  # ablaciones por tasa solo en el punto central
        key = f"{cname}|hw={hw}|x0={x0}"
        row = {}
        for T in sorted({T for pair in TEMPS.values() for T in pair}):
            t, bands = pulse_curves(M, T, x0=x0, extra={**extra, M.phonon_wavenumber: hw / u.cm})
            row[str(T)] = {b: taus(t, bands[b]) for b in TEMPS}
        res[key] = row
        r = row
        red = [r[str(T)]["red"]["xu_fit"] * 1e6 for T in TEMPS["red"]]
        g = [r[str(T)]["541"]["xu_fit"] * 1e6 for T in TEMPS["541"]]
        print(f"{key:40s} rojo 55/230 K: {red[0]:.0f}/{red[1]:.0f} µs (×{red[0] / red[1]:.2f})  "
              f"541 353/453 K: {g[0]:.0f}/{g[1]:.0f} µs (×{g[0] / g[1]:.2f})", flush=True)
    OUT.write_text(json.dumps(res, indent=1))
