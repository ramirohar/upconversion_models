"""Barrido estándar de PLAN_NOCHE.md (T de 10 a 700 K, varias potencias).

Observables fijos por (T, P):
- líneas: verde H (rad6h1, ~525 nm), verde S (rad6s1, ~545 nm), 2H9/2->4I13/2 (rad82,
  ~556 nm, cae dentro de la ventana de integración de 4S3/2 de Zhou 2013:
  534,5-567 nm), rojo (rad51), azul (rad81);
- FIR_pure = rad6h1 / rad6s1 y FIR_band = rad6h1 / (rad6s1 + rad82);
- cociente rojo/verde;
- pendientes log-log d ln I / d ln P (diferencia centrada con P·1.1 y P/1.1);
- tiempos de vida bajo excitación directa débil (autovalores del jacobiano a bombeo
  nulo restringido a {Er6s, Er6h} y a {Er5});
- ΔE efectivo de la pendiente de ln FIR contra 1/T en ventanas de la literatura.
"""

from __future__ import annotations

import json
import pathlib
import time

import numpy as np
from pint import get_application_registry

from upconversion_models.jablonski_patch import EmissionTransform, Simulator
from upconversion_models.steady import _jacobian, steady_state

u = get_application_registry()

TEMPS = np.r_[np.arange(10, 100, 10), np.arange(100, 701, 20)].astype(float)
POWERS = [0.1, 0.8, 2.0, 10.0, 100.0, 290.0]  # W/cm2; 0.8 Xu 2024, 2 Tong 2015, 290 Yu 2014
DP = 1.1
K_CM = 0.6950348  # k_B / (h c) in cm^-1 / K
WINDOWS = {"Zhou 160-300": (160, 300), "Xu 353-453": (353, 453), "Tong 303-573": (303, 573),
           "Geit 300-700": (300, 700)}


def make_sim(model):
    return Simulator(model).with_transform(EmissionTransform("emission"), append=True)


def observables(r):
    green = r["line_rad6h1"] + r["line_rad6s1"]
    return {
        "H": r["line_rad6h1"],
        "S": r["line_rad6s1"],
        "b82": r["line_rad82"],
        "green": green,
        "red": r["line_rad51"],
        "blue": r["line_rad81"],
        "FIR_pure": r["line_rad6h1"] / r["line_rad6s1"],
        "FIR_band": r["line_rad6h1"] / (r["line_rad6s1"] + r["line_rad82"]),
        "RG": r["line_rad51"] / green,
        "Er6s": r["Er6s"],
        "Er6h": r["Er6h"],
        "Er3": r["Er3"],
        "Er2": r["Er2"],
        "Yb2": r["Yb2"],
    }


def lifetimes(sim, model, T, extra):
    prob = sim.with_values({model.T: T * u.K, model.energy_flux: 0 * u.W / u.cm**2, **extra}).create_problem()
    names = [str(v) for v in sim.compiled.variables]
    dy = np.empty_like(prob.y)
    f = lambda y: np.array(prob.rhs(0.0, y, prob.p, dy)).copy()
    J = _jacobian(f, prob.y, floor=1e-6)
    i6 = [names.index("Er6s"), names.index("Er6h")]
    ev = np.linalg.eigvals(J[np.ix_(i6, i6)])
    i5 = names.index("Er5")
    return {"tau_green": float(1 / np.min(np.abs(ev))), "tau_red": float(1 / abs(J[i5, i5]))}


def run(model, extra=None, temps=TEMPS, powers=POWERS, out=None):
    extra = extra or {}
    sim = make_sim(model)
    t0 = time.time()
    data = {"temps": list(map(float, temps)), "powers": list(powers), "obs": {}, "slopes": {}, "tau": {}}
    for T in temps:
        data["tau"][str(T)] = lifetimes(sim, model, T, extra)
        for P in powers:
            vals = {}
            for tag, fac in (("c", 1.0), ("hi", DP), ("lo", 1 / DP)):
                r = steady_state(sim, {model.T: T * u.K, model.energy_flux: P * fac * u.W / u.cm**2, **extra})
                vals[tag] = observables(r)
            key = f"{T}|{P}"
            data["obs"][key] = vals["c"]
            data["slopes"][key] = {
                k: float(np.log(vals["hi"][k] / vals["lo"][k]) / np.log(DP**2))
                for k in ("H", "S", "green", "red", "blue")
            }
    data["seconds"] = time.time() - t0
    data["dE_eff"] = fits(data)
    if out:
        pathlib.Path(out).write_text(json.dumps(data, indent=1))
    return data


def series(data, name, P):
    return np.array([data["obs"][f"{T}|{P}"][name] for T in data["temps"]])


def fits(data):
    """ΔE efectivo (cm^-1) y prefactor del ajuste ln FIR = ln A - ΔE/(kT) en cada ventana."""
    T = np.array(data["temps"])
    res = {}
    for P in data["powers"]:
        for fir in ("FIR_pure", "FIR_band"):
            y = np.log(series(data, fir, P))
            for wname, (a, b) in WINDOWS.items():
                m = (T >= a) & (T <= b)
                slope, icpt = np.polyfit(1 / T[m], y[m], 1)
                res[f"{P}|{fir}|{wname}"] = {"dE": float(-slope * K_CM), "A": float(np.exp(icpt))}
    return res
