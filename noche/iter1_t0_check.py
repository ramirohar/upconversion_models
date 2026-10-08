"""Iteración 1: cuánto se aparta el modelo por defecto (Suta k_nr0, compañeros ascendentes)
de Anderson a T0 = 300 K, en las líneas principales.

    pixi run python noche/iter1_t0_check.py
"""

import json
import pathlib

from pint import get_application_registry

from upconversion_models.anderson_thermal import AndersonThermal as M
from upconversion_models.jablonski_patch import EmissionTransform, Simulator
from upconversion_models.steady import steady_state

u = get_application_registry()
base = json.loads((pathlib.Path(__file__).parent / "data" / "baseline_anderson.json").read_text())
sim = Simulator(M).with_transform(EmissionTransform("emission"), append=True)
for P in ("0.1", "10", "1000"):
    for label, extra in (("defecto", {}), ("sin compañeros ascendentes", {M.detailed_balance: 0})):
        r = steady_state(sim, {M.energy_flux: float(P) * u.W / u.cm**2, M.T: 300 * u.K, **extra})
        ref = base[P]["values"]
        g = (r["line_rad6h1"] + r["line_rad6s1"]) / ref["line_rad61"] - 1
        red = r["line_rad51"] / ref["line_rad51"] - 1
        blue = r["line_rad81"] / ref["line_rad81"] - 1
        print(f"P={P:>5} {label:28s} verde {g:+.3%}  rojo {red:+.3%}  azul {blue:+.3%}")
