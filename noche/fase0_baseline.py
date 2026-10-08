"""Fase 0: línea base de AndersonModel en estado estacionario.

Guarda poblaciones y líneas de emisión (tasas `line_*`, en s^-1 por sitio) para una
grilla de potencias, calculadas con `upconversion_models.steady.steady_state`
(BDF rtol=1e-10, atol=1e-20 + pulido de Newton). Uso:

    pixi run python noche/fase0_baseline.py
"""

import json
import pathlib
import time

import numpy as np
from pint import get_application_registry

from upconversion_models import AndersonModel
from upconversion_models.jablonski_patch import EmissionTransform, Simulator
from upconversion_models.steady import steady_state

u = get_application_registry()

POWERS = [0.1, 0.3, 1, 3, 10, 30, 100, 300, 1000]  # W/cm2
OUT = pathlib.Path(__file__).parent / "data" / "baseline_anderson.json"


def simulator():
    return Simulator(AndersonModel).with_transform(
        EmissionTransform("emission"), append=True
    )


def compute(powers=POWERS):
    sim = simulator()
    rows = {}
    for P in powers:
        t = time.time()
        r, info = steady_state(
            sim, {AndersonModel.energy_flux: P * u.W / u.cm**2}, return_info=True
        )
        rows[str(P)] = {"values": r, "info": info, "seconds": time.time() - t}
    return rows


if __name__ == "__main__":
    rows = compute()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(rows, indent=1))
    for P, row in rows.items():
        v = row["values"]
        print(
            f"P={P:>6} W/cm2  t={row['seconds']:.2f}s  newton={row['info']['newton_iter']}"
            f"  G/R={v['line_rad61'] / v['line_rad51']:.3f}"
        )
