"""Iteración 5, exploración **post hoc** (no cuenta para el criterio de emergencia):
decaimiento al apagar la excitación continua de 0,8 W/cm² (la potencia de los espectros de Xu 2024)
para el modelo final con ħω = 300, 359 y 450 cm⁻¹, y la ablación thermal = 0.

    pixi run python noche/iter5_cwoff.py
"""

import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))

from dynamics import cw_off_curves, fit_single_exp  # noqa: E402
from pint import get_application_registry  # noqa: E402

from upconversion_models.anderson_thermal import AndersonThermal as M  # noqa: E402

u = get_application_registry()
OUT = pathlib.Path(__file__).parent / "data" / "iter5_cwoff.json"

if __name__ == "__main__":
    res = {}
    for name, extra in (
        ("final hw=300", {M.phonon_wavenumber: 300 / u.cm}),
        ("final hw=359", {}),
        ("final hw=450", {M.phonon_wavenumber: 450 / u.cm}),
        ("ablación thermal=0", {M.thermal: 0}),
    ):
        row = {}
        for T in (55.0, 230.0, 353.0, 453.0):
            t, b = cw_off_curves(M, T, 0.8, extra=extra)
            row[str(T)] = {k: fit_single_exp(t, b[k]) for k in ("red", "541", "520")}
        res[name] = row
        f = lambda T, k: row[str(T)][k] * 1e6
        print(f"{name:20s} rojo 55/230 K: {f(55.0, 'red'):.0f}/{f(230.0, 'red'):.0f} µs "
              f"(×{f(55.0, 'red') / f(230.0, 'red'):.2f}) | 541 353/453 K: {f(353.0, '541'):.0f}/{f(453.0, '541'):.0f} µs "
              f"(×{f(353.0, '541') / f(453.0, '541'):.2f}) | 520 353/453 K: {f(353.0, '520'):.0f}/{f(453.0, '520'):.0f} µs",
              flush=True)
    OUT.write_text(json.dumps(res, indent=1))
