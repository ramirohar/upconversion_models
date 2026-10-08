"""Tabla de τ_R y τ_D (µs) tras un pulso, por banda y configuración (salida de iter4_dynamics.py).

    pixi run python noche/dyn_table.py data/iter4_dynamics.json
"""

import json
import pathlib
import sys

d = json.loads((pathlib.Path(__file__).parent / sys.argv[1]).read_text())
for key, res in d.items():
    temps = sorted(res, key=float)
    print(f"== {key}")
    for band in ("green", "red", "blue"):
        for q in ("tau_r", "tau_d"):
            print(f"  {band:5s} {q}: " + " ".join(f"{float(T):.0f}:{res[T][band][q] * 1e6:.0f}" for T in temps))
