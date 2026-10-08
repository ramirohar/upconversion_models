"""Iteración 4, O11: tiempos de subida y decaimiento tras un pulso contra T.

    pixi run python noche/iter4_dynamics.py
"""

import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))

from dynamics import rise_decay_vs_T  # noqa: E402

from upconversion_models.anderson_thermal import AndersonThermal as M  # noqa: E402

TEMPS = [10, 15, 20, 25, 30, 40, 50, 60, 80, 100, 150, 200, 250, 300, 350, 400]
CONFIGS = {
    "r0.3": {M.r_res: 0.3},
    "r0.1": {M.r_res: 0.1},
    "r0.6": {M.r_res: 0.6},
    "ablation nr_thermal=0": {M.nr_thermal: 0},
}
OUT = pathlib.Path(__file__).parent / "data" / "iter4_dynamics.json"

if __name__ == "__main__":
    res = {}
    for x0 in (0.01, 0.1):
        for name, extra in CONFIGS.items():
            key = f"{name}|x0={x0}"
            res[key] = rise_decay_vs_T(M, TEMPS, x0=x0, extra=extra)
            g = res[key]
            print(key, " ".join(f"{T}:{g[float(T)]['green']['tau_r'] * 1e6:.1f}" for T in TEMPS), flush=True)
    OUT.write_text(json.dumps(res, indent=1))
