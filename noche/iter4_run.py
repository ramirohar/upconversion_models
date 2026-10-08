"""Iteración 4: barrido estándar con 1<->3 y 3<->7 asistidas por un fonón (Li 2014) para
r = 0, 0.1, 0.3 y 0.6 (ħω = 359 cm⁻¹). La ablación (nr_thermal = 0) es iter3_sweep.

    pixi run python noche/iter4_run.py 0.3      # un valor de r por proceso
"""

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))

from sweep import run  # noqa: E402

from upconversion_models.anderson_thermal import AndersonThermal as M  # noqa: E402

DATA = pathlib.Path(__file__).parent / "data"

if __name__ == "__main__":
    for r in map(float, sys.argv[1:] or ["0", "0.1", "0.3", "0.6"]):
        d = run(M, extra={M.r_res: r}, out=DATA / f"iter4_r{r:g}.json")
        print(f"r={r:g}", f"{d['seconds']:.0f} s", flush=True)
