"""Iteración 1: barrido estándar de AndersonThermal (parámetros por defecto).

    pixi run python noche/iter1_run.py
"""

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))

from sweep import run  # noqa: E402

from upconversion_models.anderson_thermal import AndersonThermal  # noqa: E402

OUT = pathlib.Path(__file__).parent / "data" / "iter1_sweep.json"

if __name__ == "__main__":
    data = run(AndersonThermal, out=OUT)
    print(f"{data['seconds']:.0f} s")
