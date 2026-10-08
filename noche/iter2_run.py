"""Iteración 2: barrido estándar con CR6 asistida por fonones (ħω = 359 cm⁻¹) y variantes
de ħω = 300 y 450 cm⁻¹ (rango de literatura) para la robustez.

    pixi run python noche/iter2_run.py
"""

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))

from pint import get_application_registry  # noqa: E402
from sweep import run  # noqa: E402

from upconversion_models.anderson_thermal import AndersonThermal as M  # noqa: E402

u = get_application_registry()
DATA = pathlib.Path(__file__).parent / "data"

if __name__ == "__main__":
    for name, extra in (
        ("iter2_sweep", {}),
        ("iter2_hw300", {M.phonon_wavenumber: 300 / u.cm}),
        ("iter2_hw450", {M.phonon_wavenumber: 450 / u.cm}),
        ("iter1_hw300", {M.phonon_wavenumber: 300 / u.cm, M.cr_thermal: 0}),
        ("iter1_hw450", {M.phonon_wavenumber: 450 / u.cm, M.cr_thermal: 0}),
    ):
        d = run(M, extra=extra, out=DATA / f"{name}.json")
        print(name, f"{d['seconds']:.0f} s", flush=True)
