"""Iteración 3: barrido estándar con transferencias Yb<->Er asistidas por fonones
(ħω = 359 cm⁻¹) y variantes ħω = 300 y 450 cm⁻¹. La ablación (et_thermal = 0) son los
barridos de la iteración 2 (iter2_sweep, iter2_hw300, iter2_hw450).

    pixi run python noche/iter3_run.py
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
        ("iter3_sweep", {}),
        ("iter3_hw300", {M.phonon_wavenumber: 300 / u.cm}),
        ("iter3_hw450", {M.phonon_wavenumber: 450 / u.cm}),
    ):
        d = run(M, extra=extra, out=DATA / f"{name}.json")
        print(name, f"{d['seconds']:.0f} s", flush=True)
