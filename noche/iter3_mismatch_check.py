"""Robustez del signo de R/G(T) de la iteración 3 frente al desajuste supuesto del ETU 6→9.

El estado final de Er9 (26100 cm⁻¹) representa el manifold 4G,2K (Anderson 2013: "into the 4G,2K
manifold, above 2H9/2"); con un estado final más alto el desajuste baja. Se prueban desajustes de
1000/1650 cm⁻¹ (6s/6h) además de los de baricentros (2400/3050).

    pixi run python noche/iter3_mismatch_check.py
"""

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))

from pint import get_application_registry  # noqa: E402
from sweep import make_sim, observables  # noqa: E402

from upconversion_models.anderson_thermal import AndersonThermal as M  # noqa: E402
from upconversion_models.steady import steady_state  # noqa: E402

u = get_application_registry()
sim = make_sim(M)
for label, extra in (
    ("baricentros 2400/3050", {}),
    ("reducido 1000/1650", {M.etu6s9._mismatch_wavenumber: 1000 / u.cm, M.etu6h9._mismatch_wavenumber: 1650 / u.cm}),
):
    for P in (0.8, 290.0):
        rg = {T: observables(steady_state(sim, {M.T: T * u.K, M.energy_flux: P * u.W / u.cm**2, **extra}))["RG"]
              for T in (10.0, 303.0, 400.0, 483.0)}
        print(f"{label:22s} P={P:5g}: R/G 400K/10K = {rg[400.0] / rg[10.0]:.2f}, 483K/303K = {rg[483.0] / rg[303.0]:.2f}")
