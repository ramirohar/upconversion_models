"""Iteración 1: robustez del máximo del rojo contra T a alta potencia (no pre-registrado).

Posición del máximo de rojo, verde y azul a 100 y 290 W/cm² para ħω = 300, 359 y 450 cm⁻¹ y
con/sin compañeros ascendentes; y la ablación global thermal = 0.

    pixi run python noche/iter1_redmax.py
"""

import numpy as np
from pint import get_application_registry

from upconversion_models.anderson_thermal import AndersonThermal as M
from upconversion_models.jablonski_patch import EmissionTransform, Simulator
from upconversion_models.steady import steady_state

u = get_application_registry()
sim = Simulator(M).with_transform(EmissionTransform("emission"), append=True)
TEMPS = np.arange(20, 501, 20.0)

for P in (100.0, 290.0):
    for hw in (300, 359, 450):
        for db in (1, 0):
            rows = [
                steady_state(sim, {M.energy_flux: P * u.W / u.cm**2, M.T: T * u.K,
                                   M.phonon_wavenumber: hw / u.cm, M.detailed_balance: db})
                for T in TEMPS
            ]
            red = np.array([r["line_rad51"] for r in rows])
            green = np.array([r["line_rad6h1"] + r["line_rad6s1"] for r in rows])
            blue = np.array([r["line_rad81"] for r in rows])
            i300 = int(np.argmin(np.abs(TEMPS - 300)))
            print(f"P={P:5.0f} ħω={hw} db={db}: T_max rojo {TEMPS[red.argmax()]:.0f} K "
                  f"(máx/300K = {red.max() / red[i300]:.2f}, 20K/máx = {red[0] / red.max():.2f}); "
                  f"T_max verde {TEMPS[green.argmax()]:.0f} K; T_max azul {TEMPS[blue.argmax()]:.0f} K")
