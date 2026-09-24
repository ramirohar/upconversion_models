from jablonski import Parameter, SpectroscopicSystem, SingletState, assign, initial
from jablonski.transitions import MassAction
from upconversion_models.utils import k, h, c
import numpy as np

import pint

u = pint.get_application_registry()

class TemperatureDependentIC(SpectroscopicSystem):
    source: SingletState = initial(0.0, default=0)
    target: SingletState = initial(0.0, default=0)

    reference_rate: Parameter = assign(default=0 / u.s)   
    T: Parameter = assign(default=300 * u.K)
    T_A: Parameter = assign(default=300 * u.K)             
    phonon_wavenumber: Parameter = assign(default=500 / u.cm)  

    phonon_E = h * c * phonon_wavenumber   

    dE = source.energy - target.energy     
    p = np.ceil(dE / phonon_E)             

    # [1 + n(T)] / [1 + n(T_A)], written stably
    ratio = np.expm1(-phonon_E / (k * T_A)) / np.expm1(-phonon_E / (k * T))

    ic = MassAction(
        reactants=[source],
        products=[target],
        rate=reference_rate * ratio ** p,
    )

class EnergyTransferUpconversion(SpectroscopicSystem):
    sensitizer_high: SingletState = initial(0.0, default=0)
    sensitizer_low: SingletState = initial(0.0, default=0)

    activator_high: SingletState = initial(0.0, default=0)
    activator_low: SingletState = initial(0.0, default=0)

    rate: Parameter = assign(default=0 / u.s)

    upconversion = MassAction(
        reactants=[sensitizer_high, activator_low],
        products=[sensitizer_low, activator_high],
        rate=rate,
    )