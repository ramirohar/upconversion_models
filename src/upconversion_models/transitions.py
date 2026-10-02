from jablonski import Parameter, SpectroscopicSystem, SingletState, assign, initial
from jablonski.transitions import MassAction
from poincare import Constant
from upconversion_models.utils import k, h, c
import numpy as np
from symbolite.abstract import real
import pint

u = pint.get_application_registry()


class Base(SpectroscopicSystem):
    T: Parameter = assign(default=300 * u.K)
    source: SingletState = initial(0.0, default=0)
    target: SingletState = initial(0.0, default=0)
    phonon_wavenumber: Parameter = assign(default=500 / u.cm)
    _gap_wavenumber: Parameter = assign(default=0 / u.cm)
    f: Parameter = assign(default=h * c / k)  # hc / k_B
    p = (_gap_wavenumber / phonon_wavenumber)

    def __init__(self, *args, **kwargs):
        src, tgt = kwargs.get("source"), kwargs.get("target")
        if "_gap_wavenumber" not in kwargs and src is not None and tgt is not None:
            with u.context("sp"):
                kwargs["_gap_wavenumber"] = abs(src.energy - tgt.energy).to(1 / u.cm)
        super().__init__(*args, **kwargs)


class TemperatureDependentInternalRelaxationRef(Base):
    """Downward phonon emission: k(T) = k(T_A) · [(1 + n(T)) / (1 + n(T_A))]^p"""

    reference_rate: Parameter = assign(default=0 / u.s)
    T_A: Parameter = assign(default=300 * u.K)

    factor = real.expm1(-Base.f * Base.phonon_wavenumber / T_A) / real.expm1(
        -Base.f * Base.phonon_wavenumber / Base.T
    )

    ic = MassAction(
        reactants=[Base.source],
        products=[Base.target],
        rate=reference_rate * factor**Base.p,
    )


class TemperatureDependentInternalExcitationRef(Base):
    """Upward phonon absorption: k(T) = k(T_A) · [n(T) / n(T_A)]^p"""

    reference_rate: Parameter = assign(default=0 / u.s)
    T_A: Parameter = assign(default=300 * u.K)

    factor = real.expm1(Base.f * Base.phonon_wavenumber / T_A) / real.expm1(
        Base.f * Base.phonon_wavenumber / Base.T
    )

    ic = MassAction(
        reactants=[Base.source],
        products=[Base.target],
        rate=reference_rate * factor**Base.p,
    )


class TemperatureDependentInternalRelaxation(Base):
    """Downward phonon emission: k(T) = rate · (1 + n(T))^p"""

    rate: Parameter = assign(default=0 / u.s)

    factor = 1 + 1 / real.expm1(Base.f * Base.phonon_wavenumber / Base.T)

    ic = MassAction(
        reactants=[Base.source],
        products=[Base.target],
        rate=rate * factor**Base.p,
    )

class TemperatureDependentInternalExcitation(Base):
    """Upward phonon absorption: k(T) = rate · n(T)^p"""

    rate: Parameter = assign(default=0 / u.s)

    factor = 1 / real.expm1(Base.f * Base.phonon_wavenumber / Base.T)

    ic = MassAction(
        reactants=[Base.source],
        products=[Base.target],
        rate=rate * factor**Base.p,
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
