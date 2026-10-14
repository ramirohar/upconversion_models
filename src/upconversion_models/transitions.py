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
    p = _gap_wavenumber / phonon_wavenumber

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


class TemperatureDependentInternalConversionRef(Base):
    """Downward phonon emission plus its upward partner by detailed balance:
    k_down(T) = k(T_A) · [(1 + n(T)) / (1 + n(T_A))]^p
    k_up(T) = k(T_A) · (g_source / g_target) · [n(T) / (1 + n(T_A))]^p
    so that k_up / k_down = (g_source / g_target) · exp(-ΔE / kT) at every T.
    `detailed_balance` (default 1) scales k_up: 0 removes the upward partner,
    e.g. to check the reduction to a model without it."""

    reference_rate: Parameter = assign(default=0 / u.s)
    T_A: Parameter = assign(default=300 * u.K)
    g_source: Parameter = assign(default=1)
    g_target: Parameter = assign(default=1)
    detailed_balance: Parameter = assign(default=1)

    factor_down = real.expm1(-Base.f * Base.phonon_wavenumber / T_A) / real.expm1(
        -Base.f * Base.phonon_wavenumber / Base.T
    )
    factor_up = -real.expm1(-Base.f * Base.phonon_wavenumber / T_A) / real.expm1(
        Base.f * Base.phonon_wavenumber / Base.T
    )

    down = MassAction(
        reactants=[Base.source],
        products=[Base.target],
        rate=reference_rate * factor_down**Base.p,
    )
    up = MassAction(
        reactants=[Base.target],
        products=[Base.source],
        rate=detailed_balance * reference_rate * g_source / g_target * factor_up**Base.p,
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


class PhononAssistedEnergyTransfer(SpectroscopicSystem):
    """Energy transfer whose mismatch ΔE = (E_sh - E_sl) - (E_ah - E_al) is
    compensated by p phonons of energy ε (Miyakawa–Dexter; Auzel 2003, §2.1):
    |ΔE| >= ħω: multiphonon, ε = ħω and p = |ΔE| / ħω
    |ΔE| <  ħω: a single phonon of the mismatch energy, ε = |ΔE| and p = 1
    With n(T) the occupation of phonons of energy ε,
    ΔE > 0 (emission):   k(T) = k(T_A) · [(1 + n(T)) / (1 + n(T_A))]^p
    ΔE < 0 (absorption): k(T) = k(T_A) · [r + (1 - r) · (n(T) / n(T_A))^p]
    r is a resonant fraction (Stark components, inhomogeneous broadening)
    that keeps endothermic transfers alive as T -> 0.
    ΔE is computed from the state energies unless `mismatch` (signed) is given,
    e.g. a value measured between Stark components instead of barycenters.

    When the transfer is resonant only from the upper crystal-field component |1>
    of the sensitizer, Δ01 above |0> (Suyver 2005; Yu 2014), the rate is further
    weighted by the thermal population of |1>, p1(T) = 1 / (1 + exp(Δ01 / kT)):
    k(T) -> k(T) · [η + (1 - η) · p1(T)] / [η + (1 - η) · p1(T_A)]
    with η the relative efficiency of the transfer from |0>. Δ01 = 0 or η = 1
    leave the rate unchanged."""

    sensitizer_high: SingletState = initial(0.0, default=0)
    sensitizer_low: SingletState = initial(0.0, default=0)

    activator_high: SingletState = initial(0.0, default=0)
    activator_low: SingletState = initial(0.0, default=0)

    rate: Parameter = assign(default=0 / u.s)
    T: Parameter = assign(default=300 * u.K)
    T_A: Parameter = assign(default=300 * u.K)
    phonon_wavenumber: Parameter = assign(default=350 / u.cm)
    resonant_fraction: Parameter = assign(default=0)
    _mismatch_wavenumber: Parameter = assign(default=0 / u.cm)
    _emits: Parameter = assign(default=1)
    _one_phonon: Parameter = assign(default=0)
    f: Parameter = assign(default=h * c / k)  # hc / k_B
    phonon = _one_phonon * _mismatch_wavenumber + (1 - _one_phonon) * phonon_wavenumber
    p = _one_phonon + (1 - _one_phonon) * _mismatch_wavenumber / phonon_wavenumber

    emission = real.expm1(-f * phonon / T_A) / real.expm1(-f * phonon / T)
    absorption = real.expm1(f * phonon / T_A) / real.expm1(f * phonon / T)
    factor = _emits * emission**p + (1 - _emits) * (
        resonant_fraction + (1 - resonant_fraction) * absorption**p
    )

    sensitizer_stark_gap: Parameter = assign(default=0 / u.cm)
    sensitizer_ground_efficiency: Parameter = assign(default=1)
    p1 = 1 / (1 + real.exp(f * sensitizer_stark_gap / T))
    p1_A = 1 / (1 + real.exp(f * sensitizer_stark_gap / T_A))
    stark = (sensitizer_ground_efficiency + (1 - sensitizer_ground_efficiency) * p1) / (
        sensitizer_ground_efficiency + (1 - sensitizer_ground_efficiency) * p1_A
    )

    upconversion = MassAction(
        reactants=[sensitizer_high, activator_low],
        products=[sensitizer_low, activator_high],
        rate=rate * factor * stark,
    )

    def __init__(self, *args, mismatch=None, **kwargs):
        states = [
            kwargs.get(s)
            for s in (
                "sensitizer_high",
                "sensitizer_low",
                "activator_high",
                "activator_low",
            )
        ]
        if mismatch is None and all(s is not None for s in states):
            sh, sl, ah, al = states
            # ΔE / hc directly: the "sp" context goes through 1/E and fails at ΔE = 0
            mismatch = (
                ((sh.energy - sl.energy) - (ah.energy - al.energy)) / (h * c)
            ).to(1 / u.cm)
        if mismatch is not None and "_mismatch_wavenumber" not in kwargs:
            hw = kwargs.get("phonon_wavenumber", 350 / u.cm)
            hw = getattr(hw, "default", hw).to(1 / u.cm).magnitude
            gap = abs(mismatch.to(1 / u.cm).magnitude)
            kwargs["_mismatch_wavenumber"] = gap / u.cm
            kwargs["_emits"] = 1 if mismatch.magnitude >= 0 else 0
            kwargs["_one_phonon"] = 1 if 0 < gap < hw else 0
        super().__init__(*args, **kwargs)
