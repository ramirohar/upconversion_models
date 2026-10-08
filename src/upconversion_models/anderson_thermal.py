"""AndersonModel with the green pair split and temperature-dependent multiphonon rates.

Iteración 1 of PLAN_NOCHE.md. Starting from `AndersonModel` (Anderson 2013, Table 1):

- Er6 (2H11/2 + 4S3/2) is split into Er6s (4S3/2, g = 4, 18300 cm-1) and Er6h
  (2H11/2, g = 12, 650 cm-1 above; Suta 2025, Fig. 1a). They are coupled by Suta's
  two-phonon process (eq. 6): down = g_S·k_nr0·(1 + n)^p, up = g_H·k_nr0·n^p, with
  n the occupation of phonons of energy ΔE/p (p = 2), so up/down = (g_H/g_S)e^(-ΔE/kT).
- Radiative rates: A(2H11/2)/A(4S3/2) = C = 5.52/3 (Suta 2025, LIR prefactor), with
  f_S(T0)·k_r6s + f_H(T0)·k_r6h = Anderson's k_r6 at T0, f_S, f_H the Boltzmann
  fractions. Branching ratios as in Anderson. A(T) is constant.
- Every other Er6 process (k_nr6, k_ET6-9, k_CR6) acts with Anderson's constant on
  both sublevels, so the sum Er6s + Er6h obeys Anderson's equation whenever the
  radiative split averages to k_r6.
- Multiphonon relaxations k_nrX(T) = k_nrX(T0)·[(1 + n(T))/(1 + n(T0))]^p,
  p = gap/ħω (TemperatureDependentInternalConversionRef), with the upward partner by
  detailed balance. Er7 relaxes into Er6h (the closest level).

Switches (same model, no separate class):
- `thermal` = 1 (default) uses T; `thermal` = 0 evaluates every temperature-dependent
  rate at T0, so the model reduces to Anderson at any T (with fast thermalization).
- `k_nr0` sets the 6s <-> 6h coupling; a large value gives instantaneous
  thermalization.
"""

from jablonski import Parameter, SingletState, SpectroscopicSystem, assign, initial
from jablonski.transitions import Absorption, Fluorescence, InternalConversion
from pint import get_application_registry
from symbolite.abstract import real

from .anderson import photon_energy
from .transitions import (
    EnergyTransferUpconversion,
    TemperatureDependentInternalConversionRef,
    TemperatureDependentInternalExcitation,
    TemperatureDependentInternalRelaxation,
)
from .utils import hc, k

u = get_application_registry()

E6S = 18300  # cm-1, Anderson's Er6
GAP_HS = 650  # cm-1, Suta 2025 Fig. 1a
G_S, G_H = 4, 12
C_HS = 5.52 / 3  # A_H / A_S, Suta 2025


class AndersonThermal(SpectroscopicSystem):
    # ---- temperature ----
    T0: Parameter = assign(default=300 * u.K)
    T: Parameter = assign(default=300 * u.K)
    thermal: Parameter = assign(default=1)
    # 1: each multiphonon relaxation has its upward partner by detailed balance;
    # 0: no upward partners (Anderson's equations, used for the reduction tests)
    detailed_balance: Parameter = assign(default=1)
    T_eff = T0 + thermal * (T - T0)

    # Effective phonon energy for multiphonon relaxation: highest lattice Raman mode
    # of β-NaYF4 (Dubey 2023: 253, 307, 359 cm-1; Suta 2025 cutoff 450 cm-1)
    phonon_wavenumber: Parameter = assign(default=359 / u.cm)
    # 2H11/2 <-> 4S3/2 coupling, Suta 2025 eq. 6: p = 2 phonons of ΔE/2
    k_nr0: Parameter = assign(default=2.29e6 / u.s)
    thermalization_phonon: Parameter = assign(default=GAP_HS / 2 / u.cm)

    # ---- Anderson 2013, Table 1 (new UC model) ----
    yb_cross_section: Parameter = assign(default=1e-20 * u.cm**2)  # no value reported
    yb_rad: Parameter = assign(default=393 / u.s)
    yb_nr: Parameter = assign(default=220 / u.s)

    k_13: Parameter = assign(default=5.9 * 2e-16 / u.s * u.cm**3)
    k_25: Parameter = assign(default=0 / u.s * u.cm**3)
    k_37: Parameter = assign(default=1.54e-15 / u.s * u.cm**3)
    k_58: Parameter = assign(default=1.76e-15 / u.s * u.cm**3)
    k_69: Parameter = assign(default=6.07e-15 / u.s * u.cm**3)

    k_31: Parameter = assign(default=2e-16 / u.s * u.cm**3)
    k_52: Parameter = assign(default=0 / u.s * u.cm**3)
    k_73: Parameter = assign(default=2.04e-16 / u.s * u.cm**3)
    k_95: Parameter = assign(default=2.84e-16 / u.s * u.cm**3)

    k_nr9: Parameter = assign(default=1.76e6 / u.s)
    k_nr8: Parameter = assign(default=43_450 / u.s)
    k_nr7: Parameter = assign(default=1e6 / u.s)
    k_nr6: Parameter = assign(default=26 / u.s)
    k_nr5: Parameter = assign(default=0 / u.s)
    k_nr4: Parameter = assign(default=22_120 / u.s)
    k_nr3: Parameter = assign(default=61 / u.s)

    k_r8: Parameter = assign(default=2_330 / u.s)
    k_r6: Parameter = assign(default=1_510 / u.s)
    k_r5: Parameter = assign(default=2_039 / u.s)
    k_r3: Parameter = assign(default=73 / u.s)
    k_r2: Parameter = assign(default=110 / u.s)

    k_cr6: Parameter = assign(default=2.79e-17 / u.s * u.cm**3)
    k_cr4: Parameter = assign(default=8.04e-19 / u.s * u.cm**3)
    k_uc2: Parameter = assign(default=2.31e-17 / u.s * u.cm**3)

    site_density: Parameter = assign(default=1.38e22 / u.cm**3)
    isolated_percentage: Parameter = assign(default=0.05)
    yb_concetration: Parameter = assign(default=0.18)
    er_concetration: Parameter = assign(default=0.02)

    energy_flux: Parameter = assign(default=0 * u.W / u.cm**2)
    photon_flux = energy_flux / photon_energy

    # ---- radiative split of k_r6 at T0 ----
    hs_ratio: Parameter = assign(default=C_HS)
    f: Parameter = assign(default=hc / k)  # hc / k_B
    rho0 = (G_H / G_S) * real.exp(-f * (GAP_HS / u.cm) / T0)  # n_H / n_S at T0
    k_r6s = k_r6 * (1 + rho0) / (1 + hs_ratio * rho0)
    k_r6h = hs_ratio * k_r6s

    # ---- levels (Anderson's barycenters; Er6h from Suta 2025) ----
    Yb1: SingletState = initial(
        energy=0 * hc / u.cm, default=yb_concetration * (1 - isolated_percentage)
    )
    Yb2: SingletState = initial(energy=10200 * hc / u.cm, default=0)
    Er1: SingletState = initial(
        energy=0 * hc / u.cm, default=er_concetration * (1 - isolated_percentage)
    )
    Er2: SingletState = initial(energy=6500 * hc / u.cm, default=0)
    Er3: SingletState = initial(energy=10200 * hc / u.cm, default=0)
    Er4: SingletState = initial(energy=12500 * hc / u.cm, default=0)
    Er5: SingletState = initial(energy=15000 * hc / u.cm, default=0)
    Er6s: SingletState = initial(energy=E6S * hc / u.cm, default=0)
    Er6h: SingletState = initial(energy=(E6S + GAP_HS) * hc / u.cm, default=0)
    Er7: SingletState = initial(energy=20500 * hc / u.cm, default=0)
    Er8: SingletState = initial(energy=24500 * hc / u.cm, default=0)
    Er9: SingletState = initial(energy=26100 * hc / u.cm, default=0)

    # ---- Yb ----
    abs = Absorption(ground=Yb1, excited=Yb2, rate=yb_cross_section, pump=photon_flux)
    rad = Fluorescence(ground=Yb1, excited=Yb2, rate=yb_rad)
    norad = InternalConversion(low=Yb1, high=Yb2, rate=yb_nr)

    # ---- Yb -> Er ETU ----
    etu13 = EnergyTransferUpconversion(
        sensitizer_high=Yb2, sensitizer_low=Yb1, activator_low=Er1, activator_high=Er3,
        rate=k_13 * site_density,
    )
    etu25 = EnergyTransferUpconversion(
        sensitizer_high=Yb2, sensitizer_low=Yb1, activator_low=Er2, activator_high=Er5,
        rate=k_25 * site_density,
    )
    etu37 = EnergyTransferUpconversion(
        sensitizer_high=Yb2, sensitizer_low=Yb1, activator_low=Er3, activator_high=Er7,
        rate=k_37 * site_density,
    )
    etu58 = EnergyTransferUpconversion(
        sensitizer_high=Yb2, sensitizer_low=Yb1, activator_low=Er5, activator_high=Er8,
        rate=k_58 * site_density,
    )
    etu6s9 = EnergyTransferUpconversion(
        sensitizer_high=Yb2, sensitizer_low=Yb1, activator_low=Er6s, activator_high=Er9,
        rate=k_69 * site_density,
    )
    etu6h9 = EnergyTransferUpconversion(
        sensitizer_high=Yb2, sensitizer_low=Yb1, activator_low=Er6h, activator_high=Er9,
        rate=k_69 * site_density,
    )

    # ---- Er -> Yb back transfer ----
    etu31 = EnergyTransferUpconversion(
        sensitizer_high=Er3, sensitizer_low=Er1, activator_low=Yb1, activator_high=Yb2,
        rate=k_31 * site_density,
    )
    etu52 = EnergyTransferUpconversion(
        sensitizer_high=Er5, sensitizer_low=Er2, activator_low=Yb1, activator_high=Yb2,
        rate=k_52 * site_density,
    )
    etu73 = EnergyTransferUpconversion(
        sensitizer_high=Er7, sensitizer_low=Er3, activator_low=Yb1, activator_high=Yb2,
        rate=k_73 * site_density,
    )
    etu95 = EnergyTransferUpconversion(
        sensitizer_high=Er9, sensitizer_low=Er5, activator_low=Yb1, activator_high=Yb2,
        rate=k_95 * site_density,
    )

    # ---- 2H11/2 <-> 4S3/2 thermalization (Suta 2025, eq. 6) ----
    therm_down = TemperatureDependentInternalRelaxation(
        T=T_eff,
        source=Er6h,
        target=Er6s,
        rate=G_S * k_nr0,
        phonon_wavenumber=thermalization_phonon,
    )
    therm_up = TemperatureDependentInternalExcitation(
        T=T_eff,
        source=Er6s,
        target=Er6h,
        rate=G_H * k_nr0,
        phonon_wavenumber=thermalization_phonon,
    )

    # ---- multiphonon relaxation anchored at Anderson's values at T0 ----
    # g = 2J+1: Er2 14, Er3 12, Er4 10, Er5 10, Er6s 4, Er6h 12, Er7 8, Er8 10, Er9 12
    norad9 = TemperatureDependentInternalConversionRef(
        source=Er9, target=Er8, reference_rate=k_nr9, T=T_eff, T_A=T0,
        phonon_wavenumber=phonon_wavenumber, g_source=12, g_target=10,
        detailed_balance=detailed_balance,
    )
    norad8 = TemperatureDependentInternalConversionRef(
        source=Er8, target=Er7, reference_rate=k_nr8, T=T_eff, T_A=T0,
        phonon_wavenumber=phonon_wavenumber, g_source=10, g_target=8,
        detailed_balance=detailed_balance,
    )
    norad7 = TemperatureDependentInternalConversionRef(
        source=Er7, target=Er6h, reference_rate=k_nr7, T=T_eff, T_A=T0,
        phonon_wavenumber=phonon_wavenumber, g_source=8, g_target=12,
        detailed_balance=detailed_balance,
    )
    norad6h = TemperatureDependentInternalConversionRef(
        source=Er6h, target=Er5, reference_rate=k_nr6, T=T_eff, T_A=T0,
        phonon_wavenumber=phonon_wavenumber, g_source=12, g_target=10,
        detailed_balance=detailed_balance,
    )
    norad6s = TemperatureDependentInternalConversionRef(
        source=Er6s, target=Er5, reference_rate=k_nr6, T=T_eff, T_A=T0,
        phonon_wavenumber=phonon_wavenumber, g_source=4, g_target=10,
        detailed_balance=detailed_balance,
    )
    norad5 = TemperatureDependentInternalConversionRef(
        source=Er5, target=Er4, reference_rate=k_nr5, T=T_eff, T_A=T0,
        phonon_wavenumber=phonon_wavenumber, g_source=10, g_target=10,
        detailed_balance=detailed_balance,
    )
    norad4 = TemperatureDependentInternalConversionRef(
        source=Er4, target=Er3, reference_rate=k_nr4, T=T_eff, T_A=T0,
        phonon_wavenumber=phonon_wavenumber, g_source=10, g_target=12,
        detailed_balance=detailed_balance,
    )
    norad3 = TemperatureDependentInternalConversionRef(
        source=Er3, target=Er2, reference_rate=k_nr3, T=T_eff, T_A=T0,
        phonon_wavenumber=phonon_wavenumber, g_source=12, g_target=14,
        detailed_balance=detailed_balance,
    )

    # ---- radiative decays (Anderson's branching ratios) ----
    rad81 = Fluorescence(excited=Er8, ground=Er1, rate=0.40 * k_r8)  # blue ~410 nm
    rad82 = Fluorescence(excited=Er8, ground=Er2, rate=0.42 * k_r8)  # 2H9/2->4I13/2 ~556 nm
    rad83 = Fluorescence(excited=Er8, ground=Er3, rate=0.14 * k_r8)
    rad85 = Fluorescence(excited=Er8, ground=Er5, rate=0.04 * k_r8)

    rad6h1 = Fluorescence(excited=Er6h, ground=Er1, rate=0.70 * k_r6h)  # 2H11/2 ~525 nm
    rad6h2 = Fluorescence(excited=Er6h, ground=Er2, rate=0.25 * k_r6h)
    rad6h3 = Fluorescence(excited=Er6h, ground=Er3, rate=0.05 * k_r6h)
    rad6s1 = Fluorescence(excited=Er6s, ground=Er1, rate=0.70 * k_r6s)  # 4S3/2 ~545 nm
    rad6s2 = Fluorescence(excited=Er6s, ground=Er2, rate=0.25 * k_r6s)
    rad6s3 = Fluorescence(excited=Er6s, ground=Er3, rate=0.05 * k_r6s)

    rad51 = Fluorescence(excited=Er5, ground=Er1, rate=0.90 * k_r5)  # red ~660 nm
    rad52 = Fluorescence(excited=Er5, ground=Er2, rate=0.05 * k_r5)
    rad53 = Fluorescence(excited=Er5, ground=Er3, rate=0.05 * k_r5)

    rad31 = Fluorescence(excited=Er3, ground=Er1, rate=0.81 * k_r3)  # 1 um
    rad32 = Fluorescence(excited=Er3, ground=Er2, rate=0.19 * k_r3)
    rad2 = Fluorescence(excited=Er2, ground=Er1, rate=k_r2)  # 1.5 um

    # ---- Er-Er cross relaxation and ETU (constant, as in Anderson) ----
    cross6s = EnergyTransferUpconversion(
        sensitizer_high=Er6s, sensitizer_low=Er3, activator_low=Er1, activator_high=Er2,
        rate=k_cr6 * site_density,
    )
    cross6h = EnergyTransferUpconversion(
        sensitizer_high=Er6h, sensitizer_low=Er3, activator_low=Er1, activator_high=Er2,
        rate=k_cr6 * site_density,
    )
    cross4 = EnergyTransferUpconversion(
        sensitizer_high=Er4, sensitizer_low=Er2, activator_low=Er1, activator_high=Er2,
        rate=k_cr4 * site_density,
    )
    uc2 = EnergyTransferUpconversion(
        sensitizer_high=Er2, sensitizer_low=Er1, activator_low=Er2, activator_high=Er4,
        rate=k_uc2 * site_density,
    )
