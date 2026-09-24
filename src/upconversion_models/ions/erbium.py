from jablonski import (
    Parameter,
    SingletState,
    SpectroscopicSystem,
    assign,
    initial,
)
from jablonski.transitions import (
    Fluorescence,
    InternalConversion,
)
from pint import get_application_registry

from ..transitions import (
    EnergyTransferUpconversion,
)
from ..utils import hc

u = get_application_registry()


class ErIon(SpectroscopicSystem):
    ## Parameter values taken from 2013-Anderson's model

    concetration: Parameter = assign(default=0.02)
    isolated_percentage: Parameter = assign(default=0.05)

    Er1: SingletState = initial(
        energy=0 * hc / u.cm,
        spin_multiplicity="singlet",
        default=concetration * (1 - isolated_percentage),  # type: ignore
    )

    Er2: SingletState = initial(
        energy=6500 * hc / u.cm,
        spin_multiplicity="singlet",
        default=0,
    )

    Er3: SingletState = initial(
        energy=10200 * hc / u.cm,
        spin_multiplicity="singlet",
        default=0,
    )

    Er4: SingletState = initial(
        energy=12500 * hc / u.cm,
        spin_multiplicity="singlet",
        default=0,
    )

    Er5: SingletState = initial(
        energy=15000 * hc / u.cm,
        spin_multiplicity="singlet",
        default=0,
    )

    Er6s: SingletState = initial(
        energy=18300 * hc / u.cm,
        spin_multiplicity="singlet",
        default=0,
    )

    Er6h: SingletState = initial(
        energy=19200 * hc / u.cm,
        spin_multiplicity="singlet",
        default=0,
    )

    Er7: SingletState = initial(
        energy=20500 * hc / u.cm,
        spin_multiplicity="singlet",
        default=0,
    )

    Er8: SingletState = initial(
        energy=24500 * hc / u.cm,
        spin_multiplicity="singlet",
        default=0,
    )

    Er9: SingletState = initial(
        energy=26100 * hc / u.cm,
        spin_multiplicity="singlet",
        default=0,
    )

    # Non radiative rates
    k_nr9: Parameter = assign(default=1.76e6 / u.s)
    k_nr8: Parameter = assign(default=43_450 / u.s)
    k_nr7: Parameter = assign(default=1e6 / u.s)
    k_nr6: Parameter = assign(default=26 / u.s)
    k_nr6s: Parameter = assign(default=26 / u.s)
    k_nr6h: Parameter = assign(default=26 / u.s)
    k_nr5: Parameter = assign(default=0 / u.s)
    k_nr4: Parameter = assign(default=22_120 / u.s)
    k_nr3: Parameter = assign(default=61 / u.s)

    # Radiative rates
    k_r8: Parameter = assign(default=2_330 / u.s)
    k_r6: Parameter = assign(default=1_510 / u.s)
    k_r6s: Parameter = assign(default=1_510 / u.s)
    k_r6h: Parameter = assign(default=1_510 / u.s)
    k_r5: Parameter = assign(default=2_039 / u.s)
    k_r3: Parameter = assign(default=73 / u.s)
    k_r2: Parameter = assign(default=110 / u.s)

    # Cross relaxation
    k_cr6: Parameter = assign(default=2.79e-17 / u.s * u.cm**3)
    k_cr6s: Parameter = assign(default=2.79e-17 / u.s * u.cm**3)
    k_cr6h: Parameter = assign(default=2.79e-17 / u.s * u.cm**3)
    k_cr4: Parameter = assign(default=8.04e-19 / u.s * u.cm**3)

    # Additional parameters
    k_uc2: Parameter = assign(default=2.31e-17 / u.s * u.cm**3)

    # Level 6 poblation parameters
    fs: Parameter = assign(default=0.9)  # To be complete
    fh: Parameter = assign(default=1 - fs)

    # Manifold rates
    k_therm: Parameter = assign(default=1e13 / u.s)
    k_manifold_up: Parameter = assign(default=fh * k_therm)
    k_manifold_down: Parameter = assign(default=fs * k_therm)

    # Manifold rates
    k_therm: Parameter = assign(default=1e13 / u.s)
    k_manifold_up: Parameter = assign(default=fh * k_therm)
    k_manifold_down: Parameter = assign(default=fs * k_therm)

    # Parameters introduced (not reported in Anderson's)
    site_density: Parameter = assign(default=1.38e22 / u.cm**3)

    # Additional parameters
    k_uc2: Parameter = assign(default=2.31e-17 / u.s * u.cm**3)

    # Non radiative decays

    norad9 = InternalConversion(high=Er9, low=Er8, rate=k_nr9)
    norad8 = InternalConversion(high=Er8, low=Er7, rate=k_nr8)
    norad7 = InternalConversion(high=Er7, low=Er6h, rate=k_nr7)
    norad6h = InternalConversion(high=Er6h, low=Er5, rate=k_nr6s)
    norad6s = InternalConversion(high=Er6s, low=Er5, rate=k_nr6h)
    norad5 = InternalConversion(high=Er5, low=Er4, rate=k_nr5)
    norad4 = InternalConversion(high=Er4, low=Er3, rate=k_nr4)
    norad3 = InternalConversion(high=Er3, low=Er2, rate=k_nr3)

    # Radiative decays

    rad81 = Fluorescence(excited=Er8, ground=Er1, rate=0.40 * k_r8)  # Blue
    rad82 = Fluorescence(excited=Er8, ground=Er2, rate=0.42 * k_r8)  # Blue
    rad83 = Fluorescence(excited=Er8, ground=Er3, rate=0.14 * k_r8)  # Blue
    rad85 = Fluorescence(excited=Er8, ground=Er5, rate=0.04 * k_r8)  # Blue

    rad6h1 = Fluorescence(excited=Er6h, ground=Er1, rate=0.70 * k_r6h)  # Green
    rad6h2 = Fluorescence(excited=Er6h, ground=Er2, rate=0.25 * k_r6h)  # Green
    rad6h3 = Fluorescence(excited=Er6h, ground=Er3, rate=0.05 * k_r6h)  # Green

    rad6s1 = Fluorescence(excited=Er6s, ground=Er1, rate=0.70 * k_r6s)  # Yellow-Green
    rad6s2 = Fluorescence(excited=Er6s, ground=Er2, rate=0.25 * k_r6s)  # Yellow-Green
    rad6s3 = Fluorescence(excited=Er6s, ground=Er3, rate=0.05 * k_r6s)  # Yellow-Green

    rad51 = Fluorescence(excited=Er5, ground=Er1, rate=0.90 * k_r5)  # Red
    rad52 = Fluorescence(excited=Er5, ground=Er2, rate=0.05 * k_r5)  # Red
    rad53 = Fluorescence(excited=Er5, ground=Er3, rate=0.05 * k_r5)  # Red

    rad31 = Fluorescence(excited=Er3, ground=Er1, rate=0.81 * k_r3)  # 1um
    rad32 = Fluorescence(excited=Er3, ground=Er2, rate=0.19 * k_r3)  # 1um

    rad2 = Fluorescence(excited=Er2, ground=Er1, rate=k_r2)  # 1.5um

    # Er-Er Cross-Relaxation

    cross6h = EnergyTransferUpconversion(
        sensitizer_high=Er6h,
        sensitizer_low=Er3,
        activator_low=Er1,
        activator_high=Er2,
        rate=k_cr6h * site_density,
    )
    cross6s = EnergyTransferUpconversion(
        sensitizer_high=Er6s,
        sensitizer_low=Er3,
        activator_low=Er1,
        activator_high=Er2,
        rate=k_cr6s * site_density,
    )
    cross4 = EnergyTransferUpconversion(
        sensitizer_high=Er4,
        sensitizer_low=Er2,
        activator_low=Er1,
        activator_high=Er2,
        rate=k_cr4 * site_density,
    )

    # Er -> Er ETU
    uc2 = EnergyTransferUpconversion(
        sensitizer_high=Er2,
        sensitizer_low=Er1,
        activator_low=Er2,
        activator_high=Er4,
        rate=k_uc2 * site_density,
    )
