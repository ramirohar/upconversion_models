from jablonski import (
    SingletState,
    SpectroscopicSystem,
    Parameter,
    initial,
    assign,
)
from jablonski.transitions import (
    Absorption,
    Fluorescence,
    MassAction,
    InternalConversion,
)
from pint import get_application_registry

u = get_application_registry()

h = 6.62607015e-34 * u.J * u.s  # Planks constant
c = 299_792_458 * u.m / u.s  # Speed of light
hc = h * c


# preguntas: que es el isolated precentaje y de donde sale?


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


class YbIon(SpectroscopicSystem):
    ## Parameter values taken from 2013-Anderson's model

    Yb_cross_section: Parameter = assign(default=1e-20 * u.cm**2)  # no value reported
    Yb_rad: Parameter = assign(default=393 / u.s)
    Yb_nr: Parameter = assign(default=220 / u.s)

    isolated_percentage: Parameter = assign(default=0.05)
    Yb_concetration: Parameter = assign(default=0.18)

    Yb1: SingletState = initial(
        energy=(0 * hc / u.cm),
        spin_multiplicity="singlet",
        default=Yb_concetration * (1 - isolated_percentage),  # type: ignore
    )

    Yb2: SingletState = initial(
        energy=10200 * hc / u.cm,
        spin_multiplicity="singlet",
        default=0,
    )

    abs = Absorption(
        ground=Yb1,
        excited=Yb2,
        rate=Yb_cross_section,
        # pump=# photon_flux,
    )

    rad = Fluorescence(ground=Yb1, excited=Yb2, rate=Yb_rad)
    norad = InternalConversion(low=Yb1, high=Yb2, rate=Yb_nr)


class ErIon(SpectroscopicSystem):
    ## Parameter values taken from 2013-Anderson's model

    er_concetration: Parameter = assign(default=0.02)
    isolated_percentage: Parameter = assign(default=0.05)

    Er1: SingletState = initial(
        energy=0 * hc / u.cm,
        spin_multiplicity="singlet",
        default=er_concetration * (1 - isolated_percentage),  # type: ignore
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
    k_nr5: Parameter = assign(default=0 / u.s)
    k_nr4: Parameter = assign(default=22_120 / u.s)
    k_nr3: Parameter = assign(default=61 / u.s)

    # Radiative rates
    k_r8: Parameter = assign(default=2_330 / u.s)
    k_r6: Parameter = assign(default=1_510 / u.s)
    k_r5: Parameter = assign(default=2_039 / u.s)
    k_r3: Parameter = assign(default=73 / u.s)
    k_r2: Parameter = assign(default=110 / u.s)

    # Cross relaxation
    k_cr6: Parameter = assign(default=2.79e-17 / u.s * u.cm**3)
    k_cr4: Parameter = assign(default=8.04e-19 / u.s * u.cm**3)

    #  6s 6h thermalization
    fs: Parameter = assign(default=0.5)  # To be complete
    fh: Parameter = assign(default=1 - fs)

    # Level 6 partition
    a_r: Parameter = assign(default=1)
    b_r = (1 - fs * a_r) / (1 - fs)

    a_nr: Parameter = assign(default=1)
    b_nr = (1 - fs * a_nr) / (1 - fs)

    a_et9: Parameter = assign(default=1)
    b_et9 = (1 - fs * a_et9) / (1 - fs)

    a_cr: Parameter = assign(default=1)
    b_cr = (1 - fs * a_cr) / (1 - fs)

    # Parameters introduced (not reported in Anderson's)
    site_density: Parameter = assign(default=1.38e22 / u.cm**3)

    # Additional parameters
    k_uc2: Parameter = assign(default=2.31e-17 / u.s * u.cm**3)

    # Non radiative decays

    norad9 = InternalConversion(high=Er9, low=Er8, rate=k_nr9)
    norad8 = InternalConversion(high=Er8, low=Er7, rate=k_nr8)
    norad7h = InternalConversion(high=Er7, low=Er6h, rate=k_nr7)
    # norad7s = InternalConversion(high=Er7, low=Er6s, rate=a_nr*k_nr7)
    norad6h = InternalConversion(high=Er6h, low=Er5, rate=b_nr * k_nr6)
    norad6s = InternalConversion(high=Er6s, low=Er5, rate=a_nr * k_nr6)
    norad5 = InternalConversion(high=Er5, low=Er4, rate=k_nr5)
    norad4 = InternalConversion(high=Er4, low=Er3, rate=k_nr4)
    norad3 = InternalConversion(high=Er3, low=Er2, rate=k_nr3)

    # Radiative decays

    rad81 = Fluorescence(excited=Er8, ground=Er1, rate=0.40 * k_r8)  # Blue
    rad82 = Fluorescence(excited=Er8, ground=Er2, rate=0.42 * k_r8)  # Blue
    rad83 = Fluorescence(excited=Er8, ground=Er3, rate=0.14 * k_r8)  # Blue
    rad85 = Fluorescence(excited=Er8, ground=Er5, rate=0.04 * k_r8)  # Blue

    rad6h1 = Fluorescence(excited=Er6h, ground=Er1, rate=0.70 * b_r * k_r6)  # Green
    rad6h2 = Fluorescence(excited=Er6h, ground=Er2, rate=0.25 * b_r * k_r6)  # Green
    rad6h3 = Fluorescence(excited=Er6h, ground=Er3, rate=0.05 * b_r * k_r6)  # Green

    rad6s1 = Fluorescence(
        excited=Er6s, ground=Er1, rate=0.70 * a_r * k_r6
    )  # Yellow-Green
    rad6s2 = Fluorescence(
        excited=Er6s, ground=Er2, rate=0.25 * a_r * k_r6
    )  # Yellow-Green
    rad6s3 = Fluorescence(
        excited=Er6s, ground=Er3, rate=0.05 * a_r * k_r6
    )  # Yellow-Green

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
        rate=k_cr6 * b_cr * site_density,
    )
    cross6s = EnergyTransferUpconversion(
        sensitizer_high=Er6s,
        sensitizer_low=Er3,
        activator_low=Er1,
        activator_high=Er2,
        rate=k_cr6 * a_cr * site_density,
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


class UCNP(SpectroscopicSystem):
    Yb: SpectroscopicSystem = YbIon()
    Er: SpectroscopicSystem = ErIon()

    # Level 6 poblation parameters
    fs: Parameter = assign(default=0.5)  # To be complete
    fh: Parameter = assign(default=1 - fs)

    a_et9: Parameter = assign(default=1)
    b_et9 = (1 - fs * a_et9) / (1 - fs)

    # Parameters introduced (not reported in Anderson's)
    site_density: Parameter = assign(default=1.38e22 / u.cm**3)

    # Yb -> Er etu rates
    k_13: Parameter = assign(default=5.9 * 2e-16 / u.s * u.cm**3)
    k_25: Parameter = assign(default=0 / u.s * u.cm**3)
    k_37: Parameter = assign(default=1.54e-15 / u.s * u.cm**3)
    k_58: Parameter = assign(default=1.76e-15 / u.s * u.cm**3)
    k_69: Parameter = assign(default=6.07e-15 / u.s * u.cm**3)

    # Er -> Yb etu rates
    k_31: Parameter = assign(default=2e-16 / u.s * u.cm**3)
    k_52: Parameter = assign(default=0 / u.s * u.cm**3)
    k_73: Parameter = assign(default=2.04e-16 / u.s * u.cm**3)
    k_95: Parameter = assign(default=2.84e-16 / u.s * u.cm**3)

    # Yb -> Er transitions

    etu13 = EnergyTransferUpconversion(
        sensitizer_high=Yb.Yb2,
        sensitizer_low=Yb.Yb1,
        activator_low=Er.Er1,
        activator_high=Er.Er3,
        rate=k_13 * site_density,
    )

    etu25 = EnergyTransferUpconversion(
        sensitizer_high=Yb.Yb2,
        sensitizer_low=Yb.Yb1,
        activator_low=Er.Er2,
        activator_high=Er.Er5,
        rate=k_25 * site_density,
    )

    etu37 = EnergyTransferUpconversion(
        sensitizer_high=Yb.Yb2,
        sensitizer_low=Yb.Yb1,
        activator_low=Er.Er3,
        activator_high=Er.Er7,
        rate=k_37 * site_density,
    )

    etu58 = EnergyTransferUpconversion(
        sensitizer_high=Yb.Yb2,
        sensitizer_low=Yb.Yb1,
        activator_low=Er.Er5,
        activator_high=Er.Er8,
        rate=k_58 * site_density,
    )

    etu6s9 = EnergyTransferUpconversion(
        sensitizer_high=Yb.Yb2,
        sensitizer_low=Yb.Yb1,
        activator_low=Er.Er6s,
        activator_high=Er.Er9,
        rate=k_69 * a_et9 * site_density,
    )

    etu6h9 = EnergyTransferUpconversion(
        sensitizer_high=Yb.Yb2,
        sensitizer_low=Yb.Yb1,
        activator_low=Er.Er6h,
        activator_high=Er.Er9,
        rate=k_69 * b_et9 * site_density,
    )

    # Er -> Yb transitions

    etu31 = EnergyTransferUpconversion(
        sensitizer_high=Er.Er3,
        sensitizer_low=Er.Er1,
        activator_low=Yb.Yb1,
        activator_high=Yb.Yb2,
        rate=k_31 * site_density,
    )

    etu52 = EnergyTransferUpconversion(
        sensitizer_high=Er.Er5,
        sensitizer_low=Er.Er2,
        activator_low=Yb.Yb1,
        activator_high=Yb.Yb2,
        rate=k_52 * site_density,
    )

    etu73 = EnergyTransferUpconversion(
        sensitizer_high=Er.Er7,
        sensitizer_low=Er.Er3,
        activator_low=Yb.Yb1,
        activator_high=Yb.Yb2,
        rate=k_73 * site_density,
    )

    etu95 = EnergyTransferUpconversion(
        sensitizer_high=Er.Er9,
        sensitizer_low=Er.Er5,
        activator_low=Yb.Yb1,
        activator_high=Yb.Yb2,
        rate=k_95 * site_density,
    )
