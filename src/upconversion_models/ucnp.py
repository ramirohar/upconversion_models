from jablonski import (
    Parameter,
    SpectroscopicSystem,
    assign,
)

from pint import get_application_registry

from .transitions import (
    EnergyTransferUpconversion,
)
from .ions.erbium import ErIon
from .ions.yterbium import YbIon

from .utils import c, h

u = get_application_registry()

hc = h * c


class UCNP(SpectroscopicSystem):
    Yb: SpectroscopicSystem = YbIon()
    Er: SpectroscopicSystem = ErIon()

    # Level 6 poblation parameters
    fs: Parameter = assign(default=0.5)  # To be complete
    fh: Parameter = assign(default=1 - fs)

    # Parameters introduced (not reported in Anderson's)
    site_density: Parameter = assign(default=1.38e22 / u.cm**3)

    # Yb -> Er etu rates
    k_13: Parameter = assign(default=5.9 * 2e-16 / u.s * u.cm**3)
    k_25: Parameter = assign(default=0 / u.s * u.cm**3)
    k_37: Parameter = assign(default=1.54e-15 / u.s * u.cm**3)
    k_58: Parameter = assign(default=1.76e-15 / u.s * u.cm**3)
    k_69: Parameter = assign(default=6.07e-15 / u.s * u.cm**3)
    k_6s9: Parameter = assign(default=6.07e-15 / u.s * u.cm**3)
    k_6h9: Parameter = assign(default=6.07e-15 / u.s * u.cm**3)

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
        rate=k_6s9 * site_density,
    )

    etu6h9 = EnergyTransferUpconversion(
        sensitizer_high=Yb.Yb2,
        sensitizer_low=Yb.Yb1,
        activator_low=Er.Er6h,
        activator_high=Er.Er9,
        rate=k_6h9 * site_density,
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
