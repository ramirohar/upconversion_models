from jablonski import (
    Parameter,
    SingletState,
    SpectroscopicSystem,
    assign,
    initial,
)
from jablonski.transitions import (
    Absorption,
    Fluorescence,
    InternalConversion,
)
from pint import get_application_registry

from ..utils import hc

u = get_application_registry()

class YbIon(SpectroscopicSystem):
    ## Parameter values taken from 2013-Anderson's model

    cross_section: Parameter = assign(default=1e-20 * u.cm**2)  # no value reported
    k_r: Parameter = assign(default=393 / u.s)
    k_nr: Parameter = assign(default=220 / u.s)

    isolated_percentage: Parameter = assign(default=0.05)
    concetration: Parameter = assign(default=0.18)

    Yb1: SingletState = initial(
        energy=(0 * hc / u.cm),
        spin_multiplicity="singlet",
        default=concetration * (1 - isolated_percentage),  # type: ignore
    )

    Yb2: SingletState = initial(
        energy=10200 * hc / u.cm,
        spin_multiplicity="singlet",
        default=0,
    )

    abs = Absorption(
        ground=Yb1,
        excited=Yb2,
        rate=cross_section,
        # pump = photon_flux,
    )

    rad = Fluorescence(ground=Yb1, excited=Yb2, rate=k_r)
    norad = InternalConversion(low=Yb1, high=Yb2, rate=k_nr)

