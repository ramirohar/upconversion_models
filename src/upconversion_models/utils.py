import pint
from pint import get_application_registry

u = get_application_registry()

h = 6.62607015e-34 * u.J * u.s  # Planks constant
c = 299_792_458 * u.m / u.s  # Speed of light
k = 8.617333262e-5 * u.eV / u.K # Boltzman's constant

hc = h*c

def energy_from_wavelength(wavelength: pint.Quantity):
    return h * c / wavelength
