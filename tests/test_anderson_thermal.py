"""Iteración 1: AndersonThermal se reduce a AndersonModel.

- termalización rápida y T = T0: se recupera la línea base;
- dependencia térmica apagada (thermal = 0) y termalización rápida: se recupera
  Anderson a cualquier T.
"""

import json
import pathlib

import numpy as np
import pytest
from pint import get_application_registry

from upconversion_models.anderson_thermal import AndersonThermal as M
from upconversion_models.jablonski_patch import EmissionTransform, Simulator
from upconversion_models.steady import steady_state

u = get_application_registry()
ROOT = pathlib.Path(__file__).resolve().parents[1]
BASELINE = json.loads((ROOT / "noche" / "data" / "baseline_anderson.json").read_text())
SIM = Simulator(M).with_transform(EmissionTransform("emission"), append=True)
FAST = 1e10 / u.s  # 4e10-1.2e11 s^-1 per direction; 1e11 makes Newton ill-conditioned


def lumped(r):
    """Map AndersonThermal outputs onto AndersonModel names."""
    out = dict(r)
    out["Er6"] = r["Er6s"] + r["Er6h"]
    for j in (1, 2, 3):
        out[f"line_rad6{j}"] = r[f"line_rad6s{j}"] + r[f"line_rad6h{j}"]
    return out


def compare(power, values):
    r = lumped(steady_state(SIM, {M.energy_flux: float(power) * u.W / u.cm**2, **values}))
    ref = BASELINE[power]["values"]
    for key, value in ref.items():
        assert np.isclose(r[key], value, rtol=1e-6, atol=0), (key, r[key], value)


@pytest.mark.parametrize("power", ["0.1", "10", "1000"])
def test_fast_thermalization_at_T0(power):
    compare(power, {M.k_nr0: FAST, M.T: 300 * u.K, M.detailed_balance: 0})


@pytest.mark.parametrize("T", [20, 120, 500, 700])
@pytest.mark.parametrize("power", ["1", "100"])
def test_thermal_off_any_T(power, T):
    compare(power, {M.k_nr0: FAST, M.T: T * u.K, M.thermal: 0, M.detailed_balance: 0})


def test_boltzmann_partition_fast():
    """Fast thermalization: Er6h/Er6s = (g_H/g_S) exp(-ΔE/kT) at any T."""
    for T in (300, 600):  # at 100 K the 4F7/2 feed of 2H11/2 already breaks equilibrium
        r = steady_state(
            SIM, {M.energy_flux: 10 * u.W / u.cm**2, M.k_nr0: FAST, M.T: T * u.K}
        )
        expected = 3 * np.exp(-650 * 1.4387769 / T)
        assert np.isclose(r["Er6h"] / r["Er6s"], expected, rtol=1e-4)
