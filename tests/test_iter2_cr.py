"""Iteración 2: con cr_thermal = 0 el modelo es el de la iteración 1 a cualquier T."""

import json
import pathlib
import sys

import numpy as np
import pytest
from pint import get_application_registry

from upconversion_models.anderson_thermal import AndersonThermal as M

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "noche"))
from sweep import make_sim, observables  # noqa: E402

from upconversion_models.steady import steady_state  # noqa: E402

u = get_application_registry()
ITER1 = json.loads((ROOT / "noche" / "data" / "iter1_sweep.json").read_text())
SIM = make_sim(M)


@pytest.mark.parametrize("T", [20.0, 160.0, 300.0, 500.0, 700.0])
@pytest.mark.parametrize("P", [0.1, 290.0])
def test_cr_thermal_off_is_iteration1(T, P):
    r = observables(steady_state(SIM, {M.T: T * u.K, M.energy_flux: P * u.W / u.cm**2, M.cr_thermal: 0, M.et_thermal: 0}))
    ref = ITER1["obs"][f"{T}|{P}"]
    for k, v in ref.items():
        assert np.isclose(r[k], v, rtol=1e-6), (k, r[k], v)


def test_cr6_mismatch_from_levels():
    assert np.isclose(M.cross6s._mismatch_wavenumber.default.to(1 / u.cm).magnitude, 1600)
    assert np.isclose(M.cross6h._mismatch_wavenumber.default.to(1 / u.cm).magnitude, 2250)
    assert M.cross6s._emits.default == 1 and M.cross6s._one_phonon.default == 0
