"""Iteración 4: con nr_thermal = 0 el modelo es el de la iteración 3 a cualquier T."""

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
ITER3 = json.loads((ROOT / "noche" / "data" / "iter3_sweep.json").read_text())
SIM = make_sim(M)


@pytest.mark.parametrize("T", [20.0, 160.0, 300.0, 500.0, 700.0])
@pytest.mark.parametrize("P", [0.1, 290.0])
def test_nr_thermal_off_is_iteration3(T, P):
    r = observables(steady_state(SIM, {M.T: T * u.K, M.energy_flux: P * u.W / u.cm**2, M.nr_thermal: 0}))
    ref = ITER3["obs"][f"{T}|{P}"]
    for k, v in ref.items():
        assert np.isclose(r[k], v, rtol=1e-6), (k, r[k], v)


@pytest.mark.parametrize("name,gap,emits", [("etu13", 40, 0), ("etu37", 90, 0), ("etu31", 40, 1), ("etu73", 90, 1)])
def test_li2014_mismatch(name, gap, emits):
    block = getattr(M, name)
    assert np.isclose(block._mismatch_wavenumber.default.to(1 / u.cm).magnitude, gap)
    assert block._emits.default == emits and block._one_phonon.default == 1
