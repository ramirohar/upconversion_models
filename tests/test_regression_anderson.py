"""Regresión de AndersonModel contra la línea base de la Fase 0 (rtol = 1e-6)."""

import json
import pathlib
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "noche"))

from fase0_baseline import OUT, compute  # noqa: E402

BASELINE = json.loads(OUT.read_text())


@pytest.mark.parametrize("power", list(BASELINE))
def test_steady_state_matches_baseline(power):
    new = compute([float(power)])[str(float(power))]["values"]
    ref = BASELINE[power]["values"]
    for key, value in ref.items():
        assert np.isclose(new[key], value, rtol=1e-6, atol=0), key
