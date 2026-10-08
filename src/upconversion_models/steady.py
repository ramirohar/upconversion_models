"""Precise steady states: time integration followed by a Newton polish.

`jablonski_patch.SteadyState` stops the integration when the derivatives are
small (SmallDerivatives(atol=1e-12, rtol=1e-6)), which leaves a relative error of
order 1e-7 in the slow populations. Regression tests at rtol = 1e-6 need more
margin, so here the integrated state is refined by Newton iterations on
f(y) = 0 constrained by the conservation laws of the system (left null space of
the Jacobian, i.e. total Yb and total Er populations).

Fixed numerics (Fase 0 of PLAN_NOCHE.md):
- integration: scipy LSODA, rtol=1e-10, atol=1e-16, up to t_end = 10 s or until
  max|dy/dt| / max(y) < 1e-7 s^-1 (the Newton polish does the rest);
- Newton: central-difference Jacobian, stop at max|f_i| / (|y_i| · k_i) < 1e-10
  (k_i: total loss rate of state i) or 30 iterations.
"""

from __future__ import annotations

import re
from typing import Mapping

import numpy as np
import pint
from scipy.integrate import solve_ivp

u = pint.get_application_registry()

RTOL_INT = 1e-10
ATOL_INT = 1e-16
T_END = 10.0
NEWTON_TOL = 1e-10
NEWTON_MAXITER = 30


def _jacobian(f, y, floor=1e-12):
    n = y.size
    J = np.empty((n, n))
    for i in range(n):
        h = 1e-6 * max(abs(y[i]), floor)
        yp, ym = y.copy(), y.copy()
        yp[i] += h
        ym[i] -= h
        J[:, i] = (f(yp) - f(ym)) / (2 * h)
    return J


def _conservation(names, f, y):
    """Conserved totals: one row per ion species, grouping states by the
    alphabetic prefix of their name (Yb1, Yb2 -> Yb; Er1 ... Er9, Er6s -> Er).

    An SVD of the Jacobian misreads the rank when rates span ~10 decades (slow
    modes look null), so the groups are taken from the names and checked: each
    total must have zero time derivative."""
    prefixes = [re.match(r"[A-Za-z]+?(?=\d)", n.split(".")[-1]) for n in names]
    groups = sorted({m.group(0) for m in prefixes if m})
    W = np.array(
        [[1.0 if m and m.group(0) == g else 0.0 for m in prefixes] for g in groups]
    )
    d = W @ f(y)
    assert np.all(np.abs(d) <= 1e-6 * np.max(np.abs(y))), d  # roundoff of fast rates
    return W


def steady_state(sim, values: Mapping = {}, *, return_info: bool = False):
    """Steady state of `sim` (a jablonski_patch.Simulator) with `values` applied.

    Returns a dict {output name: magnitude} with every transform output of the
    simulator (state populations and, if present, `line_*` emission rates).
    """
    prob = sim.with_values(values).create_problem()
    p = prob.p
    dy = np.empty_like(prob.y)

    def f(y):
        return np.array(prob.rhs(0.0, y, p, dy), dtype=float).copy()

    y0 = prob.y.copy()

    def small(t, y):
        d = f(y)
        return np.max(np.abs(d)) / max(np.max(np.abs(y)), 1e-300) - 1e-7

    small.terminal = True
    sol = solve_ivp(
        lambda t, y: f(y),
        (0.0, T_END),
        y0,
        method="LSODA",
        rtol=RTOL_INT,
        atol=ATOL_INT,
        events=small,
    )
    y = sol.y[:, -1].copy()

    # Newton polish with the conservation laws of the initial condition
    J = _jacobian(f, y)
    W = _conservation([str(v) for v in sim.compiled.variables], f, y)
    c = W @ y0
    it = 0
    for it in range(NEWTON_MAXITER):
        F = f(y)
        J = _jacobian(f, y)
        loss = np.maximum(np.abs(np.diag(J)), 1e-300)
        err = np.max(np.abs(F) / np.maximum(np.abs(y) * loss, 1e-300))
        if err < NEWTON_TOL:
            break
        # row scaling by the loss rate and column scaling by the population keep
        # the system well conditioned when rates span many decades
        s = np.maximum(np.abs(y), 1e-30)
        A = np.vstack([J / loss[:, None], W]) * s[None, :]
        b = -np.concatenate([F / loss, W @ y - c])
        step, *_ = np.linalg.lstsq(A, b, rcond=None)
        y = y + s * step
    F = f(y)
    loss = np.maximum(np.abs(np.diag(_jacobian(f, y))), 1e-300)
    err = float(np.max(np.abs(F) / np.maximum(np.abs(y) * loss, 1e-300)))

    out = np.empty((len(prob.scale), 1))
    out = prob.transform(np.array([sol.t[-1]]), y[:, None], p, out)
    names = [str(k) for k in sim.transform.output.keys()]
    res = {}
    for name, s, val in zip(names, prob.scale, np.asarray(out).reshape(len(names), -1)[:, 0]):
        mag = s.magnitude if isinstance(s, pint.Quantity) else s
        res[name] = float(val * mag)
    if return_info:
        return res, {"t_int": float(sol.t[-1]), "newton_iter": it, "residual": err}
    return res


def steady_sweep(sim, variable, values, extra: Mapping = {}):
    """Steady states over `values` of `variable`; returns {output: np.array}."""
    rows = [steady_state(sim, {**extra, variable: v}) for v in values]
    return {k: np.array([r[k] for r in rows]) for k in rows[0]}
