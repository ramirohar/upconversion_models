"""Dinámica tras un pulso corto de 980 nm (O11: tiempos de subida y decaimiento, Yu 2014 Fig. 7).

El pulso de Yu (10 ns) es mucho más corto que toda la dinámica (µs–ms), así que se modela como
una excitación instantánea de una fracción `x0` de los Yb. Se integra con LSODA y, para cada banda,
se obtiene τ_D de la pendiente logarítmica tardía y τ_R del tiempo del máximo,
t_peak = τ_R·τ_D/(τ_D − τ_R)·ln(τ_D/τ_R) (forma de Vial, Yu 2014 eq. 3).
"""

from __future__ import annotations

import numpy as np
from pint import get_application_registry
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

from upconversion_models.jablonski_patch import EmissionTransform, Simulator

u = get_application_registry()

BANDS = {"green": ("line_rad6h1", "line_rad6s1"), "red": ("line_rad51",), "blue": ("line_rad81",)}


def pulse_curves(model, T, x0=0.01, extra=None, t_end=5e-3, n=4000):
    sim = Simulator(model).with_transform(EmissionTransform("emission"), append=True)
    prob = sim.with_values({model.T: T * u.K, model.energy_flux: 0 * u.W / u.cm**2, **(extra or {})}).create_problem()
    names = [str(v) for v in sim.compiled.variables]
    y0 = prob.y.copy()
    iyb1, iyb2 = names.index("Yb1"), names.index("Yb2")
    exc = x0 * y0[iyb1]
    y0[iyb1] -= exc
    y0[iyb2] += exc
    dy = np.empty_like(y0)
    t = np.r_[0, np.geomspace(1e-8, t_end, n)]
    sol = solve_ivp(lambda tt, y: np.array(prob.rhs(tt, y, prob.p, dy)).copy(), (0, t_end), y0,
                    method="LSODA", t_eval=t, rtol=1e-8, atol=1e-18)
    out = np.empty((len(prob.scale), t.size))
    out = prob.transform(sol.t, sol.y, prob.p, out)
    keys = [str(k) for k in sim.transform.output.keys()]
    lines = {k: np.asarray(out)[i] for i, k in enumerate(keys)}
    return sol.t, {b: sum(lines[k] for k in ks) for b, ks in BANDS.items()}


def rise_decay(t, I):
    ip = int(np.argmax(I))
    tp = t[ip]
    # decay from the log-slope between the times where I falls to 30 % and 5 % of the peak
    after = np.arange(ip, t.size)
    i1 = after[np.argmax(I[after] < 0.3 * I[ip])]
    i2 = after[np.argmax(I[after] < 0.05 * I[ip])]
    if i2 <= i1:
        return np.nan, np.nan, tp
    tau_d = -(t[i2] - t[i1]) / np.log(I[i2] / I[i1])
    f = lambda tr: tr * tau_d / (tau_d - tr) * np.log(tau_d / tr) - tp
    try:
        tau_r = brentq(f, 1e-4 * tp, 0.999 * tau_d) if f(0.999 * tau_d) > 0 else np.nan
    except ValueError:
        tau_r = np.nan
    return tau_r, tau_d, tp


def rise_decay_vs_T(model, temps, x0=0.01, extra=None):
    res = {}
    for T in temps:
        t, bands = pulse_curves(model, T, x0=x0, extra=extra)
        res[float(T)] = {b: dict(zip(("tau_r", "tau_d", "t_peak"), map(float, rise_decay(t, I)))) for b, I in bands.items()}
    return res
