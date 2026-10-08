"""Predicción falsable (dirección abierta 1): cociente del decaimiento de la UC contra la fluencia
del pulso de 980 nm, modelo final.

x0 (fracción de Yb excitada) ↔ fluencia F ≈ x0·hν/σ, con σ_Yb(980 nm) = 1,2e-20 cm² (citado por
Langping 2023) y hν = 2,03e-19 J. Válido mientras x0 ≪ 1 (sin saturación del Yb).

    pixi run python noche/prediction_pulse.py
"""

import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).parent))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from dynamics import fit_single_exp, pulse_curves  # noqa: E402
from pint import get_application_registry  # noqa: E402

from upconversion_models.anderson_thermal import AndersonThermal as M  # noqa: E402

u = get_application_registry()
HNU, SIGMA = 2.03e-19, 1.2e-20  # J, cm2
X0 = np.geomspace(1e-3, 0.1, 9)
OUT = pathlib.Path(__file__).parent
res = {}
for hw in (300, 359, 450):
    rows = []
    for x0 in X0:
        tau = {}
        for T in (55.0, 230.0, 353.0, 453.0):
            t, b = pulse_curves(M, T, x0=float(x0), extra={M.phonon_wavenumber: hw / u.cm})
            tau[T] = {k: fit_single_exp(t, b[k]) for k in ("red", "541")}
        rows.append({"x0": float(x0), "F_mJ_cm2": float(x0 * HNU / SIGMA * 1e3),
                     "red_ratio": tau[55.0]["red"] / tau[230.0]["red"], "red_230_us": tau[230.0]["red"] * 1e6,
                     "g_ratio": tau[353.0]["541"] / tau[453.0]["541"], "g_353_us": tau[353.0]["541"] * 1e6})
        print(hw, {k: round(v, 3) for k, v in rows[-1].items()}, flush=True)
    res[hw] = rows
(OUT / "data" / "prediction_pulse.json").write_text(json.dumps(res, indent=1))

fig, axes = plt.subplots(1, 2, figsize=(8, 3.4))
for hw, c in zip(res, ("#bdbdbd", "#737373", "#252525")):
    F = [r["F_mJ_cm2"] for r in res[hw]]
    axes[0].plot(F, [r["red_ratio"] for r in res[hw]], "o-", color=c, lw=2, ms=4, label=f"ħω = {hw} cm⁻¹")
    axes[1].plot(F, [r["g_ratio"] for r in res[hw]], "o-", color=c, lw=2, ms=4, label=f"ħω = {hw} cm⁻¹")
for ax, vals, col, title in ((axes[0], (1.17, 1.22), "#d62728", "rojo: τ(55 K)/τ(230 K)"),
                              (axes[1], (1.60,), "#00a03c", "541 nm: τ(353 K)/τ(453 K)")):
    for v in vals:
        ax.axhline(v, color=col, ls=":", lw=1.5)
    ax.axhline(1, color="#888888", lw=1)
    ax.set_xscale("log")
    ax.set_xlabel("fluencia del pulso de 980 nm (mJ/cm²)")
    ax.set_title(title + " (punteado: Xu 2024)", fontsize=9)
    ax.grid(alpha=0.2)
axes[0].set_ylabel("cociente de τ de la UC")
axes[0].legend(fontsize=7, frameon=False)
fig.suptitle("Predicción: dependencia del cociente de decaimientos con la fluencia (modelo final)", fontsize=10)
fig.tight_layout()
fig.savefig(OUT / "figs" / "final" / "prediction_pulse.png", dpi=150)
