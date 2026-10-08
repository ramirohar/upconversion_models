"""Figuras clave de la sesión (modelo final = iteración 3, más las comparaciones relevantes).

    pixi run python noche/figs_final.py

Salida en noche/figs/final/:
- delta_e.png: ΔE efectivo de ln FIR contra potencia (FIR puro y con 2H9/2→4I13/2 en la banda de
  4S3/2) frente a los valores publicados (O2).
- red_green.png: rojo/verde contra T a 0,8 y 290 W/cm², iteraciones 1, 2 y 3 (O6).
- intensities.png: verde, rojo y azul contra T a 0,1 y 290 W/cm², modelo final (O4).
- uc_decay.png: decaimiento de la UC tras pulso (rojo y 541 nm) contra la intensidad del pulso,
  modelo final, frente a Xu 2024 (iteración 5).
"""

import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = pathlib.Path(__file__).parent
DATA = HERE / "data"
OUT = HERE / "figs" / "final"
OUT.mkdir(parents=True, exist_ok=True)
GREEN, RED, BLUE, INK, MUTED = "#00a03c", "#d62728", "#1f5fd6", "#222222", "#888888"
SEQ = ["#bdbdbd", "#737373", "#252525"]


def load(n):
    return json.loads((DATA / f"{n}.json").read_text())


def ser(d, k, P):
    return np.array([d["obs"][f"{T}|{P}"][k] for T in d["temps"]])


fin = load("iter3_sweep")
T = np.array(fin["temps"])

# 1. ΔE efectivo contra potencia
fig, ax = plt.subplots(figsize=(5.2, 3.6))
P = np.array(fin["powers"])
for fir, ls, lab in (("FIR_pure", "-", "modelo, FIR puro"), ("FIR_band", "--", "modelo, banda 541 con 2H9/2→4I13/2")):
    y = [fin["dE_eff"][f"{p}|{fir}|Tong 303-573"]["dE"] for p in P]
    ax.plot(P, y, ls, color=INK, lw=2, marker="o", ms=5, label=lab)
lit = [("Zhou13 (752)", 752, None), ("Geit17, Tong15 (713–714)", 713.5, None), ("Dubey23, Xu24 (816–817)", 816.5, None)]
for name, v, p in lit:
    ax.axhline(v, color=MUTED, lw=1, alpha=0.6)
    ax.annotate(name, (0.12, v), fontsize=7, color=MUTED, xytext=(0, 2), textcoords="offset points")
ax.axhline(650, color=GREEN, lw=1.5, ls=":")
ax.annotate("650 cm⁻¹ (Suta, excitación 77 K)", (0.12, 650), fontsize=7, color=GREEN, xytext=(0, -9), textcoords="offset points")
ax.set_xscale("log")
ax.set_xlabel("P (W/cm²)")
ax.set_ylabel("ΔE efectivo, ajuste 303–573 K (cm⁻¹)")
ax.set_title("O2: pendiente de ln FIR (modelo final)", fontsize=10)
ax.legend(fontsize=7, frameon=False, loc="center right")
fig.tight_layout()
fig.savefig(OUT / "delta_e.png", dpi=150)
plt.close(fig)

# 2. rojo/verde contra T
fig, axes = plt.subplots(1, 2, figsize=(8, 3.4), sharey=True)
for ax, p in zip(axes, (0.8, 290.0)):
    for (n, lab), c in zip((("iter1_sweep", "iter. 1"), ("iter2_sweep", "iter. 2"), ("iter3_sweep", "iter. 3 (final)")), SEQ):
        d = load(n)
        ax.plot(d["temps"], ser(d, "RG", p), color=c, lw=2, label=lab)
    ax.set_yscale("log")
    ax.set_xlabel("T (K)")
    ax.set_title(f"{p:g} W/cm² ({'Xu 2024: R/G sube' if p == 0.8 else 'Yu 2014: R/G sube'})", fontsize=9)
    ax.grid(alpha=0.2)
axes[0].set_ylabel("rojo / verde")
axes[0].legend(fontsize=7, frameon=False)
fig.suptitle("O6: cociente rojo/verde", fontsize=10)
fig.tight_layout()
fig.savefig(OUT / "red_green.png", dpi=150)
plt.close(fig)

# 3. intensidades del modelo final
fig, axes = plt.subplots(1, 2, figsize=(8, 3.4), sharey=True)
i300 = int(np.argmin(np.abs(T - 300)))
for ax, p in zip(axes, (0.1, 290.0)):
    for k, c in (("green", GREEN), ("red", RED), ("blue", BLUE)):
        y = ser(fin, k, p)
        ax.plot(T, y / y[i300], color=c, lw=2, label=k)
    ax.axvspan(100, 160, color=MUTED, alpha=0.15, lw=0)
    ax.set_yscale("log")
    ax.set_xlabel("T (K)")
    ax.set_title(f"{p:g} W/cm²", fontsize=9)
    ax.grid(alpha=0.2)
axes[0].set_ylabel("I(T) / I(300 K)")
axes[0].legend(fontsize=7, frameon=False)
fig.suptitle("O4: intensidades del modelo final (gris: máximos publicados, Yu y Langping)", fontsize=10)
fig.tight_layout()
fig.savefig(OUT / "intensities.png", dpi=150)
plt.close(fig)

# 4. decaimiento de la UC contra intensidad del pulso
dec = load("iter5_decays")
fig, axes = plt.subplots(1, 2, figsize=(8, 3.4))
x0s = [0.003, 0.01, 0.1]
for hw, c in zip((300, 359, 450), SEQ):
    red = [dec[f"final (iter. 3)|hw={hw}|x0={x}"] for x in x0s]
    axes[0].plot(x0s, [r["55.0"]["red"]["xu_fit"] / r["230.0"]["red"]["xu_fit"] for r in red], "o-", color=c, lw=2, label=f"ħω = {hw} cm⁻¹")
    axes[1].plot(x0s, [r["353.0"]["541"]["xu_fit"] / r["453.0"]["541"]["xu_fit"] for r in red], "o-", color=c, lw=2, label=f"ħω = {hw} cm⁻¹")
for ax, v, lab in ((axes[0], (1.17, 1.22), "rojo: τ(55 K)/τ(230 K)"), (axes[1], (1.60,), "541 nm: τ(353 K)/τ(453 K)")):
    for vv in v:
        ax.axhline(vv, color=RED if ax is axes[0] else GREEN, lw=1.5, ls=":")
    ax.axhline(1, color=MUTED, lw=1)
    ax.set_xscale("log")
    ax.set_xlabel("fracción de Yb excitada por el pulso, x0")
    ax.set_title(lab + " (punteado: Xu 2024)", fontsize=9)
    ax.grid(alpha=0.2)
axes[0].set_ylabel("cociente de τ de la UC")
axes[0].legend(fontsize=7, frameon=False)
fig.suptitle("Iteración 5: decaimiento de la UC tras pulso de 980 nm (modelo final)", fontsize=10)
fig.tight_layout()
fig.savefig(OUT / "uc_decay.png", dpi=150)
plt.close(fig)
print("ok")
