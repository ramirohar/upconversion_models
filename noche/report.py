"""Resumen numérico y figuras de un barrido estándar (salida de sweep.run).

    pixi run python noche/report.py data/iter1_sweep.json figs/iter1
"""

from __future__ import annotations

import json
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

HERE = pathlib.Path(__file__).parent
GREEN, RED, BLUE, INK, MUTED = "#00a03c", "#d62728", "#1f5fd6", "#222222", "#888888"
# literatura (refs.bib): observables con número en el texto
LIT_TAU_GREEN = [(40, 430, "Langping23"), (300, 166, "Langping23"), (353, 226, "Xu24"), (453, 141, "Xu24")]
LIT_TAU_RED = [(230, 497, "Xu24"), (55, 583, "Xu24")]
LIT_DE = {"Zhou13": 752, "Geit17": 714, "Tong15": 713, "Dubey23": 816, "Xu24": 817}


def load(path):
    return json.loads(pathlib.Path(path).read_text())


def series(d, name, P):
    return np.array([d["obs"][f"{T}|{P}"][name] for T in d["temps"]])


def slope(d, name, P):
    return np.array([d["slopes"][f"{T}|{P}"][name] for T in d["temps"]])


def power_colors(powers):
    cm = plt.get_cmap("Greys")
    return {P: cm(0.35 + 0.65 * i / max(len(powers) - 1, 1)) for i, P in enumerate(powers)}


def summary(d):
    T = np.array(d["temps"])
    lines = []
    lines.append("ΔE_eff (cm-1) [FIR_pure / FIR_band]:")
    for P in d["powers"]:
        row = []
        for w in ("Zhou 160-300", "Tong 303-573", "Xu 353-453", "Geit 300-700"):
            a = d["dE_eff"][f"{P}|FIR_pure|{w}"]["dE"]
            b = d["dE_eff"][f"{P}|FIR_band|{w}"]["dE"]
            row.append(f"{w}: {a:6.0f}/{b:6.0f}")
        A = d["dE_eff"][f"{P}|FIR_pure|Geit 300-700"]["A"]
        lines.append(f"  P={P:>6}: " + " | ".join(row) + f" | A(pure, 300-700)={A:.2f}")
    lines.append("T del máximo (K) de verde / rojo / azul / visible, y cociente máx/300 K del verde:")
    i300 = int(np.argmin(np.abs(T - 300)))
    for P in d["powers"]:
        g, r, b = series(d, "green", P), series(d, "red", P), series(d, "blue", P)
        vis = g + r + b
        lines.append(
            f"  P={P:>6}: {T[g.argmax()]:.0f} / {T[r.argmax()]:.0f} / {T[b.argmax()]:.0f} / {T[vis.argmax()]:.0f}"
            f"   verde(10K)/verde(300K)={g[0] / g[i300]:.2f}  verde(500K)/verde(300K)={g[np.argmin(np.abs(T - 500))] / g[i300]:.2f}"
        )
    lines.append("Rojo/verde a 100 / 300 / 500 K:")
    for P in d["powers"]:
        rg = series(d, "RG", P)
        idx = [int(np.argmin(np.abs(T - t))) for t in (100, 300, 500)]
        lines.append(f"  P={P:>6}: " + " / ".join(f"{rg[i]:.3f}" for i in idx))
    lines.append("Pendientes log-log (verde S 541 / verde H 520 / rojo) a 300 y 480 K:")
    for P in d["powers"]:
        idx = [int(np.argmin(np.abs(T - t))) for t in (300, 480)]
        s = [f"{slope(d, 'S', P)[i]:.2f}/{slope(d, 'H', P)[i]:.2f}/{slope(d, 'red', P)[i]:.2f}" for i in idx]
        lines.append(f"  P={P:>6}: 300 K {s[0]}  480 K {s[1]}")
    tg = np.array([d["tau"][str(t)]["tau_green"] for t in d["temps"]]) * 1e6
    tr = np.array([d["tau"][str(t)]["tau_red"] for t in d["temps"]]) * 1e6
    lines.append("τ verde / τ rojo (µs), excitación directa débil:")
    for t in (40, 60, 220, 300, 360, 440):
        i = int(np.argmin(np.abs(T - t)))
        lines.append(f"  T={T[i]:.0f} K: {tg[i]:.0f} / {tr[i]:.0f}")
    return "\n".join(lines)


def figures(d, outdir, label=""):
    outdir = pathlib.Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    T = np.array(d["temps"])
    powers = d["powers"]
    pc = power_colors(powers)

    # 1. ln FIR contra 1/T
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6), sharey=True)
    for ax, fir in zip(axes, ("FIR_pure", "FIR_band")):
        for P in powers:
            ax.plot(1000 / T, np.log(series(d, fir, P)), color=pc[P], lw=2, label=f"{P:g} W/cm²")
        boltz = np.log(5.52) - 650 / 0.6950348 / T
        ax.plot(1000 / T, boltz, ls="--", color=GREEN, lw=1.5, label="Boltzmann 650 cm⁻¹")
        ax.set_xlabel("1000 / T (K⁻¹)")
        ax.set_title(fir.replace("_", " "), fontsize=10)
        ax.set_xlim(1.2, 12)
        ax.set_ylim(-12, 1)
        ax.grid(alpha=0.2)
    axes[0].set_ylabel("ln FIR (520/541)")
    axes[1].legend(fontsize=7, frameon=False)
    fig.suptitle(f"FIR {label}", fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / "fir.png", dpi=150)
    plt.close(fig)

    # 2. intensidades contra T (normalizadas a 300 K), una columna por potencia
    sel = [P for P in powers if P in (0.1, 0.8, 10.0, 290.0)]
    fig, axes = plt.subplots(1, len(sel), figsize=(3 * len(sel), 3.2), sharey=True)
    i300 = int(np.argmin(np.abs(T - 300)))
    for ax, P in zip(np.atleast_1d(axes), sel):
        for name, col in (("green", GREEN), ("red", RED), ("blue", BLUE)):
            y = series(d, name, P)
            ax.plot(T, y / y[i300], color=col, lw=2, label=name)
        ax.set_yscale("log")
        ax.set_title(f"{P:g} W/cm²", fontsize=10)
        ax.set_xlabel("T (K)")
        ax.grid(alpha=0.2)
        ax.axvspan(100, 160, color=MUTED, alpha=0.15, lw=0)
    np.atleast_1d(axes)[0].set_ylabel("I(T) / I(300 K)")
    np.atleast_1d(axes)[0].legend(fontsize=7, frameon=False)
    fig.suptitle(f"Intensidades {label} (franja gris: máximos publicados 100–160 K)", fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / "intensities.png", dpi=150)
    plt.close(fig)

    # 3. rojo/verde contra T
    fig, ax = plt.subplots(figsize=(4.5, 3.4))
    for P in powers:
        ax.plot(T, series(d, "RG", P), color=pc[P], lw=2, label=f"{P:g} W/cm²")
    ax.set_yscale("log")
    ax.set_xlabel("T (K)")
    ax.set_ylabel("rojo / verde")
    ax.legend(fontsize=7, frameon=False)
    ax.grid(alpha=0.2)
    ax.set_title(f"Rojo/verde {label}", fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / "red_green.png", dpi=150)
    plt.close(fig)

    # 4. pendientes log-log a 0.8 W/cm2 (Xu 2024)
    P = 0.8 if 0.8 in powers else powers[0]
    fig, ax = plt.subplots(figsize=(4.5, 3.4))
    for name, col, lab in (("H", GREEN, "520 (2H11/2)"), ("S", "#7a9a00", "541 (4S3/2)"), ("red", RED, "653 (rojo)")):
        ax.plot(T, slope(d, name, P), color=col, lw=2, label=lab)
    for x, ys in ((303, (1.53, 1.22, 1.30)), (483, (1.82, 1.73, 1.69))):
        for yv, col in zip(ys, (GREEN, "#7a9a00", RED)):
            ax.plot(x, yv, "o", ms=8, mfc="white", mec=col, mew=2)
    ax.set_xlabel("T (K)")
    ax.set_ylabel("n = d ln I / d ln P")
    ax.set_title(f"Pendientes a {P:g} W/cm² (círculos: Xu 2024) {label}", fontsize=9)
    ax.legend(fontsize=7, frameon=False)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(outdir / "slopes.png", dpi=150)
    plt.close(fig)

    # 5. tiempos de vida
    tg = np.array([d["tau"][str(t)]["tau_green"] for t in d["temps"]]) * 1e6
    tr = np.array([d["tau"][str(t)]["tau_red"] for t in d["temps"]]) * 1e6
    fig, ax = plt.subplots(figsize=(4.5, 3.4))
    ax.plot(T, tg, color=GREEN, lw=2, label="τ 4S3/2 (modelo)")
    ax.plot(T, tr, color=RED, lw=2, label="τ 4F9/2 (modelo)")
    for x, y, src in LIT_TAU_GREEN:
        ax.plot(x, y, "o", ms=8, mfc="white", mec=GREEN, mew=2)
        ax.annotate(src, (x, y), fontsize=6, color=MUTED, xytext=(3, 3), textcoords="offset points")
    for x, y, src in LIT_TAU_RED:
        ax.plot(x, y, "s", ms=8, mfc="white", mec=RED, mew=2)
        ax.annotate(src, (x, y), fontsize=6, color=MUTED, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel("T (K)")
    ax.set_ylabel("τ (µs)")
    ax.legend(fontsize=7, frameon=False)
    ax.grid(alpha=0.2)
    ax.set_title(f"Vidas medias {label}", fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / "lifetimes.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    d = load(HERE / sys.argv[1])
    out = HERE / sys.argv[2]
    print(summary(d))
    figures(d, out, label=sys.argv[3] if len(sys.argv) > 3 else "")
