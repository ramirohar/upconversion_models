"""Tabla de contraste de varios barridos contra los números publicados del mapa.

    pixi run python noche/compare.py iter1_sweep iter2_sweep ...
"""

import json
import pathlib
import sys

import numpy as np

DATA = pathlib.Path(__file__).parent / "data"


def load(name):
    return json.loads((DATA / f"{name}.json").read_text())


def ser(d, key, P):
    return np.array([d["obs"][f"{T}|{P}"][key] for T in d["temps"]])


def at(d, y, T):
    return float(np.interp(T, d["temps"], y))


def tau(d, which):
    return np.array([d["tau"][str(t)][which] for t in d["temps"]]) * 1e6


def metrics(d):
    T = np.array(d["temps"])
    m = {}
    for P in (0.1, 0.8, 290.0):
        g = ser(d, "green", P)
        vis = g + ser(d, "red", P) + ser(d, "blue", P)
        m[f"T_max verde @{P}"] = T[g.argmax()]
        m[f"verde máx/300K @{P}"] = g.max() / at(d, g, 300)
        m[f"verde 10K/300K @{P}"] = g[0] / at(d, g, 300)
        m[f"T_max visible @{P}"] = T[vis.argmax()]
        m[f"visible 160K/300K @{P}"] = at(d, vis, 160) / at(d, vis, 300)
        m[f"visible 40K/300K @{P}"] = at(d, vis, 40) / at(d, vis, 300)
        for band in ("red", "blue"):
            y = ser(d, band, P)
            m[f"T_max {band} @{P}"] = T[y.argmax()]
            m[f"{band} 10K/300K @{P}"] = y[0] / at(d, y, 300)
        rg = ser(d, "RG", P)
        m[f"R/G 483K/303K @{P}"] = at(d, rg, 483) / at(d, rg, 303)
        m[f"R/G 400K/10K @{P}"] = at(d, rg, 400) / rg[0]
    tg = tau(d, "tau_green")
    m["tau verde 40K (µs)"] = at(d, tg, 40)
    m["tau verde 300K (µs)"] = at(d, tg, 300)
    m["tau verde 40K/300K"] = at(d, tg, 40) / at(d, tg, 300)
    m["tau verde 353K/453K"] = at(d, tg, 353) / at(d, tg, 453)
    s = lambda k, P, Tq: at(d, np.array([d["slopes"][f"{t}|{P}"][k] for t in d["temps"]]), Tq)
    for k in ("H", "S", "red"):
        m[f"n_{k} 303K @0.8"] = s(k, 0.8, 303)
        m[f"n_{k} 483K @0.8"] = s(k, 0.8, 483)
    return m


LIT = {
    "visible 160K/300K @0.1": "2.2 (Langping, máx ~160 K)",
    "T_max visible @290.0": "~100 bulk / ~150 nano (Yu)",
    "T_max verde @290.0": "~100 bulk / ~150 nano (Yu)",
    "T_max red @290.0": "~100 bulk / ~150 nano (Yu)",
    "T_max blue @290.0": "~100 bulk / ~150 nano (Yu)",
    "tau verde 40K/300K": "2.59 (Langping 430/166)",
    "tau verde 353K/453K": "1.60 (Xu 226/141)",
    "R/G 483K/303K @0.8": "> 1 (Xu: rojo pasa a ser el pico mayor)",
    "R/G 400K/10K @290.0": "> 1 (Yu bulk, monótono)",
    "n_H 303K @0.8": "1.53 (Xu)", "n_S 303K @0.8": "1.22 (Xu)", "n_red 303K @0.8": "1.30 (Xu)",
    "n_H 483K @0.8": "1.82 (Xu)", "n_S 483K @0.8": "1.73 (Xu)", "n_red 483K @0.8": "1.69 (Xu)",
}

if __name__ == "__main__":
    names = sys.argv[1:]
    ms = [metrics(load(n)) for n in names]
    keys = list(ms[0])
    w = max(map(len, keys))
    print(" " * w, *[f"{n:>13s}" for n in names], "  literatura")
    for k in keys:
        print(f"{k:<{w}}", *[f"{m[k]:13.3g}" for m in ms], "  " + LIT.get(k, ""))
