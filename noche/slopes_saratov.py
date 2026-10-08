"""O7b: cambio relativo de la pendiente log-log del verde entre 300 y 320 K a 2 y 10 W/cm²
(Sagaidachnaya y Kochubey 2020: hasta +16 % a alta intensidad y +6 % a baja entre 30 y 50 °C,
1,5–9,4 W/cm²).

    pixi run python noche/slopes_saratov.py iter1_sweep iter2_sweep ...
"""

import json
import pathlib
import sys

DATA = pathlib.Path(__file__).parent / "data"

for name in sys.argv[1:]:
    d = json.loads((DATA / f"{name}.json").read_text())
    row = []
    for P in (2.0, 10.0):
        for band in ("green", "red"):
            a = d["slopes"][f"300.0|{P}"][band]
            b = d["slopes"][f"320.0|{P}"][band]
            row.append(f"{band}@{P:g}: {a:.2f}->{b:.2f} ({(b / a - 1) * 100:+.1f} %)")
    print(f"{name:14s}", " | ".join(row))
