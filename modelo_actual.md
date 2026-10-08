# Modelo actual: `AndersonThermal`

Archivo: `src/upconversion_models/anderson_thermal.py`. Base: `AndersonModel` (Anderson 2013,
Tabla 1, columna "new UC model"). Citas textuales de cada número en `noche/refs.bib`.

## Niveles

| estado | término | energía (cm⁻¹) | g = 2J+1 | fuente |
|---|---|---|---|---|
| Yb1 | ²F7/2 | 0 | — | Anderson (código) |
| Yb2 | ²F5/2 | 10200 | — | Anderson (código) |
| Er1 | ⁴I15/2 | 0 | 16 | Anderson |
| Er2 | ⁴I13/2 | 6500 | 14 | Anderson (código) |
| Er3 | ⁴I11/2 | 10200 | 12 | Anderson (código) |
| Er4 | ⁴I9/2 | 12500 | 10 | Anderson (código) |
| Er5 | ⁴F9/2 | 15000 | 10 | Anderson (código) |
| Er6s | ⁴S3/2 | 18300 | 4 | Er6 de Anderson (código); g de Suta 2025 |
| Er6h | ²H11/2 | 18950 | 12 | ΔE = 650 cm⁻¹, Suta 2025 Fig. 1a; g de Suta 2025 |
| Er7 | ⁴F7/2 | 20500 | 8 | Anderson (código) |
| Er8 | ²H9/2 | 24500 | 10 | Anderson (código) |
| Er9 | ⁴G11/2 | 26100 | 12 | Anderson (código) |

Las energías de Anderson están en el código heredado (`anderson.py`); el paper no las tabula. Sirven
solo para los gaps de las tasas multifonónicas y para las energías de las líneas.

## Parámetros

| parámetro | valor | fuente |
|---|---|---|
| σ_Yb | 1e-20 cm² | sin valor publicado (supuesto del código heredado) |
| k_Yb,rad / k_Yb,NR | 393 / 220 s⁻¹ | Anderson Tabla 1 (k_Yb = 613, k_Yb_rad = 393) |
| k_ET1−3 | 5,9 × 2e-16 cm³/s | Anderson Tabla 1 (k_ET1−3/k_ET3−1 = 5,9; k_ET3−1 = 2e-16) |
| k_ET3−7, k_ET5−8, k_ET6−9 | 1,54e-15, 1,76e-15, 6,07e-15 cm³/s | Anderson Tabla 1 |
| k_ET2−5, k_ET5−2 | 0 | Anderson Tabla 1 |
| k_ET3−1, k_ET7−3, k_ET9−5 | 2e-16, 2,04e-16, 2,84e-16 cm³/s | Anderson Tabla 1 |
| k_NR9, k_NR8, k_NR7, k_NR6, k_NR5, k_NR4, k_NR3 (a T0) | 1,76e6, 43450, 1e6, 26, 0, 22120, 61 s⁻¹ | Anderson Tabla 1 |
| k_R8, k_R6, k_R5, k_R3, k_R2 | 2330, 1510, 2039, 73, 110 s⁻¹ | Anderson Tabla 1 |
| ramificaciones radiativas | 0,40/0,42/0,14/0,04 (8); 0,70/0,25/0,05 (6); 0,90/0,05/0,05 (5); 0,81/0,19 (3) | Anderson, ecuaciones de la p. 39–40 |
| k_CR6, k_CR4, k_UC2 | 2,79e-17, 8,04e-19, 2,31e-17 cm³/s | Anderson Tabla 1 |
| fracción aislada β | 0,05 | Anderson Tabla 1 |
| Yb, Er | 18 %, 2 % | Anderson (muestra) |
| densidad de sitios | 1,38e22 cm⁻³ | código heredado (no publicado por Anderson) |
| T0 | 300 K | supuesto: temperatura de las medidas de Anderson (no informada explícitamente) |
| ħω (multifonónica) | 359 cm⁻¹ | mayor modo Raman de red de β-NaYF4, Dubey 2023 (485 y 628 cm⁻¹ atribuidos a orgánicos); rango de prueba 300–450 (450 = corte, Suta 2025) |
| k_nr(0) (²H11/2 ↔ ⁴S3/2) | 2,29e6 s⁻¹ | Suta 2025, eq. 6 |
| p del acople 6h ↔ 6s | 2 (ħω_eff = 325 cm⁻¹) | Suta 2025 |
| A(²H11/2)/A(⁴S3/2) | 5,52/3 = 1,84 | Suta 2025 (C·g2/g1 = 5,52) |

## Transiciones

- Absorción Yb, decaimiento Yb radiativo y no radiativo.
- ETU Yb→Er: 1→3, 2→5 (0), 3→7, 5→8, 6s→9 y 6h→9 (ambos con k_ET6−9).
- Retrotransferencia Er→Yb: 3→1, 5→2 (0), 7→3, 9→5.
- **Asistidas por fonones (iteración 3)**: 5→8 (ΔE = +700 cm⁻¹), 6s→9 (+2400), 6h→9 (+3050) y la
  retrotransferencia 9→5 (+900), con k(T) = k(T0)·[(1+n(T))/(1+n(T0))]^(ΔE/ħω). Los desajustes salen de
  los baricentros del modelo. 1↔3 y 3↔7 quedan constantes por defecto; con `nr_thermal` = 1 absorben
  o emiten un fonón de 40 y 90 cm⁻¹ (Li 2014) con fracción resonante `r_res` (iteración 4, rechazada).
- Termalización 6h ↔ 6s: baja g_S·k_nr(0)·(1+n)², sube g_H·k_nr(0)·n², n a 325 cm⁻¹.
- Multifonónicas 9→8, 8→7, 7→6h, 6h→5, 6s→5, 5→4 (0), 4→3, 3→2:
  k(T) = k(T0)·[(1+n(T))/(1+n(T0))]^(gap/ħω), y su compañera ascendente por balance detallado
  k↑ = k(T0)·(g_alto/g_bajo)·[n(T)/(1+n(T0))]^(gap/ħω).
- Radiativas de 8, 6h, 6s, 5, 3, 2 con las ramificaciones de Anderson; k_r6s y k_r6h con
  f_S(T0)·k_r6s + f_H(T0)·k_r6h = k_R6 y k_r6h = 1,84·k_r6s.
- Relajación cruzada 6s,6h + 1 → 3 + 2 (k_CR6) **asistida por fonones** (iteración 2): desajuste
  1600 cm⁻¹ (6s) y 2250 cm⁻¹ (6h) emitido como fonones, k(T) = k_CR6(T0)·[(1+n(T))/(1+n(T0))]^(ΔE/ħω)
  (Langping 2023 para el mecanismo; la ley es la misma de las multifonónicas).
- 4 + 1 → 2 + 2 (k_CR4) y ETU Er–Er 2 + 2 → 4 + 1 (k_UC2), constantes.

## Interruptores

| parámetro | efecto | valor que apaga |
|---|---|---|
| `thermal` | evalúa todas las tasas dependientes de T a T0 | 0 |
| `detailed_balance` | compañeros ascendentes de las multifonónicas | 0 |
| `k_nr0` | acople 6h ↔ 6s | (grande = termalización instantánea) |
| `cr_thermal` | CR6 asistida por fonones (iteración 2) | 0 |
| `et_thermal` | transferencias Yb↔Er asistidas por fonones (iteración 3) | 0 |
| `nr_thermal` | 1↔3 y 3↔7 con un fonón (Li 2014) y fracción resonante `r_res` (iteración 4, mal condicionada) | 0 (**por defecto**) |
