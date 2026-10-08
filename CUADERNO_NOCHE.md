# Cuaderno de la sesión nocturna 2026-10-08

Rama `noche/2026-10-08` (desde `unfolding` @ 4eed439). Plan: `PLAN_NOCHE.md`.
Fuentes y citas textuales: `noche/refs.bib`. Scripts, datos y figuras: `noche/`.

## Bitácora de reloj

- 02:09 inicio. Fase 0 (instalación, corrida de AndersonModel, línea base, regresión) hasta ~02:50,
  commit 960a350.
- Entre ~02:50 y ~11:30 la máquina entró en reposo (pmset: "Maintenance Sleep" con "DarkWake" cada
  15 min, a batería); la sesión quedó prácticamente detenida.
- 11:36 se reanuda con la máquina despierta; el presupuesto de ~8 h se reinicia desde aquí (cierre ~19:30).

## Fase 0: base

### 0.1 Entorno y corrida
`pixi install` sin problemas (pixi 0.75). Se agregó `pytest` al entorno. `AndersonModel` corre;
un estado estacionario tarda ~0,2 s (barrido de 180 puntos: ~40 s, muy por debajo del tope de 20 min).

### 0.2 Er6 es el par verde
Anderson 2013 describe el nivel emisor verde como "the thermalized emitting-state manifold,
2H11/2,4S3/2" y el ETU nuevo como "out of the green-emitting 2H11/2,4S3/2 states" con la flecha
"kET6−9" (refs.bib: anderson2013 quote2, quote3). En el código `Er6` está a 18300 cm⁻¹, que es la
posición de 4S3/2: el nivel 6 de Anderson es el par 2H11/2 + 4S3/2 tratado como un solo estado.

Nota: `AndersonModel` escribe las energías con `energy_factor = (h·c).magnitude · eV`, que da números
sin sentido físico (10200 "cm⁻¹" → 2,03e-21 eV). No afecta a Anderson (sin dependencia térmica) pero el
modelo con T usa `hc / u.cm`.

### 0.3 Solver, tolerancias y línea base
- `src/upconversion_models/steady.py`: integración BDF (rtol 1e-10, atol 1e-20) hasta que
  max|dy/dt|/max(y) < 1e-9 s⁻¹, seguida de un pulido de Newton sobre f(y) = 0 con las leyes de
  conservación (espacio nulo izquierdo del jacobiano: Yb total y Er total). Residuo relativo < 1e-10.
- Contra `jablonski_patch.SteadyState` (criterio SmallDerivatives(atol=1e-12, rtol=1e-6)) la diferencia
  máxima es 2e-8 relativa (0,1–1000 W/cm²), así que ambos sirven, pero la regresión a 1e-6 usa el pulido.
- Línea base: `noche/fase0_baseline.py` → `noche/data/baseline_anderson.json`, grilla
  0,1–1000 W/cm² (9 puntos), poblaciones y tasas `line_*`.
- Regresión: `tests/test_regression_anderson.py` (rtol 1e-6, 9 casos, pasan). Correr con
  `pixi run python -m pytest -q tests`.

### 0.4 Mapa de observaciones

Todas sobre β-NaYF4:Yb,Er con ~18–20 % Yb y 2 % Er salvo indicación. "Calib." = ¿usado en la
calibración de Anderson 2013? Anderson calibró solo a temperatura ambiente con curvas de decaimiento
y cocientes de intensidad bajo excitación pulsada (943 nm, 2–200 mJ/cm²) y con la dependencia del
cociente verde/rojo con la concentración de Er (refs: anderson2013 quote1, quote5). Ninguna
observación en función de T se usó para calibrar.

| # | observable | sistema | condiciones | valor o tendencia | figura/tabla | acceso | calib. |
|---|---|---|---|---|---|---|---|
| O1 | FIR = I(2H11/2)/I(4S3/2) contra T | varios (ver O2) | 980 nm cw | crece monótono; ln FIR lineal en 1/T | Yu14 Fig 5c; Zhou13 Fig 3b; Geit17 Fig 3b, 6b; Tong15 Fig 6; Dubey23 Fig 7b; Xu24 Fig 6a | sí | no |
| O2 | ΔE efectivo de la pendiente de ln FIR | micro 20Yb2Er (Zhou); 20 nm 18Yb2Er (Geit.); micro/nano 20Yb2Er (Tong); 29 nm (Dubey); micro 20Yb2Er (Xu) | 160–300 K (Zhou, 29 mW); 300–900 K (Geit.); 303–573 K (Tong, 2 W/cm²); 300–650 K (Dubey); 353–453 K (Xu, 0,8 W/cm²) | 752; 714 / 716; 713; 816; 817 cm⁻¹, contra 650 ± 10 cm⁻¹ de excitación a 77 K (Suta) | Zhou13 eq. 4; Geit17 Fig 3b, 6b; Tong15 Fig 6a; Dubey23 Fig 7b; Xu24 Fig 6a | sí | no |
| O3 | prefactor C·g_H/g_S del ajuste de FIR | idem | idem | 8,06 (Zhou); 9,13 (Tong); 5,52 (Suta, con término de offset) | Zhou13 eq. 4; Tong15 Fig 6a | sí | no |
| O4 | intensidad UC total contra T por debajo de ambiente | bulk y 25/45 nm 18Yb2Er (Yu); core-shell 20 nm 20Yb2Er (Langping) | 290 W/cm² (Yu); 0,1 W/cm² (Langping) | no monótona: máximo ~100 K (bulk) y ~150 K (nano) en Yu; máximo ~160 K, 2,2× el valor de 300 K en Langping; todas las bandas (azul, verde, rojo) en Yu | Yu14 Fig 4; Langping23 Fig 2b,c | sí | no |
| O4b | ídem, Suyver 2005 (bulk, 5–200 K) | bulk | — | aumento de 5 a 50 K (solo vía cita de Li 2014) | — | inaccesible | no |
| O5 | intensidad UC sobre ambiente | >32 nm (Li); 29 nm (Dubey); micro (Xu) | 300–650 K | decrece con T (apagado térmico); en <32 nm aumenta 25–85 °C (Li, efecto de tamaño) | Li14 Fig 5a; Dubey23 Fig 7a; Xu24 Fig 2b | sí | no |
| O6 | cociente rojo/verde contra T | bulk 18Yb2Er (Yu); micro 20Yb2Er (Xu) | 290 W/cm² (Yu); 0,8 W/cm² (Xu) | crece monótono con T (bulk); a 483 K el pico rojo 653 nm supera al verde 541 nm, al revés que a 303 K | Yu14 Fig 5d; Xu24 Fig 3a,d | sí | no |
| O7 | pendientes log-log I ∝ P^n contra T | micro 20Yb2Er (Xu) | 303 y 483 K | 520/541/653 nm: 1,53/1,22/1,30 a 303 K y 1,82/1,73/1,69 a 483 K (crecen con T) | Xu24 Fig 3c,f | sí | no |
| O8 | pendientes log-log a T ambiente | Yu: 4 tamaños; Zhou: micro | 300 K | verde 1,38–1,92 (Yu); 1,95 a baja densidad y decreciente a alta (Zhou) | Yu14 Fig 3; Zhou13 Fig 2 | sí | no (Anderson usó pulsos) |
| O9 | vida media de 4S3/2 contra T | core-shell 20Yb2Er (Langping, exc. 488 nm); micro (Xu, exc. 980 nm) | 40–300 K; 353–453 K | 430 → 166 µs (40 → 300 K); 226 → 141 µs (353 → 453 K) | Langping23 Fig 3e; Xu24 Fig 5b | sí | no (Anderson: decaimientos a T ambiente) |
| O10 | vida media de 4F9/2 contra T | micro (Xu) | 230 → 55 K | 497 → 583 µs (653 nm), 496 → 605 µs (661 nm) | Xu24 Fig 5 | sí | no |
| O11 | tiempos de subida UC contra T | 25, 45 nm y bulk (Yu) | 10–400 K, pulsado 980 nm | suben de 10 a 30 K y bajan por encima de 30 K, en todas las transiciones | Yu14 Fig 7a,c,e | sí | no |
| O12 | emisión azul (2H9/2, 407 nm) contra T | micro (Xu) | 0,8 W/cm² | decrece con T de 55 a 195 K y de 353 a 453 K | Xu24 Fig 4c,d | sí | no |
| O13 | contaminación de la banda verde por 2H9/2 → 4I13/2 (557 nm) | micro (Zhou) | alta potencia | aparece a alta potencia (proceso de tres fotones) | Zhou13 Fig 5a | sí | no |
| O14 | calentamiento por el láser | micro (Zhou) | 500 mW, 180/240/300 K | ΔT = 38/30/18 K | Zhou13 Fig 5b | sí | no (no es observable del modelo) |
| O15 | cociente verde/rojo contra energía de pulso, curvas de decaimiento a 300 K | micro 18Yb2Er | 943 nm pulsado | — | Anderson13 Fig 3 | sí | **sí** |

Parámetros de literatura (no son datos de contraste):
ΔE(2H11/2–4S3/2) = 650 ± 10 cm⁻¹, g_H = 12, g_S = 4, C = A_H/A_S = 5,52/3 = 1,84,
k_nr(0) = 2,29 µs⁻¹ con p = 2 (ħω_eff = 325 cm⁻¹) (Suta 2025); fonones de red de β-NaYF4 hasta
359 cm⁻¹ en Raman (Dubey 2023; 485 y 628 cm⁻¹ atribuidos a orgánicos), corte 450 cm⁻¹ (Suta 2025).

Observaciones repetidas por varios grupos y no usadas en la calibración: O1, O2 (cinco grupos),
O4 (dos grupos con acceso, un tercero vía cita), O5, O6 (dos grupos), O9 (dos grupos).
Inaccesibles (Unpaywall is_oa = False): Suyver 2005 y 2006 (J. Lumin.), Renero-Lecuna 2011
(Chem. Mater.), Pollnau 2000 (PRB).

## Iteración 1 (fija): desdoblamiento y dependencia térmica

### Modelo
`src/upconversion_models/anderson_thermal.py` → `AndersonThermal` (ver `modelo_actual.md`).
- Er6 → Er6s (4S3/2, g = 4, 18300 cm⁻¹) y Er6h (2H11/2, g = 12, +650 cm⁻¹; Suta 2025 Fig. 1a).
- Acople 6h ↔ 6s de Suta 2025 eq. 6: baja = g_S·k_nr(0)·(1+n)², sube = g_H·k_nr(0)·n², n a
  ħω = 325 cm⁻¹, k_nr(0) = 2,29 µs⁻¹ → sube/baja = 3·e^(−650 cm⁻¹/kT).
- A_H/A_S = 5,52/3 (Suta 2025), con f_S(T0)·k_r6s + f_H(T0)·k_r6h = 1510 s⁻¹ (Anderson).
  k_nr6, k_ET6−9 y k_CR6 de Anderson, iguales para ambos subniveles. A(T) constante.
- Multifonónicas k(T) = k(T0)·[(1+n(T))/(1+n(T0))]^p, p = gap/ħω, ħω = 359 cm⁻¹ (mayor modo Raman
  de red de β-NaYF4, Dubey 2023), con el compañero ascendente por balance detallado (g = 2J+1).
- Interruptores del mismo modelo: `thermal` (0 → todas las tasas a T0), `detailed_balance`
  (0 → sin compañeros ascendentes), `k_nr0` (acople 6h–6s).

### Tests (`tests/test_anderson_thermal.py`, 21 tests con los de la Fase 0, pasan)
- termalización rápida (k_nr0 = 1e10 s⁻¹), T = T0 y `detailed_balance` = 0: se recupera la línea base
  a rtol 1e-6 (0,1, 10 y 1000 W/cm²);
- `thermal` = 0, `detailed_balance` = 0, termalización rápida: se recupera Anderson a 20, 120, 500 y
  700 K (1 y 100 W/cm²);
- partición de Boltzmann 6h/6s con termalización rápida a 300 y 600 K (rtol 1e-4). A 100 K la
  alimentación de 2H11/2 desde 4F7/2 ya rompe el equilibrio aun con k_nr0 = 1e10 s⁻¹.

Nota: con los compañeros ascendentes encendidos el modelo por defecto no es idéntico a Anderson a T0:
Anderson ajustó tasas efectivas sin ellos. Se cuantifica en los resultados.

Cambio numérico: `steady.py` pasó de BDF a LSODA (BDF tardaba > 20 s con acoples de 1e9 s⁻¹), con
umbral de integración 1e-7 s⁻¹ y t_end = 10 s, y el pulido de Newton ahora escala filas y columnas y
toma las leyes de conservación por especie a partir de los nombres (el SVD confundía modos lentos con
modos nulos cuando las tasas abarcan 10 décadas). La regresión de la Fase 0 sigue pasando a 1e-6.

### Pre-registro (antes del barrido estándar)
Antes de este pre-registro corrí una prueba de humo: 10 W/cm² a 20, 100, 300 y 500 K (verde total,
rojo, azul, rad82, poblaciones de Er6). Lo que vi ahí condiciona las expectativas E2, E3 y E5, y lo
declaro.

Barrido: `noche/sweep.py` + `noche/iter1_run.py`; T = 10–90 K cada 10 y 100–700 K cada 20;
P = 0,1, 0,8, 2, 10, 100 y 290 W/cm².

| fila | expectativa | signo/magnitud esperados |
|---|---|---|
| O1 | FIR_pure crece con T; por debajo de ~150–200 K queda por encima de Boltzmann (alimentación no térmica de 2H11/2; es el T_on conocido, no cuenta) | — |
| O2 | ΔE_eff de FIR_pure en ventanas ≥ 300 K ≈ 650 ± 20 cm⁻¹, es decir, **no** reproduce los 713–817 cm⁻¹ publicados. Con FIR_band (rad82 dentro de la banda de 4S3/2) ΔE_eff > 650 y creciente con P: unos +40 cm⁻¹ a 10 W/cm² (estimado de la prueba de humo), pocos cm⁻¹ a 0,8 W/cm², más de +100 cm⁻¹ a 290 W/cm². Como esto ya lo anticipo, no podrá contar como emergente salvo que la magnitud sea muy distinta | ver texto |
| O3 | prefactor de FIR_pure ≈ 3·1,84 = 5,5 en ventanas altas | — |
| O4 | verde, rojo y azul crecen monótonamente al enfriar (sin máximo a 100–160 K) a todas las potencias: **contradice** Yu 2014 y Langping 2023 | monótono |
| O5 | verde y rojo bajan con T por encima de 300 K: reproduce | — |
| O6 | rojo/verde **baja** con T (prueba de humo: 1,40 / 0,71 / 0,37 a 100 / 300 / 500 K): **contradice** Yu 2014 (bulk) y Xu 2024 | signo opuesto |
| O7 | pendientes log-log crecen con T (más pérdidas lineales de 4I11/2 y 4I9/2): reproduce el signo de Xu 2024; magnitud incierta | + |
| O9 | τ(4S3/2) por excitación directa baja poco con T: ~10 % entre 40 y 300 K y ~15 % entre 353 y 453 K, contra ×2,6 (Langping) y ×1,6 (Xu); valor absoluto ~650 µs (Anderson) contra 166–430 µs: **subestima** | débil |
| O10 | τ(4F9/2) constante (k_nr5 = 0) ≈ 490 µs; Xu ve 497 → 583–605 µs al enfriar 230 → 55 K: **no reproduce** | 0 |
| O11 | tiempos de subida: no se evalúan en esta iteración | — |
| O12 | azul baja con T: reproduce | — |

Sorpresas posibles: un máximo de la intensidad verde contra T a alguna potencia; rojo/verde creciente
con T; ΔE_eff de FIR_pure lejos de 650 cm⁻¹ por encima de 300 K.

### Resultados (barrido `noche/data/iter1_sweep.json`, 217 s; figuras `noche/figs/iter1/`)
Generados con `pixi run python noche/iter1_run.py` y `pixi run python noche/report.py data/iter1_sweep.json figs/iter1`.

**A T0 el modelo por defecto se aparta de Anderson** en ≤ 0,24 % en el verde, ≤ 2 % en el rojo y
≤ 0,75 % en el azul (0,1–1000 W/cm²); sin compañeros ascendentes, ≤ 0,24 % (`noche/iter1_t0_check.py`).

ΔE efectivo (cm⁻¹) de la pendiente de ln FIR, FIR_pure / FIR_band:

| P (W/cm²) | 160–300 K | 303–573 K | 353–453 K | 300–700 K |
|---|---|---|---|---|
| 0,1 | 630 / 631 | 649 / 649 | 648 / 649 | 649 / 649 |
| 0,8 | 629 / 637 | 648 / 654 | 648 / 654 | 649 / 653 |
| 2 | 627 / 646 | 648 / 660 | 648 / 661 | 648 / 658 |
| 10 | 619 / 682 | 648 / 685 | 648 / 690 | 648 / 681 |
| 100 | 591 / 764 | 646 / 776 | 645 / 792 | 646 / 759 |
| 290 | 566 / 794 | 644 / 849 | 643 / 876 | 644 / 822 |

Prefactor de FIR_pure (300–700 K): 5,45–5,50 (es el dato de entrada de Suta).

Intensidades (T del máximo en K: verde / rojo / azul; verde(10 K)/verde(300 K)):
0,1 W/cm²: 320 / 500 / 10; 0,72. 0,8: 340 / 10 / 10; 0,73. 10: 420 / 10 / 10; 0,69.
290: 540 / 180 / 10; 0,56.

Rojo/verde a 100 / 300 / 500 K: 0,1 W/cm²: 0,050 / 0,045 / 0,214; 0,8: 0,312 / 0,175 / 0,238;
10: 1,40 / 0,71 / 0,37; 290: 2,65 / 1,37 / 0,55.

Pendientes log-log a 0,8 W/cm² (541 / 520 / rojo): 300 K 1,75 / 1,75 / 2,46; 480 K 1,89 / 1,89 / 2,05.

τ(4S3/2) por excitación directa: 115 µs de 40 a 220 K, 112 µs a 300 K, 108 µs a 360 K, 87 µs a 440 K.
τ(4F9/2) = 490 µs a toda T. El τ verde corto no es un error: k_CR6·densidad·n(Er1) =
2,79e-17 × 1,38e22 × 0,019 ≈ 7300 s⁻¹ domina sobre k_R6 = 1510 s⁻¹ (en el pre-registro supuse que
CR6 era despreciable; me equivoqué en la cuenta).

### Evaluación contra el mapa

| fila | resultado del modelo | estado |
|---|---|---|
| O1 | FIR crece con T; meseta por debajo de ~170 K (T_on, sube con P: resultado conocido, no cuenta) | reproduce |
| O2 | FIR_pure: 643–649 cm⁻¹ en ventanas ≥ 300 K → **no** da 713–817. FIR_band llega a 760–876 solo a 100–290 W/cm²; a las potencias publicadas (0,8–2 W/cm²) da 654–661 | contradice (lo de FIR_band estaba pre-registrado y no alcanza a las potencias reales) |
| O3 | prefactor 5,5 contra 8,06 (Zhou) y 9,13 (Tong) | contradice (es la entrada de Suta) |
| O4 | verde **baja** al enfriar por debajo de 320–540 K (0,56–0,73 a 10 K); azul sube ×4–5; rojo con meseta o máximo muy chato (ver abajo). Langping ve 2,2× a 160 K respecto de 300 K a 0,1 W/cm²; el modelo da 0,85× | contradice; además el signo del verde al enfriar contradice mi pre-registro |
| O5 | verde baja sobre 300 K a P ≤ 2 W/cm²; a ≥ 100 W/cm² sube hasta ~500 K | reproduce a baja potencia |
| O6 | rojo/verde sube con T sobre 300 K a P ≤ 0,8 W/cm² (no pre-registrado: dije que bajaba) y baja con T a P ≥ 10 W/cm². Xu (0,8 W/cm²) ve el rojo pasar a ser el pico más alto a 483 K; el modelo da R/G 0,175 → 0,238, lejos de ese orden. Yu (290 W/cm²) ve R/G creciente: el modelo da el signo opuesto | signo correcto a baja P, magnitud no; contradice a alta P |
| O7 | pendientes del verde suben con T (+0,14 entre 300 y 480 K contra +0,29/+0,51 en Xu); la del rojo baja (2,46 → 2,05) contra +0,39 en Xu; el modelo da pendientes iguales para 520 y 541 nm y Xu no (1,53 contra 1,22 a 303 K, señal de calentamiento por el láser, O14) | signo del verde sí; rojo contradice |
| O9 | τ verde casi constante hasta 300 K (115 µs) contra 430 → 166 µs (Langping) y 410 µs a 77 K (Suta, no contraste); 353 → 453 K: −21 % contra −38 % (Xu) | contradice |
| O10 | τ rojo constante (490 µs; coincide con 497 µs de Xu a 230 K) sin el aumento a 583–605 µs a 55 K | contradice la tendencia |
| O12 | azul baja con T | reproduce |

**Máximo del rojo a alta potencia (no pre-registrado).** A 290 W/cm² (la potencia de Yu 2014) el rojo
tiene un máximo a 180 K. Prueba de robustez (`noche/iter1_redmax.py`): con ħω = 300 cm⁻¹ el máximo
está en 200–220 K, con 359 cm⁻¹ en 120–180 K y con 450 cm⁻¹ desaparece; la caída por debajo del
máximo es de 0–11 %. En la misma figura de Yu el verde y el azul también tienen máximo, y el modelo
no los da (verde máximo a 500 K, azul monótono). No cumple el criterio 3 (no se sostiene en el rango
de ħω) y reproducirlo solo en un canal sería elegir a mano: **no es emergente**.

**Mecanismo del verde que baja al enfriar.** A temperatura ambiente parte de la población que el ETU
6→9 saca del par verde vuelve por 9→8→7→6 (k_NR8 = 43450 s⁻¹ frente a k_R8 = 2330 s⁻¹). Al enfriar,
k_NR8 cae ~9× (gap de 4000 cm⁻¹, p = 11) y 2H9/2 pasa a emitir en azul y en 2H9/2→4I13/2: el
reciclado se corta y el verde pierde. Es una consecuencia directa del k_NR8 de Anderson con la ley
del gap; no reproduce ningún dato.

Estado de la iteración 1: completada. El modelo reproduce O1, O5 (baja P) y O12, y contradice O2,
O3, O4, O6 (alta P), O7 (rojo), O9 y O10.

## Iteración 2: relajación cruzada 4S3/2 asistida por fonones

### Pregunta
O9: la vida media de 4S3/2 cae con T en varios trabajos (Langping 2023: 430 → 166 µs entre 40 y
300 K, core-shell 20Yb2Er, excitación directa; Xu 2024: 226 → 141 µs entre 353 y 453 K, micro
20Yb2Er; Yu 2014 Fig. 7d: los tiempos de decaimiento bajan con T). No se usó en la calibración. El
modelo de la iteración 1 da 115 µs casi constantes hasta 300 K, porque el canal dominante de pérdida
del par verde es la relajación cruzada CR6 de Anderson (≈ 7300 s⁻¹ contra k_R6 = 1510 s⁻¹) y es
independiente de T.

### Fuentes
- Langping 2023 (texto completo): "for Er3+-rich core–shell systems at the relatively high temperatures,
  the lattice vibration (i.e., phonons) facilitates to fill the energy gaps in various CR processes"
  (p. 4) y "the cascade phonon-assisted Er3+-Er3+ CR plays a dominant role in the upconversion energy
  loss at room temperature. A cryogenic environment can suppress the harmful CR" (conclusión).
  Para su muestra de bajo dopaje atribuyen la caída de τ a relajaciones no radiativas (p. 4); el
  mecanismo lo tomo de su discusión de alto dopaje.
- Suta 2025: "both at higher Er 3+ and Yb 3+ contents in β-NaYF4, the 4S3/2 level decays faster due to
  additional cross-relaxation based on an interaction with neighbouring Er 3+ or Yb 3+ ions" (p. 7096).
- Ley de la transferencia asistida por fonones (emisión de p fonones, factor (1+n)^p): Miyakawa y
  Dexter 1970 (PRB 1, 2961) y Auzel 2004 (Chem. Rev. 104, 139): ambos **inaccesibles** (Unpaywall
  is_oa = False). La forma es la misma que la de la relajación multifonónica de Riseberg y Moos 1968
  (texto completo), que es lo que se usa.

### Propuestas
1. **CR6 asistida por fonones (elegida).** La CR6 de Anderson (4S3/2 + 4I15/2 → 4I11/2 + 4I13/2)
   tiene un desajuste de (18300 − 10200) − 6500 = 1600 cm⁻¹ desde 4S3/2 y 2250 cm⁻¹ desde 2H11/2 que
   se emite como fonones: k_CR6(T) = k_CR6(T0)·[(1+n(T))/(1+n(T0))]^(ΔE/ħω). Anclada al valor de
   Anderson a T0, con el mismo ħω de las multifonónicas. **Parámetros libres nuevos: 0.** Bloque
   `PhononAssistedEnergyTransfer` ya existente en `transitions.py`. Interruptor: `cr_thermal`
   (1 encendido, 0 = CR6 constante como en la iteración 1).
2. Retrotransferencia 4S3/2 → Yb (4S3/2 + ²F7/2 → 4I13/2 + ²F5/2, desajuste 1600 cm⁻¹; Suta y
   Langping la mencionan). Descartada: canal que Anderson no tiene, su tasa a T0 sería un parámetro
   libre sin anclaje.
3. Todas las transferencias (ETU Yb→Er, retrotransferencias, CR4, UC2) asistidas por fonones.
   Descartada por ahora: varios desajustes con baricentros son menores que ħω (3→7: −100 cm⁻¹,
   7→3: +100 cm⁻¹) y en esos casos la ley necesita una fracción resonante libre por transferencia.

### Pre-registro (antes de implementar y simular)
Cuentas a mano con ħω = 359 cm⁻¹ (n(300 K) = 0,218):
- τ verde: a T → 0 el factor es (1/1,218)^4,46 = 0,415, así que k_CR6 ≈ 3030 s⁻¹ y
  τ(40 K) ≈ 1/(1375 + 3030) ≈ 227 µs; τ(300 K) ≈ 112 µs; cociente ≈ 2,0 (Langping: 2,6). Con
  ħω = 300 cm⁻¹ el factor a T → 0 es 0,24 y el cociente ≈ 2,7. Entre 353 y 453 K espero τ ≈ 88 → 54 µs
  (cociente ≈ 1,6; Xu: 1,60). Los valores absolutos quedan por debajo de los publicados (Anderson
  tiene una CR mucho más fuerte que esas muestras). Esto es consecuencia directa de una sola tasa:
  aunque coincida con los datos, no puede contar como emergente (criterio 4).
- Intensidad verde contra T: la CR más débil al enfriar compensa en parte la pérdida por el corte
  del reciclado 9→8→7→6 de la iteración 1. Espero verde(10 K)/verde(300 K) entre 1,0 y 1,5 a baja
  potencia y que el verde caiga más rápido que en la iteración 1 por encima de 300 K. Es posible
  que aparezca un máximo del verde a temperatura intermedia, pero no sé estimar dónde; si aparece
  entre 100 y 200 K y con un cociente máx/300 K de 1,5–2,5 a 0,1 W/cm² (Langping 2,2× a ~160 K) lo
  consideraré una magnitud no pre-registrada.
- Rojo/verde: al subir T la CR más fuerte quita verde, así que R/G debería crecer con T más que en
  la iteración 1 (posible cambio de signo a alta P; no lo sé).
- Azul: sigue bajando con T. FIR (O1–O3) y τ rojo (O10): sin cambios.
- Pendientes log-log: sin cambio apreciable (la CR6 es lineal en la población del verde).
Sorprendente: un cambio en FIR o en τ rojo; un máximo del verde muy marcado (> 3×).

### Resultados
Implementación: commit 9a761e7 (después del pre-registro 82c3998). Test de reducción
`tests/test_iter2_cr.py`: con `cr_thermal` = 0 se recupera el barrido de la iteración 1 a rtol 1e-6
(5 temperaturas × 2 potencias) y los desajustes salen de los niveles (1600 y 2250 cm⁻¹). 32 tests pasan.
Barridos: `noche/iter2_run.py` (ħω = 359, 300 y 450 cm⁻¹, con y sin CR asistida; ~5,4 min cada uno).
Tabla completa: `noche/data/compare_iter2.txt` (`noche/compare.py`). Figuras: `noche/figs/iter2/`.

| métrica | iter. 1 (359) | iter. 2 (359) | iter. 2 (300) | iter. 2 (450) | literatura |
|---|---|---|---|---|---|
| τ verde 40 K / 300 K | 1,02 | 2,01 | 2,87 | 1,46 | 2,59 (Langping) |
| τ verde 353 K / 453 K | 1,32 | 1,98 | 2,47 | 1,61 | 1,60 (Xu) |
| τ verde 40 K (µs) | 115 | 226 | 322 | 164 | 430 (Langping) |
| visible 160 K / 300 K a 0,1 W/cm² | 0,79 | 1,35 | 1,47 | 1,18 | 2,2 (Langping, máximo ~160 K) |
| T del máximo del verde a 0,1 W/cm² (K) | 320 | 10 | 140 | 10 | ~160 (Langping) |
| T del máximo del verde a 290 W/cm² (K) | 540 | 380 | 360 | 440 | ~100–150 (Yu) |
| verde 10 K / 300 K a 290 W/cm² | 0,57 | 0,61 | 0,42 | 0,80 | > 1 (Yu: máximo ~100 K) |
| R/G 483 K / 303 K a 0,8 W/cm² | 1,21 | 1,21 | 2,32 | 0,88 | > 1 (Xu) |
| R/G 400 K / 10 K a 290 W/cm² | 0,32 | 0,32 | 0,20 | 0,51 | > 1 (Yu) |

Contra el pre-registro: el cociente de τ (2,0) y la subida del verde al enfriar a baja potencia (1,3–1,4)
salieron como estaban escritos; τ entre 353 y 453 K cae más de lo previsto (1,98 contra 1,6). R/G,
pendientes, FIR y τ rojo no cambian, como se esperaba.

**Máximo del verde a ~140 K con ħω = 300 cm⁻¹.** Con ħω = 300 cm⁻¹ (el modo Raman de 307 cm⁻¹ de
Dubey) aparece un máximo del verde a 140 K (0,1 W/cm²) y 160 K (0,8 W/cm²), con máx/300 K = 1,45, y del
visible a 120 K: cerca del máximo de Langping (~160 K, 2,2×). Con 359 y 450 cm⁻¹ no hay máximo (verde
monótono hasta 10 K). Sale de la competencia entre la CR6 que se apaga al enfriar y el reciclado
9→8→7→6 que también se corta (iteración 1), pero depende de ħω dentro del rango de literatura, así
que no cumple el criterio 3, y a 290 W/cm² (Yu) sigue sin haber máximo a baja T. **No es emergente.**

Estado: el mecanismo reproduce el orden de magnitud de la caída de τ(4S3/2) (O9: "reproduce lo
conocido", consecuencia directa de una tasa) y mejora el signo de O4 a baja potencia (sube al enfriar),
sin máximo robusto. O6 (R/G a alta potencia) y O4 a 290 W/cm² siguen en contra.

## Iteración 3: transferencias Yb ↔ Er asistidas por fonones

### Pregunta
O6 a alta potencia y O4 a 290 W/cm²: Yu 2014 (bulk 18Yb2Er, 290 W/cm²) ve R/G creciente con T y un
máximo de todas las bandas cerca de 100 K; Xu 2024 (0,8 W/cm²) ve el rojo superar al verde a 483 K.
El modelo (iteraciones 1 y 2) da R/G **decreciente** con T a ≥ 10 W/cm², el azul crece ×4–5 al
enfriar y el verde a 290 W/cm² baja al enfriar por debajo de 360–440 K. Dos grupos, no usado en la
calibración.

### Fuentes
- Li 2014 (texto completo): "The energy transfer between donors (Yb3+) and acceptors (Er3+) requires
  low-energy phonons to match these energy diﬀerences." y desajustes de 40–90 cm⁻¹ entre
  ²F7/2→²F5/2 (~10260 cm⁻¹) y 4I15/2→4I11/2 (~10300) o 4I11/2→4F7/2 (~10350) (p. 3, Fig. 1b).
- Langping 2023: "the hindered Yb3+→Er3+ energy transfer (caused by retrogressive phonon activity)"
  como causa de la baja de intensidad a baja T en su muestra de bajo dopaje (p. 3).
- Anderson 2013: el ETU 6→9 lleva a Er "into the 4G,2K manifold, above 2H9/2" (p. 37): el estado
  final del modelo (Er9, 26100 cm⁻¹) es una representación del manifold.
- Ley: la misma de la iteración 2 (Miyakawa–Dexter, inaccesible; forma de Riseberg–Moos).

### Propuestas
1. **ETU y retrotransferencias con |ΔE| ≥ ħω asistidas por fonones (elegida).** Con los baricentros
   del modelo: 5→8: ΔE = 10200 − (24500 − 15000) = +700 cm⁻¹; 6s→9: +2400; 6h→9: +3050; 9→5
   (retrotransferencia): (26100 − 15000) − 10200 = +900. Todas exotérmicas (emisión de fonones),
   k(T) = k(T0)·[(1+n(T))/(1+n(T0))]^(ΔE/ħω), ancladas a Anderson. **Parámetros libres: 0.**
   Interruptor `et_thermal`. Las transferencias casi resonantes (1→3, 3→1: 0; 3→7: −100; 7→3: +100)
   quedan constantes.
2. Las transferencias casi resonantes con los desajustes de Li 2014 (−40/−90 cm⁻¹, absorción de un
   fonón). Descartada: necesita una fracción resonante libre por transferencia para que el ETU no se
   anule a T → 0, y es la explicación conocida de los autores.
3. Población térmica de la componente |1⟩ de ²F5/2 del Yb (Suyver 2005, Δ01 = 39 cm⁻¹). Descartada:
   el número solo está en un abstract (`abstract-only`), no se puede usar como parámetro.

### Pre-registro (antes de implementar y simular)
Factores a T → 0 con ħω = 359 cm⁻¹: 5→8 ×0,58; 6s→9 ×0,27; 6h→9 ×0,19; 9→5 ×0,61. A 500 K: 5→8
×1,6; 6s→9 ×5,0; 9→5 ×1,8.
- Azul: su fuente principal (6→9 →8) cae ×0,27 al enfriar, contra el aumento ×4–5 por k_NR8; espero
  que el azul a baja T quede cerca del valor de 300 K, quizás con un máximo suave; no sé el signo neto.
- Rojo (de 9→5): baja al enfriar. R/G: espero que pase a **crecer** con T también a alta potencia
  (signo de Yu y Xu), magnitud incierta.
- Verde a 290 W/cm²: al enfriar se drena menos hacia 9, así que espero que suba respecto de la
  iteración 2 a baja T; no espero un máximo cerca de 100 K simultáneo en verde, rojo y azul (eso
  sería sorprendente).
- Sobre 300 K: el drenaje 6→9 crece (×5 a 500 K) y el verde cae más rápido; el rojo cae menos.
- Pendientes log-log: el verde depende más de la potencia a alta T (más drenaje no lineal), así que
  su pendiente debería **bajar** algo a alta T respecto de la iteración 2. La del rojo: no sé.
- τ por excitación directa, FIR en equilibrio: sin cambios (no hay Yb excitado).
Sorprendente: máximo simultáneo de las tres bandas a 80–200 K a 290 W/cm²; cambio de FIR por encima
de 300 K.

### Resultados
Implementación: commit d9fbb72 (después del pre-registro 3d1ef4f). `tests/test_iter3_et.py`: con
`et_thermal` = 0 se recupera el barrido de la iteración 2 a rtol 1e-6; los desajustes salen de los
niveles (700, 2400, 3050, 900 cm⁻¹). 46 tests pasan. Barridos `noche/iter3_run.py` (~9 min cada uno);
tabla `noche/data/compare_iter3.txt`; figuras `noche/figs/iter3/`. La ablación son los barridos de
la iteración 2 (mismo modelo con `et_thermal` = 0).

| métrica | iter. 2 (359) | iter. 3 (300) | iter. 3 (359) | iter. 3 (450) | literatura |
|---|---|---|---|---|---|
| R/G 400 K / 10 K a 290 W/cm² | 0,32 | 5,27 | 2,83 | 1,72 | > 1, monótono (Yu, bulk) |
| R/G 483 K / 303 K a 0,8 W/cm² | 1,21 | 4,35 | 2,96 | 2,07 | > 1 (Xu) |
| verde 10 K / 300 K a 290 W/cm² | 0,61 | 2,54 | 1,79 | 1,32 | máximo ~100 K (Yu) |
| T del máximo del rojo a 290 W/cm² (K) | 160 | 260 | 240 | 120 | ~100 bulk / ~150 nano (Yu) |
| azul 10 K / 300 K a 290 W/cm² | 4,94 | 5,87 | 4,19 | 2,30 | máximo ~100 K (Yu) |
| T del máximo del rojo a 0,8 W/cm² (K) | 10 | 340 | 360 | 420 | — |
| pendiente del rojo a 483 K, 0,8 W/cm² | 2,07 | 2,43 | 2,45 | 2,45 | 1,69 (Xu) |

Contra el pre-registro:
- R/G pasa a crecer con T a todas las potencias y para los tres ħω, como estaba pre-registrado
  (signo de Yu y de Xu). **Reproduce la tendencia; no es emergente** (estaba en el pre-registro).
- El verde a 290 W/cm² sube al enfriar (pre-registrado) y es monótono: no hay máximo cerca de 100 K
  en ninguna banda a la vez. El azul sigue creciendo al enfriar (×2,3–5,9).
- **No pre-registrado:** a baja potencia (0,1–0,8 W/cm²) el rojo **sube** con T hasta un máximo a
  340–460 K (el pre-registro decía "el rojo cae menos"). Busqué datos de la intensidad roja absoluta
  sobre 300 K a baja potencia: Dubey 2023 (29 nm, 1,2 W) dice que la emisión a 660 nm baja
  ("the overall UC intensities decrease signi cantly", Fig. 7a), Li 2014 ve caída gradual en
  partículas > 32 nm (sin dar la potencia en el texto), Xu 2024 solo da el orden de los picos. El único
  dato comparable lo contradice: no cuenta.
- La pendiente del rojo a alta T empeora (2,45 contra 1,69 de Xu).

Estado: O6 reproducido en signo (esperado). O4 sigue en contra (sin máximo a 100–160 K). O7 (rojo)
en contra.

## Iteración 4: transferencias Yb ↔ Er casi resonantes que absorben un fonón

### Pregunta
O4: Yu 2014 (290 W/cm²) y Langping 2023 (0,1 W/cm²) ven un máximo de la intensidad UC a 100–160 K
(dos grupos con texto completo, un tercero, Suyver 2005, por cita). Después de las iteraciones 2 y 3
el modelo sube monótonamente al enfriar en verde y azul: le falta algo que se congele a baja T.

### Fuentes
- Li 2014 (texto completo, p. 3, Fig. 1b): "there are slight energy diﬀerences (40 −90 cm −1) between
  the 2F7/2 → 2F5/2 transition ( ∼10260 cm −1) of Yb 3+ ions and the 4I15/2 → 4I11/2 (∼10 300 cm −1) or
  4I11/2 → 4F7/2 (∼10 350 cm−1) transition of Er 3+ ions." y "The energy transfer between donors
  (Yb3+) and acceptors (Er3+) requires low-energy phonons to match these energy diﬀerences."
  → 1→3: ΔE = 10260 − 10300 = −40 cm⁻¹; 3→7: 10260 − 10350 = −90 cm⁻¹ (absorción de un fonón de
  esa energía); las inversas 3→1 (+40) y 7→3 (+90) emiten un fonón.
- Langping 2023: "the hindered Yb3+→Er3+ energy transfer (caused by retrogressive phonon activity)".
- Yu 2014 atribuye el máximo a la componente |1⟩ del ²F5/2 (Suyver): misma física de activación de
  la transferencia Yb→Er, con un número (39 cm⁻¹) que solo está en un abstract.

### Propuestas
1. **1↔3 y 3↔7 asistidas por un fonón con los desajustes de Li 2014 (elegida).** Absorción:
   k(T) = k(T0)·[r + (1 − r)·n(ε,T)/n(ε,T0)], ε = |ΔE|; emisión: k(T) = k(T0)·[1 + n(ε,T)]/[1 + n(ε,T0)].
   Bloque `PhononAssistedEnergyTransfer` (rama de un fonón) con el desajuste dado. **Un parámetro
   libre:** la fracción resonante r (compartida por 1→3 y 3→7), que representa transferencias
   resonantes entre componentes Stark o por ensanchamiento inhomogéneo; sin valor de literatura.
   Interruptor `nr_thermal` (0 = constantes, iteración 3).
2. Componente |1⟩ del ²F5/2 del Yb (Suyver 2005). Descartada: Δ01 solo en un abstract.
3. Sin fracción resonante (r = 0, cero parámetros libres). Descartada como modelo, pero se corre
   como caso límite: predice que la UC se anula a T → 0, y los datos (Langping: "varies little" de
   300 a 40 K) lo contradicen de entrada.

No ajusto r: barro r = 0,1, 0,3 y 0,6 con ħω = 359 cm⁻¹ y reporto los tres.

### Pre-registro (antes de implementar y simular)
- Por debajo de ~60 K (ε = 40 cm⁻¹) y ~130 K (90 cm⁻¹) los ETU 1→3 y 3→7 tienden a r·k(T0); las
  retrotransferencias 3→1 y 7→3 bajan a 0,17 y 0,35 de su valor a 300 K.
- A baja potencia el verde, el rojo y el azul pasan por un **máximo** a T intermedia porque el
  factor r + (1−r)·n/n0 de los dos ETU se multiplica contra la ganancia ×1,2–1,6 de las iteraciones
  2–3. Espero que el máximo se corra a T más bajas al subir r: con r = 0,1 cerca de 250–300 K, con
  r = 0,6 por debajo de 150 K. A 290 W/cm² (saturado) el efecto es más chico y el máximo, si aparece,
  a T más baja que a baja potencia.
- Sobre 300 K los ETU 1→3 y 3→7 crecen casi lineal con T (×~1,6 a 483 K en la parte no resonante):
  la caída térmica del verde se frena y las pendientes log-log bajan a alta T (más saturación), al
  revés que en Xu.
- R/G sigue creciendo con T. τ por excitación directa y FIR en equilibrio: sin cambios.
- Si aparece un máximo simultáneo de las tres bandas a 100–160 K, es la explicación que dan los
  propios autores (transferencia Yb→Er impedida a baja T): **reproduce lo conocido**, no emergente.
Sorprendente: que la posición del máximo no dependa de r; que R/G se vuelva no monótono; que las
pendientes suban a alta T.

### Adenda al pre-registro de la iteración 4: dinámica tras un pulso (O11), antes de calcularla
Agrego al barrido el observable dinámico O11 (Yu 2014 Fig. 7: los tiempos de subida suben de 10 a
30 K y bajan por encima de 30 K, en todas las transiciones y tamaños β). Código: `noche/dynamics.py`
(pulso instantáneo que excita 1 % de los Yb; τ_D de la pendiente tardía, τ_R del tiempo del máximo
con la forma de Vial). Lo probé solo con `thermal` = 0 a 300 K para depurarlo (verde: τ_R = 21 µs,
τ_D = 423 µs). Expectativa: τ_R está fijado sobre todo por el vaciado del Yb (613 s⁻¹ más la
transferencia 1→3); al enfriar, 1→3 tiende a r·k, así que τ_R debería **crecer** al enfriar, de forma
monótona; entre 10 y 30 K el factor de 1→3 cambia poco (n(40 cm⁻¹) pasa de 0,003 a 0,17 contra 4,7
a 300 K), así que espero τ_R casi plano en ese tramo, sin el máximo a 30 K de Yu. Sorprendente: un
τ_R no monótono con máximo a 20–60 K.

### Resultados
Implementación: commit 992b9ba (después del pre-registro 8f38501; adenda O11 en b21772d antes de
calcular la dinámica). `tests/test_iter4_nr.py`: con `nr_thermal` = 0 se recupera la iteración 3 a
rtol 1e-6; desajustes de Li 2014 verificados. 60 tests pasan. Barridos `noche/iter4_run.py` (r = 0,
0,1, 0,3, 0,6), tabla `noche/data/compare_iter4.txt`, figuras `noche/figs/iter4_r0.3/`. Dinámica:
`noche/iter4_dynamics.py` → `noche/data/iter4_dynamics.json` (`noche/dyn_table.py` la imprime).

| métrica | iter. 3 | r = 0 | r = 0,1 | r = 0,3 | r = 0,6 | literatura |
|---|---|---|---|---|---|---|
| T máx. visible a 0,1 W/cm² (K) | 10 | 300 | 280 | 200 | 10 | ~160 (Langping) |
| visible 160 K / 300 K a 0,1 W/cm² | 1,32 | 0,62 | 0,76 | 1,07 | 1,60 | 2,2 (Langping) |
| visible 40 K / 300 K a 0,1 W/cm² | 1,37 | 0,015 | 0,19 | 0,92 | 2,59 | "varies little" 300→40 K (Langping) |
| T máx. verde / rojo / azul a 290 W/cm² (K) | 10 / 240 / 10 | 200 / 340 / 100 | 160 / 320 / 70 | 10 / 260 / 60 | 10 / 10 / 80 | ~100 en las tres (Yu bulk) |

- Con r = 0 la UC se apaga al enfriar (×0,015 a 40 K): descartado por los datos, como se anticipó.
- Con r = 0,1–0,3 aparece un máximo del visible a baja potencia, pero a 200–280 K y de solo
  ×1,0–1,1, contra ~160 K y ×2,2 en Langping. A 290 W/cm² ningún r da el máximo de las tres bandas
  cerca de 100 K. La posición del máximo depende de r (lo pre-registrado).
- **Defecto de formulación:** la fracción resonante r solo entra en las transferencias que absorben
  (1→3, 3→7); las inversas (3→1, 7→3) caen a 0,17 y 0,35 a T → 0 sin fracción resonante. Para r
  grande eso rompe el balance entre ida y vuelta y produce el ×2,6–3 del verde a baja T con r = 0,6:
  es un artefacto, no física.
- O11 (dinámica, adenda): τ_R del verde **crece** al enfriar en forma monótona para todo r y para
  la ablación (r = 0,3: 21 µs a 300 K → 38 µs a 10 K; casi plano entre 10 y 30 K), como estaba
  pre-registrado: no reproduce el máximo a 30 K de Yu.
- O7b (nueva fila, Sagaidachnaya y Kochubey 2020, β-NaYF4:17Yb,3Er de 440 nm, 1,5–9,4 W/cm²,
  22–55 °C): la pendiente del verde sube hasta 6 % (baja intensidad) y 16 % (alta) entre 30 y 50 °C;
  el modelo da entre −1 % y +2 % entre 300 y 320 K en todas las iteraciones
  (`noche/slopes_saratov.py`). Contradice.

Estado de la iteración 4: **mal condicionada** (el resultado depende del parámetro libre r y la
formulación de la fracción resonante es asimétrica). Se deja en el código con `nr_thermal` y se
**apaga por defecto** en el modelo final (`nr_thermal` = 0).

## Iteración 5: O10 y O9 con el observable correcto (decaimiento de la UC tras pulso de 980 nm)

### Pregunta
Al calcular τ_D para la adenda de la iteración 4 noté que **comparé mal O10 y la parte de Xu de O9
en la iteración 1**: Xu 2024 midió los decaimientos con láser de 980 nm ("The spectrum and lifetime of
the system were recorded using a transient fluorescence spectrometer (FLS-980) and 980 nm laser",
sec. 2.1) y ajustó una exponencial y = y0 + A·e^(−t/τ) (eq. 1, Fig. 5): es el decaimiento de la
**UC**, que incluye la alimentación desde niveles intermedios, no la vida media del nivel emisor
que usé (excitación directa, 490 µs constante para 4F9/2). No hay mecanismo nuevo: es una
re-evaluación del modelo actual (iteración 3; `nr_thermal` = 0 por defecto desde 21fbefb+).

Lo que ya vi antes de este pre-registro (`noche/data/iter4_dynamics.json`, ħω = 359, τ_D por la
pendiente logarítmica entre 30 % y 5 % del máximo): con la configuración de la iteración 3 el τ_D del
rojo es ~511 µs a 230 K y 595 µs a 55 K con 1 % de Yb excitado (~460 y 553 µs con 10 %); Xu: 497 →
583 µs (653 nm) y 496 → 605 µs (661 nm). Lo declaro: el número a 359 cm⁻¹ no es una predicción ciega.

### Pre-registro del protocolo de evaluación (antes de correrlo)
Script `noche/iter5_decays.py` (se escribe después de este commit). Configuraciones:
iteración 1 (`cr_thermal` = `et_thermal` = 0), iteración 2 (`et_thermal` = 0), iteración 3 (modelo
final), ablación total (`thermal` = 0) y ablaciones por tasa (la T de un solo bloque multifonónico
fijada en T0: norad3, norad4, norad8, norad9). ħω ∈ {300, 359, 450} cm⁻¹; fracción de Yb excitada
x0 ∈ {0,003, 0,01, 0,1}; dos extracciones de τ: pendiente logarítmica 30 %→5 % y ajuste
y0 + A·e^(−t/τ) desde el máximo hasta 5 % (la forma de Xu). Temperaturas: 55 y 230 K (rojo, 653 =
rad51) y 353 y 453 K (verde 541 = rad6s1 y 520 = rad6h1).

Criterio fijado ahora para declarar emergente el τ_D(T) del rojo:
(i) en el modelo final, τ_D(55 K)/τ_D(230 K) del rojo entre 1,10 y 1,30 (Xu: 1,17 y 1,22) y τ_D(230 K)
dentro de ±20 % de 497 µs, para **todo** ħω ∈ {300, 359, 450} y x0 ∈ {0,003, 0,01, 0,1} y con las dos
extracciones; (ii) con `thermal` = 0 el cociente cae a ≤ 1,03; (iii) que no lo fije una sola tasa
puesta a mano (la ablación por tasa dice qué bloques importan; si basta una sola tasa multifonónica
y su efecto es directo sobre el nivel emisor, no cuenta; si viene de niveles alimentadores, sí).
Expectativa honesta: el cociente se mantiene en 1,1–1,3 para ħω 300–450 y x0 ≤ 0,01; con x0 = 0,1
el valor absoluto baja (~460 µs) y puede salir del ±20 %. El verde 541 entre 353 y 453 K: espero un
cociente ~1,3–1,6 (Xu: 1,60) con valores absolutos ~1,5× mayores que los de Xu.

### Resultados
`noche/iter5_decays.py` → `noche/data/iter5_decays.json` (ajuste de Xu y0 + A·e^(−t/τ)):

| configuración | ħω | x0 | rojo 55 / 230 K (µs) | cociente | 541 nm 353 / 453 K (µs) | cociente |
|---|---|---|---|---|---|---|
| final (iter. 3) | 300 | 0,003 / 0,01 / 0,1 | 1170/893 · 727/632 · 202/323 | 1,31 · 1,15 · **0,63** | 335/200 · 188/96 · 106/40 | 1,67 · 1,96 · 2,63 |
| final (iter. 3) | 359 | 0,003 / 0,01 / 0,1 | 1000/827 · 669/573 · 189/446 | 1,21 · 1,17 · **0,42** | 365/250 · 213/126 · 126/60 | 1,46 · 1,69 · 2,11 |
| final (iter. 3) | 450 | 0,003 / 0,01 / 0,1 | 848/768 · 595/513 · 401/567 | 1,11 · 1,16 · **0,71** | 395/309 · 240/167 · 147/89 | 1,28 · 1,43 · 1,66 |
| thermal = 0 | todos | todos | constante | 1,00 | constante | 1,00 |
| Xu 2024 | — | — | 583/497 (653 nm), 605/496 (661 nm) | 1,17 / 1,22 | 226/141 | 1,60 |

Ablaciones por tasa (ħω = 359, x0 = 0,01): fijar en T0 la T de norad3 o norad4 no cambia nada; la de
norad8 baja el cociente del verde de 1,69 a 1,15 y la de norad9 sube el del rojo a 1,33. La
dependencia térmica del decaimiento de la UC viene de los niveles que alimentan a los emisores
(2H9/2 → 4F7/2 → verde; 4G11/2 → 4F9/2 y 2H9/2), no del nivel emisor (Xu la atribuye a Γrad + k_nr
del propio nivel emisor, su eq. 2).

Evaluación contra el criterio pre-registrado:
- (i) **falla**: con x0 = 0,1 el cociente del rojo se invierte (0,42–0,71) y con x0 = 0,003 el valor
  absoluto a 230 K (768–893 µs) sale del ±20 % de 497 µs. Solo con x0 ≈ 0,01 el modelo da 573 → 669 µs
  (×1,17) contra 497 → 583 µs (×1,17) de Xu, y 213 → 126 µs contra 226 → 141 µs en el verde 541.
- (ii) se cumple (thermal = 0 → ×1,00). (iii) el efecto viene de niveles alimentadores.
- Xu no informa la energía del pulso, así que no hay forma de fijar x0 desde la fuente. Exploración
  **post hoc** (no cuenta, `noche/iter5_cwoff.py`): apagando una excitación continua de 0,8 W/cm²
  (régimen débil) los decaimientos son de ms (rojo ~4,5 ms, verde ~1,3 ms), un orden de magnitud más
  largos que los de Xu.

**No es emergente**: la coincidencia cuantitativa con Xu existe en una franja estrecha de intensidad de
pulso que la fuente no permite verificar. Queda como la dirección abierta más prometedora: un
experimento de decaimiento de la UC contra T a energía de pulso conocida (o los datos crudos de Xu)
decidiría si el mecanismo "alimentación térmica desde 2H9/2 y 4G11/2" es el correcto.

Corrección del mapa: la fila O10 (y la parte de Xu de O9) debe compararse con el decaimiento de la UC,
no con la vida media por excitación directa; el veredicto de la iteración 1 ("contradice") se reemplaza
por "depende de la energía del pulso, no decidible con la fuente".

## Chequeo posterior de la iteración 3 (robustez del signo de R/G)
`noche/iter3_mismatch_check.py` (post hoc). El estado final del ETU 6→9 (Er9, 26100 cm⁻¹) representa
el manifold 4G,2K; con un estado final más alto el desajuste baja. Con desajustes de 1000/1650 cm⁻¹ en
lugar de 2400/3050: a 0,8 W/cm² R/G sigue creciendo con T (400 K/10 K = 1,55; 483 K/303 K = 1,78), pero
a 290 W/cm² **vuelve a bajar** (0,83 y 0,99). El signo correcto a alta potencia (Yu) depende del
estado final supuesto para el ETU desde el par verde; a baja potencia (Xu) es robusto.

## Estado final del mapa

Modelo final: `AndersonThermal` con `thermal` = `cr_thermal` = `et_thermal` = 1 y `nr_thermal` = 0
(iteraciones 1–3; la 4 queda apagada). Figuras clave en `noche/figs/final/` (`noche/figs_final.py`).

| fila | observable | estado con el modelo final |
|---|---|---|
| O1 | FIR crece con T, Boltzmann sobre ~170 K | reproduce (meseta a baja T = T_on conocido) |
| O2 | ΔE efectivo 713–817 cm⁻¹ (5 grupos) | **no reproduce**: 647–649 cm⁻¹ con FIR puro; la contaminación por 2H9/2→4I13/2 solo lo sube a 700–800 a 50–290 W/cm², no a las potencias publicadas (fig. `delta_e.png`) |
| O3 | prefactor 8–9 | no reproduce (5,5, entrada de Suta) |
| O4 | máximo de intensidad a 100–160 K | no reproduce (verde y azul crecen al enfriar; iteración 4 lo da solo con un parámetro libre y una formulación defectuosa) |
| O5 | apagado térmico sobre ambiente | reproduce a baja potencia |
| O6 | R/G crece con T | reproduce el signo (pre-registrado); a alta potencia depende del estado final del ETU 6→9 |
| O7 | pendientes log-log crecen con T (Xu) | verde: signo sí, magnitud no (+0,1 contra +0,3/+0,5); rojo: no. Xu da n(520) ≠ n(541) a la misma T, imposible con termalización rápida sin calentamiento por el láser |
| O7b | pendiente del verde +6–16 % entre 30 y 50 °C (Saratov) | no reproduce (−1 % a +2 %) |
| O8 | pendientes a 300 K, ~2 a baja densidad y bajando | reproduce (física de saturación esperada) |
| O9 | τ(4S3/2) baja con T | reproduce el orden del cociente (×1,5–2,9 contra ×2,6), valores absolutos ~2× más cortos; consecuencia directa de la CR6 asistida |
| O10 | decaimiento de la UC del rojo sube al enfriar (Xu) | coincide (573→669 µs contra 497→583 µs) solo con ~1 % de Yb excitado por el pulso; con 10 % se invierte. No decidible con la fuente |
| O11 | τ_R con máximo a 30 K (Yu) | no reproduce (monótono) |
| O12 | azul baja con T | reproduce |
| O13 | 2H9/2→4I13/2 en la banda verde a alta potencia | reproduce |
| O14 | calentamiento por láser | fuera del modelo |

## Lo más interesante

**No apareció una propiedad emergente** que cumpla los cuatro criterios. Los dos candidatos más cerca
fallaron en el criterio 3:
1. *Decaimiento de la UC contra T (iteración 5).* Con el modelo final y ~1 % de Yb excitado, el
   decaimiento del rojo pasa de 573 a 669 µs entre 230 y 55 K (Xu: 497 → 583 µs) y el del verde 541 nm de
   213 a 126 µs entre 353 y 453 K (Xu: 226 → 141 µs), sin parámetros libres y con un mecanismo distinto
   del que proponen los autores: la dependencia térmica viene de los niveles que alimentan a los
   emisores (2H9/2 → 4F7/2 para el verde, 4G11/2 para el rojo), no del nivel emisor. Pero con 10 % de
   Yb excitado el cociente del rojo se invierte, y Xu no informa la energía del pulso.
2. *Máximo del verde a ~140 K a 0,1 W/cm² (iteración 2, ħω = 300 cm⁻¹).* Sale de la competencia entre
   la CR6 asistida por fonones, que se apaga al enfriar, y el corte del reciclado 9→8→7→6; cerca del
   máximo de Langping (~160 K) pero solo para ħω = 300 cm⁻¹ y de ×1,45 contra ×2,2.

Lo que el modelo sí muestra con firmeza: (a) en el esquema de Anderson el par verde recupera, a
temperatura ambiente, parte de lo que el ETU 6→9 le saca (9→8→7→6), y ese reciclado se corta al
enfriar porque k_NR8 cae ~9×: el azul crece al enfriar y el verde pierde esa vía; (b) el signo de
R/G(T) depende de que los ETU desde el par verde y la retrotransferencia 9→5 estén asistidos por
fonones; (c) O2 (ΔE efectivo 713–817 cm⁻¹ en cinco grupos) no se explica con ecuaciones de tasa y
A constantes: la termalización rápida fija el FIR en Boltzmann con el gap de entrada.

Direcciones abiertas más prometedoras:
1. **Decaimiento de la UC contra T a energía de pulso conocida** (o los datos crudos de Xu 2024): el
   modelo predice que el cociente τ(55 K)/τ(230 K) del rojo pasa de ~1,2 a < 1 al subir la fracción de
   Yb excitada de 1 % a 10 %; es un test directo y barato del mecanismo de alimentación térmica.
2. **Estructura Stark de 2H11/2 y 4S3/2** para O2/O3: con las energías de los subniveles (necesitan
   una fuente accesible) el FIR integrado tendría un ΔE efectivo dependiente de T; es lo único dentro
   del marco de ecuaciones de tasa que puede subir la pendiente por encima de 650 cm⁻¹ a baja potencia.
3. **Calentamiento por el láser** como parte del modelo: Xu da n(520) − n(541) = 0,31 a 303 K, que con
   termalización rápida implica ~30 K de calentamiento por unidad de ln P; eso contamina O7 y O7b y
   debería modelarse antes de usar pendientes contra T como dato de contraste.

## Pasos salteados, fallas y limitaciones
- La sesión estuvo detenida ~8,5 h (02:50–11:30) por reposo de la máquina; el presupuesto se
  reinició a las 11:36.
- Inaccesibles (Unpaywall is_oa = False): Suyver 2005 y 2006, Renero-Lecuna 2011, Pollnau 2000,
  Miyakawa–Dexter 1970, Auzel 2004. Suyver 2005 solo como abstract (no usado). `sources/tesis_juan.pdf`
  (25 MB) no se extrajo.
- Valores de figuras: los números del mapa salen del texto; las curvas (Yu Fig. 4, 5, 7; Langping
  Fig. 2) solo se usan por la tendencia descrita en el texto.
- Iteración 4: formulación asimétrica de la fracción resonante (defecto declarado) y parámetro libre.
- El T_on y su subida con la potencia aparecen en todas las iteraciones (resultado conocido, no cuenta).
- Muestras: Anderson calibró sobre polvo micrométrico de Lorad; los datos de contraste son de otras
  muestras (micro y nano, con y sin recubrimiento), así que los valores absolutos no son comparables
  y se priorizan tendencias y cocientes.
