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
