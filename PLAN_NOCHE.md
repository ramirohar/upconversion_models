# Sesión nocturna: búsqueda de una propiedad emergente contrastable

## Contexto y alcance
- Repo `upconversion_models`, rama de trabajo nueva `noche/<fecha>` desde `unfolding`. No pushear a main.
- Base única: `AndersonModel` (src/upconversion_models/anderson.py). Toda extensión se justifica desde ahí.
- Fuera de alcance: `anderson_unfolded_assited.py`, los notebooks kinetic_onset y emergent_result,
  y el resultado ya conocido "T_on sube con la potencia". Si reaparece, se anota como conocido y no cuenta como emergente.
- `anderson_unfolded.py` se puede reutilizar como punto de partida de la iteración 1 solo si pasa la regresión de la fase 0.
- No modificar jablonski ni poincare. Las extensiones viven en src/upconversion_models
  (bloques nuevos en transitions.py).

## Reglas fijas
- Fuentes: solo vías legales (Unpaywall API, arXiv, ChemRxiv, PMC, acceso abierto del editor, PDFs en sources/).
  Si un paper no es accesible, se registra como inaccesible y se sigue.
- refs.bib: cada número va con su cita textual y su ubicación (tabla, figura o ecuación).
  Si solo se leyó el abstract, se marca `abstract-only`, y un número así no puede usarse como parámetro ni como dato de contraste.
- Cada mecanismo nuevo se apaga con un parámetro del mismo modelo (tasa = 0 o factor = 1), nunca con una clase aparte.
- La dependencia A(T) de las tasas radiativas queda constante salvo que una fuente la mida para este sistema.
- El pre-registro se commitea antes de simular. Después se commitean los resultados de la iteración y su entrada en CUADERNO_NOCHE.md.
- Presupuesto: medir en la fase 0 cuánto tarda un `SteadyState.sweep`. Ningún barrido debe pasar de unos 20 min.
  Tope total de unas 8 h de reloj; las iteraciones son las que entren.

## Fase 0: base
1. `pixi install` y una corrida de `AndersonModel`.
2. Confirmar que Er6 es el par verde (2H11/2 + 4S3/2).
3. Fijar solver y tolerancias. Guardar la línea base: poblaciones en estado estacionario y líneas de emisión
   para una grilla de potencias. La regresión compara estado estacionario con tolerancia relativa de 1e-6
   (no máximos transitorios, que dependen de save_at).
4. Mapa de observaciones (es el entregable central de la fase 0). Tabla en CUADERNO_NOCHE.md con columnas:
   observable | sistema (host, fase, dopaje, tamaño) | condiciones (T, P, λ) | valor o tendencia |
   figura o tabla | ¿acceso? | ¿usado en la calibración de Anderson?
   Prioridad: observables reportados por varios grupos de forma independiente, sobre NaYF4:Yb,Er
   en condiciones cercanas a las de Anderson 2013.

## Iteración 1 (fija): desdoblamiento y dependencia térmica
- Er6 → 6H (g = 12) y 6S (g = 4), con balance detallado k↑ = k↓ (g_H/g_S) e^(−ΔE/kT). ΔE se toma de la literatura, con su fuente.
- Tasas multifonónicas dependientes de T, ancladas a los valores de Anderson a T0 (variantes *Ref de transitions.py).
- Tests: con termalización rápida y T = T0 se recupera la línea base. Con la dependencia térmica apagada
  se recupera Anderson a cualquier T.
- Correr el barrido estándar y confrontar con el mapa de observaciones. Esto define qué filas
  el modelo ya reproduce y cuáles no.

## El ciclo (iteraciones 2 en adelante)
1. Elegir una observación: una fila del mapa que el modelo no reproduce o contradice. Prioridad para las que se repiten en
   varios trabajos y no se usaron en la calibración.
2. Búsqueda dirigida de mecanismos que puedan explicarla. Registrar las fuentes en refs.bib.
3. Propuesta: de una a tres extensiones con mecanismo, fuente, cambios en el modelo, parámetros nuevos
   (valor o rango de literatura) y cantidad de parámetros libres. Elegir la de mejor cociente entre plausibilidad
   y parámetros libres, y anotar por qué se descartaron las otras.
4. Pre-registro (commit): qué se espera en la observación elegida, qué se espera en las demás filas del mapa
   y qué sería sorprendente.
5. Implementación, tests y regresión con la extensión apagada.
6. Barrido estándar (T de 10 a 700 K, varias potencias). Observables fijos: FIR(T) y ln FIR contra 1/T,
   intensidades verde y roja contra T, pendientes log-log contra P, tiempos de vida si hay dinámica.
7. Evaluación: emergente, reproduce lo conocido, contradice la literatura o mal condicionada (ver criterios).
   Actualizar el mapa.

## Criterio de emergencia (los cuatro a la vez)
1. Aparece en el mismo observable que un dato publicado, en condiciones comparables, y ese dato no se usó para calibrar.
2. No estaba en el pre-registro, o estaba con el signo o la magnitud equivocados.
3. Desaparece al apagar el mecanismo y se sostiene en los rangos de la literatura de sus parámetros.
4. No es una consecuencia directa de una sola ecuación agregada (por ejemplo, poner una A(T) y obtener un FIR(T) no cuenta).

Escala de fuerza, para decidir entre candidatos:
reproduce cuantitativamente > reproduce la tendencia o el signo; varios grupos > un paper;
explica una anomalía o una discrepancia publicada > reproduce algo esperado.

## Semillas (orientadas a observaciones, no un guion)
- Intensidades verde y roja contra T en NaYF4:Yb,Er (apagado térmico y aumento térmico en nanopartículas).
- Cociente rojo/verde contra T y su dependencia con la potencia.
- Pendientes log-log contra T.
- Tiempo de vida de 4S3/2 y de 4F9/2 contra T.
- Desviaciones del FIR respecto de Boltzmann a baja T (Suta y Meijerink 2020). Ojo: k_therm sale de Suta,
  así que Suta no sirve como contraste.
- Mecanismos candidatos: ETU asistida por fonones (desajuste Yb–Er), (n̄+1)^p en otros gaps,
  relajación cruzada Er–Er (requiere separar site_density en fracciones de Yb y Er: costo alto).

## Salida
- CUADERNO_NOCHE.md: mapa de observaciones y, por iteración, la pregunta, las fuentes, la propuesta elegida y las descartadas,
  el pre-registro (con su hash de commit), el resultado y el estado.
- modelo_actual.md: niveles, transiciones y cada parámetro con su fuente.
- Figuras clave con el script que las genera.
- Lo más interesante: la propiedad emergente si apareció, y las dos o tres direcciones abiertas más prometedoras.
- Verificación final por un agente independiente que:
  (a) chequee que cada número tenga fuente y cita textual;
  (b) verifique con git log que cada pre-registro sea anterior a sus resultados;
  (c) vuelva a correr la ablación de la propiedad emergente desde el script.
