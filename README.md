# Estudio comparativo de algoritmos de reconstrucción EIT (MATLAB + EIDORS)

Entorno de simulación reproducible que compara seis métodos de reconstrucción
diferencial: Back Projection (Sheffield), NOSER, TSVD, Total Variation,
Gauss-Newton iterativo con Tikhonov y GREIT.

## Requisitos

- MATLAB R2020a o posterior (usa `tiledlayout`, `boxchart`, `exportgraphics`).
  No necesita toolboxes adicionales.
- EIDORS 3.x con `startup.m` ejecutado (o su ruta en `cfg.eidors_startup`).

## Uso

```matlab
run('/ruta/a/eidors/startup.m');     % o cfg.eidors_startup = '...'
cfg = config_experiment();

% 1) Primero, solo la prueba piloto
cfg.run_mode = 'pilot';
S = main_EIT_comparison(cfg);

% 2) Revisar results/<fecha>/log.txt y figures/pilot.png; después:
cfg.run_mode = 'full';
S = main_EIT_comparison(cfg);
```

Todos los parámetros están en `config_experiment.m`.

## Qué revisar en el log del piloto

| Comprobación | Dónde | Qué hacer si falla |
|---|---|---|
| Funciones de EIDORS disponibles | inicio del log | instalar o actualizar EIDORS; el estudio no sustituye algoritmos |
| Mallas `e2c2` / `c2c2` válidas | "Modelo directo ... modelo inverso" | cambiar `cfg.model.fwd_mesh` / `inv_mesh` |
| Electrodos alineados entre mallas | "desfase angular máximo" | elegir otra pareja de mallas |
| Sonda de datos nulos ≈ 0 | `sonda0` por algoritmo | la salida no es una variación de conductividad |
| Geometría del muestreo | `geom` < 1,5 px | el modelo de reconstrucción no está centrado o escalado como se supone |
| Interpretación del umbral TSVD | "TSVD umbral ..." | en EIDORS 3.12 es un porcentaje de s_1 (hipótesis 'relativo %') |
| Equivalencia 3D-2D para GREIT | "Modelo 3D para GREIT ..." | errores < 5 %; si no, revisar los electrodos del modelo extruido |
| `inv_solve_gn` en modo diferencial | "Sonda inv_solve_gn" | informativo; la implementación por defecto es la explícita |
| Iteraciones de TV | `t=... ` / columnas de `elem_data` | si `keep_iterations` no está soportado, las iteraciones quedan como NaN |

## Salidas (`results/<fecha_hora>/`)

- `log.txt`: diario completo de la ejecución.
- `results.mat`: configuración, entorno, modelos, fantomas, datos, semillas,
  hiperparámetros, curvas de ajuste, reconstrucciones, métricas y tiempos.
- `figures/*.fig` y `*.png`: comparaciones visuales, RIE por escenario, RIE
  frente al SNR, Pearson, tiempos, iteraciones y curvas de ajuste.
- `tables/*.csv` y `*.tex`: tablas de los tres experimentos (LaTeX con `booktabs`).
- `data/`: semillas, tensiones limpias y definición de fantomas.
- `config.json`, `environment.txt`.

## Estructura del código

| Fichero | Función |
|---|---|
| `main_EIT_comparison.m` | flujo completo |
| `config_experiment.m` | todos los parámetros |
| `check_environment.m` | versiones y funciones de EIDORS |
| `build_models.m` | mallas directa (fina) e inversa (gruesa), verificadas; modelo 3D de GREIT |
| `extrude_fmdl.m` | extrusión de la malla gruesa a 3D (mk_GREIT_model exige 3D) |
| `generate_phantoms.m`, `make_validation_scenarios.m` | fantomas de prueba y validación |
| `simulate_data.m`, `add_measurement_noise.m`, `build_datasets.m` | datos limpios y conjuntos con ruido |
| `make_algorithm.m` | configuración de cada algoritmo |
| `gn_tikhonov_diff.m` | Gauss-Newton iterativo diferencial con Tikhonov |
| `reconstruct_one.m`, `evaluate_sets.m` | reconstrucción y evaluación |
| `run_pilot.m`, `verify_tsvd_interpretation.m` | prueba piloto y verificaciones |
| `tune_hyperparameters.m`, `fit_calibration.m` | ajuste en validación |
| `build_inverse_models.m`, `run_reconstructions.m` | evaluación definitiva |
| `make_eval_grid.m`, `sample_on_grid.m`, `truth_on_grid.m`, `compute_metrics.m` | cuadrícula común y métricas |
| `generate_results.m`, `plot_*.m`, `write_latex_table.m`, `save_results.m` | figuras, tablas y almacenamiento |

## Definición del ruido

Los niveles se indican en porcentaje (`cfg.noise.pct = [0 1 5 15]`). Con
`cfg.noise.reference = 'difference'` (por defecto), el ruido que resulta en
dV = V1 - V0 tiene un RMS igual a ese porcentaje del RMS(dV limpio): cada
adquisición recibe ruido de desviación típica p·RMS(dV)/sqrt(2). Equivale a un
SNR sobre dV de 40, 26,0 y 16,5 dB. El log y `data/seeds.csv` registran el SNR
efectivo de cada conjunto.

## Figuras de conductividad

`cfg.fig.colormap = 'eit'` (por defecto): fondo (dsigma = 0) azul, aumentos hacia
colores cálidos y disminuciones hacia azul oscuro. `'diverging'` recupera el
mapa azul-blanco-rojo.
