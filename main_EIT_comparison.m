function S = main_EIT_comparison(cfg)
%MAIN_EIT_COMPARISON  Estudio comparativo de seis algoritmos de reconstrucción
%   EIT diferencial (BP, NOSER, TSVD, TV, GN + Tikhonov, GREIT) en MATLAB +
%   EIDORS, por simulación numérica.
%
%   S = main_EIT_comparison()        ejecuta con config_experiment()
%   S = main_EIT_comparison(cfg)     ejecuta con una configuración propia
%
%   Con cfg.run_mode = 'pilot' se detiene tras la prueba piloto.
%
%   Flujo:
%     0. Entorno: versiones de MATLAB y EIDORS, funciones disponibles.
%     1. Modelos FEM: directo fino y inverso grueso, verificados.
%     2. Fantomas de prueba (4) y de validación (aleatorios, independientes).
%     3. Datos limpios (problema directo en la malla fina) y conjuntos con
%        ruido, generados UNA vez con semillas registradas.
%     4. Prueba piloto (inclusión central, sin ruido, seis algoritmos).
%     5. Ajuste de hiperparámetros en validación (mínimo RIE medio).
%     6. Modelos definitivos y entrenamiento de GREIT.
%     7. Experimentos 1-3 sobre los 124 conjuntos de prueba.
%     8. Figuras, tablas y almacenamiento.
%
%   Todos los resultados quedan en cfg.output_root/<fecha_hora>/.

if nargin < 1 || isempty(cfg)
    cfg = config_experiment();
end
if ~isfield(cfg, 'run_mode'), cfg.run_mode = 'full'; end

out_dir = fullfile(cfg.output_root, datestr(now, 'yyyymmdd_HHMMSS')); %#ok<TNOW1,DATST>
mkdir(out_dir); mkdir(fullfile(out_dir, 'figures'));
diary(fullfile(out_dir, 'log.txt'));
cleanup = onCleanup(@() diary('off'));
t_total = tic;
log_msg('=== Estudio comparativo EIT | resultados en %s', out_dir);

%% 0. Entorno
env = check_environment(cfg);

%% 1. Modelos
log_msg('--- Modelos FEM');
models = build_models(cfg);
grid = make_eval_grid(cfg);

%% 2-3. Fantomas y datos
log_msg('--- Fantomas y datos de prueba');
ph_test  = generate_phantoms(models.fwd, cfg.scenarios, cfg.model.sigma_bkg);
sim_test = simulate_data(models.fwd, ph_test, cfg.model.sigma_bkg);
ds_test  = build_datasets(sim_test, cfg.noise.snr_db, cfg.noise.labels, ...
    cfg.noise.n_realizations, cfg.noise.seed_base, cfg.noise.reference);
log_msg('Conjuntos de prueba: %d (%d escenarios)', numel(ds_test), numel(ph_test));
report_effective_snr(ds_test, cfg);

log_msg('--- Fantomas y datos de validación');
val_sc  = make_validation_scenarios(cfg);
ph_val  = generate_phantoms(models.fwd, val_sc, cfg.model.sigma_bkg);
sim_val = simulate_data(models.fwd, ph_val, cfg.model.sigma_bkg);
ds_val  = build_datasets(sim_val, cfg.noise.snr_db, cfg.noise.labels, ...
    cfg.val.n_realizations, cfg.val.seed_base, cfg.noise.reference);
log_msg('Conjuntos de validación: %d (%d fantomas)', numel(ds_val), numel(ph_val));
assert(isempty(intersect([[ds_test.seed_v0] [ds_test.seed_v1]], [[ds_val.seed_v0] [ds_val.seed_v1]])), ...
    'EIT:main:seeds', 'Las semillas de prueba y validación se solapan.');

%% 4. Prueba piloto
log_msg('--- Prueba piloto');
[pilot, grid] = run_pilot(models, sim_test, ph_test, grid, cfg, out_dir);
S = struct('cfg', cfg, 'env', env, 'models', models, 'grid', grid, 'ph_test', ph_test, ...
    'sim_test', sim_test, 'ds_test', ds_test, 'ph_val', ph_val, 'sim_val', sim_val, ...
    'ds_val', ds_val, 'pilot', pilot);
if ~pilot.passed && cfg.pilot.stop_on_fail
    save(fullfile(out_dir, 'pilot_failed.mat'), '-struct', 'S', '-v7.3');
    error('EIT:main:pilot', ['La prueba piloto no se ha superado. Revise %s/log.txt y ', ...
        'pilot_failed.mat antes de lanzar el experimento completo.'], out_dir);
end
if strcmpi(cfg.run_mode, 'pilot')
    save(fullfile(out_dir, 'pilot.mat'), '-struct', 'S', '-v7.3');
    log_msg('Modo piloto: fin (%.1f s).', toc(t_total));
    return
end

% Distribuciones reales en la cuadrícula común (máscara ya definitiva)
[T_test, Timg_test] = truth_on_grid(models.fwd, ph_test, grid);
T_val = truth_on_grid(models.fwd, ph_val, grid);

%% 5. Ajuste de hiperparámetros (solo validación)
log_msg('--- Ajuste de hiperparámetros (%d conjuntos de validación)', numel(ds_val));
[hp_final, tuning] = tune_hyperparameters(models, ds_val, T_val, grid, cfg);

%% 6. Modelos definitivos
log_msg('--- Modelos definitivos');
algs = build_inverse_models(models, hp_final, cfg, pilot);

%% 7. Experimentos 1-3
log_msg('--- Evaluación (%d conjuntos x %d algoritmos)', numel(ds_test), numel(algs));
res = run_reconstructions(algs, ds_test, T_test, grid, cfg);

%% 8. Resultados: primero se guardan los datos, después las figuras
S.grid = grid;
S.T_test = T_test; S.Timg_test = Timg_test; S.T_val = T_val;
S.hp_final = hp_final; S.tuning = tuning; S.algs = algs; S.res = res;
S.time_total = toc(t_total);
save_results(out_dir, S);
log_msg('--- Figuras y tablas');
try
    tables = generate_results(res, ds_test, ph_test, Timg_test, tuning, grid, cfg, out_dir);
    save(fullfile(out_dir, 'results.mat'), 'tables', '-append');
    S.tables = tables;
catch err
    log_msg('Error al generar figuras o tablas: %s', err.message);
    log_msg('Las reconstrucciones y métricas ya están en results.mat. Tras corregir, ejecute:');
    log_msg('  regenerate_results(''%s'')', out_dir);
    rethrow(err);
end
log_msg('=== Fin: %.1f min', S.time_total / 60);
end

% ======================================================================
function report_effective_snr(ds, cfg)
% SNR efectivo sobre dV (el que realmente ven los algoritmos).
for n = 1:numel(cfg.noise.snr_db)
    if isinf(cfg.noise.snr_db(n)), continue, end
    sel = [ds.noise_idx] == n;
    log_msg('  %s (%g dB por adquisición): SNR efectivo en dV entre %.1f y %.1f dB', ...
        cfg.noise.labels{n}, cfg.noise.snr_db(n), min([ds(sel).snr_eff_dv_db]), ...
        max([ds(sel).snr_eff_dv_db]));
end
end
