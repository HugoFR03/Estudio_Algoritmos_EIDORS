function save_results(out_dir, S)
%SAVE_RESULTS  Almacena todo lo necesario para reproducir y auditar el estudio.
%
%   out_dir/results.mat      configuración, entorno, modelos, fantomas,
%                            datos simulados, conjuntos (con semillas),
%                            hiperparámetros, ajuste, reconstrucciones,
%                            métricas y tiempos
%   out_dir/data/seeds.csv   semillas de cada conjunto de prueba y validación
%   out_dir/data/*.csv       tensiones simuladas limpias y fantomas
%   out_dir/config.json      configuración legible
%   out_dir/environment.txt  versiones y funciones de EIDORS comprobadas

ddir = fullfile(out_dir, 'data');
if ~exist(ddir, 'dir'), mkdir(ddir); end

save(fullfile(out_dir, 'results.mat'), '-struct', 'S', '-v7.3');

% Semillas
seeds = [seed_table(S.ds_test, "prueba"); seed_table(S.ds_val, "validacion")];
writetable(seeds, fullfile(ddir, 'seeds.csv'));

% Datos simulados limpios
writematrix([S.sim_test.V0, S.sim_test.V1], fullfile(ddir, 'voltages_test_clean.csv'));
writematrix([S.sim_val.V0, S.sim_val.V1], fullfile(ddir, 'voltages_validation_clean.csv'));

% Fantomas (definición geométrica)
writetable(phantom_table(S.ph_test), fullfile(ddir, 'phantoms_test.csv'));
writetable(phantom_table(S.ph_val), fullfile(ddir, 'phantoms_validation.csv'));

% Configuración
try
    txt = jsonencode(S.cfg, 'PrettyPrint', true);   % R2021a o posterior
catch
    txt = jsonencode(S.cfg);
end
fid = fopen(fullfile(out_dir, 'config.json'), 'w', 'n', 'UTF-8');
fprintf(fid, '%s', txt);
fclose(fid);

% Entorno
fid = fopen(fullfile(out_dir, 'environment.txt'), 'w', 'n', 'UTF-8');
fprintf(fid, 'Fecha: %s\nMATLAB %s (%s)\nEIDORS %s\n\n', S.env.timestamp, ...
    S.env.matlab_release, S.env.matlab_version, S.env.eidors_version);
fn = fieldnames(S.env.functions);
for i = 1:numel(fn)
    e = S.env.functions.(fn{i});
    fprintf(fid, '==== %s  [%s]\n%s\n%s\n\n', fn{i}, tern(e.available, 'disponible', 'NO disponible'), ...
        e.path, e.help_head);
end
fclose(fid);
log_msg('Resultados guardados en %s', out_dir);
end

function T = seed_table(ds, set_name)
n = numel(ds);
T = table(repmat(set_name, n, 1), string({ds.scenario_name})', string({ds.noise_label})', ...
    [ds.snr_db]', [ds.realization]', [ds.seed_v0]', [ds.seed_v1]', [ds.sigma_v0]', ...
    [ds.sigma_v1]', [ds.snr_eff_dv_db]', 'VariableNames', {'conjunto', 'escenario', 'ruido', ...
    'snr_db', 'realizacion', 'semilla_V0', 'semilla_V1', 'sigma_V0', 'sigma_V1', 'snr_efectivo_dV_db'});
end

function T = phantom_table(ph)
C = cell(0, 7);
for s = 1:numel(ph)
    for k = 1:numel(ph(s).inclusions)
        q = ph(s).inclusions(k);
        C(end+1, :) = {string(ph(s).name), k, q.center(1), q.center(2), q.radius, q.sigma, ...
            ph(s).n_elems_inclusion(k)}; %#ok<AGROW>
    end
end
T = cell2table(C, 'VariableNames', {'escenario', 'inclusion', 'x', 'y', 'radio', 'sigma_Sm', ...
    'elementos_malla_fina'});
end

function o = tern(c, a, b)
if c, o = a; else, o = b; end
end
