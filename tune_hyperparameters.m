function [hp_final, tuning] = tune_hyperparameters(models, ds, T, grid, cfg)
%TUNE_HYPERPARAMETERS  Selección de hiperparámetros en el conjunto de VALIDACIÓN.
%
%   Criterio común para todos los métodos ajustables: mínimo RIE medio
%   sobre todos los conjuntos de validación (todos los fantomas, todos los
%   niveles de ruido y realizaciones). Un candidato que produce alguna
%   reconstrucción no finita queda descartado (RIE = Inf).
%
%   Los escenarios de prueba no intervienen en ningún momento. Los valores
%   elegidos se mantienen fijos durante toda la evaluación.
%
%   Casos particulares:
%     BP    : 1) alineación de la numeración de electrodos (32 transformaciones,
%             criterio |Pearson| medio en los conjuntos de validación sin
%             ruido, invariante a la escala); 2) normalización de los datos;
%             3) calibración global de amplitud.
%     GN    : con la implementación explícita, una sola ejecución por lambda
%             con el mayor max_iter; el iterado k-ésimo es exactamente el que
%             produciría una ejecución con max_iter = k (los criterios de
%             parada no dependen de max_iter).
%     GREIT : un entrenamiento por combinación (Noise Figure, tamaño de
%             objetivo), con semilla fija, y calibración de amplitud.

tuning = struct();
hp_final = struct();
scen = [ds.scenario];

% --- Punto de control: se guarda tras cada algoritmo y, si la huella
% (rejillas, opciones y datos de validación) coincide, se reanuda.
ck = cfg.tuning.checkpoint_file;
fp = tuning_fingerprint(ds, grid, cfg);
if cfg.tuning.resume && exist(ck, 'file')
    C = load(ck);
    if isfield(C, 'fp') && isequaln(C.fp, fp)
        tuning = C.tuning; hp_final = C.hp_final;
        log_msg('Reanudando el ajuste desde %s (ya ajustados: %s)', ck, ...
            strjoin(fieldnames(hp_final)', ', '));
    else
        log_msg('Hay un punto de control en %s, pero no corresponde a esta configuración: se ignora.', ck);
    end
end

for a = 1:numel(cfg.algorithms)
    name = cfg.algorithms{a};
    if isfield(hp_final, name)
        log_msg('Ajuste de %s: recuperado del punto de control.', name);
        continue
    end
    t = tic;
    log_msg('Ajuste de %s ...', name);
    switch name
        case 'BP'
            [best, tab] = tune_bp(models, ds, T, grid, cfg);
        case 'GN'
            if strcmpi(cfg.gn.implementation, 'explicit')
                [best, tab] = tune_gn_explicit(models, ds, T, grid, cfg);
            else
                [best, tab] = tune_grid(name, models, ds, T, grid, cfg);
            end
        otherwise
            [best, tab] = tune_grid(name, models, ds, T, grid, cfg);
    end
    hp_final.(name) = best;
    tuning.(name).table = tab;
    tuning.(name).best = best;
       tuning.(name).time = toc(t);
    % Volumen de trabajo del ajuste: combinaciones evaluadas y reconstrucciones
    % realizadas (para normalizar el tiempo de ajuste entre algoritmos)
    tuning.(name).n_candidates = height(tab);
    switch name
        case 'BP'
            n_tr = 1 + (2 * cfg.model.n_elec - 1) * cfg.bp.search_orientation;
            tuning.(name).n_candidates = numel(cfg.hp.BP.normalize) * n_tr;
            tuning.(name).n_recon = numel(cfg.hp.BP.normalize) * ...
                (n_tr * nnz(isinf([ds.snr_db])) + numel(ds));
        case 'GN'
            if strcmpi(cfg.gn.implementation, 'explicit')
                tuning.(name).n_recon = numel(cfg.hp.GN.lambda) * numel(ds);
            else
                tuning.(name).n_recon = height(tab) * numel(ds);
            end
        otherwise
            tuning.(name).n_recon = height(tab) * numel(ds);
    end
    log_msg('  %s: RIE medio de validación %.4f con %s (%.1f s, %d combinaciones, %d reconstrucciones)', ...
        name, min(tab.RIE_mean), hp_to_string(best), tuning.(name).time, ...
        tuning.(name).n_candidates, tuning.(name).n_recon);
    save(ck, 'tuning', 'hp_final', 'fp', '-v7.3');
end
tuning.n_sets = numel(ds);
tuning.scenarios = unique(scen);
end

% ======================================================================
function [best, tab] = tune_grid(name, models, ds, T, grid, cfg)
G = cfg.hp.(name);
f = fieldnames(G);
vals = cellfun(@(x) G.(x), f, 'UniformOutput', false);
grids = cell(1, numel(f));
[grids{:}] = ndgrid(vals{:});
nC = numel(grids{1});
rows = cell(nC, 1);
calib = ismember(name, cfg.calibration.algorithms);
for c = 1:nC
    hp = struct();
    for i = 1:numel(f), hp.(f{i}) = grids{i}(c); end
    if strcmp(name, 'GREIT'), hp.noise_figure_used = NaN; end
    try
        A = make_algorithm(name, hp, models, cfg);
        hp = A.hp;                                   % incluye noise_figure_used en GREIT
        E = evaluate_sets(A, ds, T, grid);
        [rie, sc] = score(E, T, ds, calib);
        hp.scale = sc;
        rows{c} = row_of(hp, rie, E, A);
    catch err
        % Un candidato que falla se descarta (RIE = Inf) sin detener el ajuste.
        log_msg('  %s: candidato descartado (%s): %s', name, hp_to_string(hp), err.message);
        hp.scale = NaN;
        rows{c} = row_of(hp, Inf, [], []);
    end
end
[best, tab] = pick(rows);
end

function fp = tuning_fingerprint(ds, grid, cfg)
% Huella de todo lo que determina el resultado del ajuste.
fp.hp = cfg.hp;
fp.opts = {cfg.bp, cfg.tsvd, cfg.tv, cfg.gn, cfg.greit, cfg.calibration, cfg.model, cfg.eval, ...
           cfg.algorithms};
fp.n_sets = numel(ds);
fp.data_sum = [sum(sum([ds.vh])), sum(sum([ds.vi])), sum(sum([ds.vh].^2))];
fp.mask = find(grid.mask);
end

function [best, tab] = tune_gn_explicit(models, ds, T, grid, cfg)
lam = cfg.hp.GN.lambda;
mit = cfg.hp.GN.max_iter;
rows = {};
for i = 1:numel(lam)
    hp = struct('lambda', lam(i), 'max_iter', max(mit));
    A = make_algorithm('GN', hp, models, cfg);
    E = evaluate_sets(A, ds, T, grid, true);
    for m = mit
        Em = E;
        for k = 1:numel(ds)
            H = E.history{k};
            if isempty(H)
                d = zeros(size(A.gn.sigma0));          % ningún paso aceptado
            else
                d = H(:, min(m, size(H, 2)));
            end
            img = sample_on_grid(E.fmdl_rec, d, grid);
            r = img(grid.mask);
            mt = compute_metrics(T(:, ds(k).scenario), r);
            Em.R(:, k) = r; Em.RIE(k) = mt.RIE; Em.pearson(k) = mt.pearson;
            Em.iterations(k) = min(m, size(H, 2));
            Em.converged(k) = E.converged(k) && E.iterations(k) <= m;
        end
        [rie, ~] = score(Em, T, ds, false);
        h = struct('lambda', lam(i), 'max_iter', m, 'scale', 1);
        rows{end+1} = row_of(h, rie, Em, A); %#ok<AGROW>
    end
end
[best, tab] = pick(rows);
end

function [best, tab] = tune_bp(models, ds, T, grid, cfg)
n = cfg.model.n_elec;
clean = find(isinf([ds.snr_db]));
rows = {};
for nf = cfg.hp.BP.normalize
    if cfg.bp.search_orientation
        tr = [repmat((0:n-1)', 2, 1), [ones(n, 1); -ones(n, 1)]];
    else
        tr = [0 1];
    end
    sc = -Inf(size(tr, 1), 1);
    for i = 1:size(tr, 1)
        A = make_algorithm('BP', struct('normalize', nf, 'transform', tr(i, :), 'scale', 1), models, cfg);
        E = evaluate_sets(A, ds(clean), T, grid);
        sc(i) = abs(mean(E.pearson, 'omitnan'));
    end
    [bestcorr, ib] = max(sc);
    hp = struct('normalize', nf, 'transform', tr(ib, :), 'scale', 1);
    A = make_algorithm('BP', hp, models, cfg);
    E = evaluate_sets(A, ds, T, grid);
    [rie, s] = score(E, T, ds, true);
    hp.scale = s;
    log_msg('  BP normalize=%d: orientación k=%d s=%+d (|Pearson| %.3f), escala %.4g, RIE %.4f', ...
        nf, tr(ib, 1), tr(ib, 2), bestcorr, s, rie);
    rows{end+1} = row_of(hp, rie, E, A); %#ok<AGROW>
end
[best, tab] = pick(rows);
end

% ======================================================================
function [rie, s] = score(E, T, ds, calib)
s = 1;
if any(~isfinite(E.R(:)))
    rie = Inf; return
end
if calib
    [s, rie] = fit_calibration(T, E.R, [ds.scenario]);
else
    rie = mean(E.RIE);
end
end

function r = row_of(hp, rie, E, A)
r = hp;
if isfield(r, 'transform')
    r.transform_k = r.transform(1); r.transform_s = r.transform(2);
    r = rmfield(r, 'transform');
end
r.RIE_mean = rie;
if isempty(E)                                   % candidato fallido
    r.pearson_mean = NaN; r.time_recon_median = NaN;
    r.iterations_mean = NaN; r.time_train = NaN;
    return
end
r.pearson_mean = mean(E.pearson, 'omitnan');
r.time_recon_median = median(E.time, 'omitnan');
r.iterations_mean = mean(E.iterations, 'omitnan');
r.time_train = A.time_train;
end

function [best, tab] = pick(rows)
tab = struct2table([rows{:}]);
[~, i] = min(tab.RIE_mean);
if ~isfinite(tab.RIE_mean(i))
    error('EIT:tune:none', 'Ningún candidato produce reconstrucciones finitas.');
end
best = table2struct(tab(i, :));
if isfield(best, 'transform_k')
    best.transform = [best.transform_k best.transform_s];
end
end

function s = hp_to_string(hp)
f = fieldnames(hp);
p = {};
for i = 1:numel(f)
    v = hp.(f{i});
    if isnumeric(v) && isscalar(v) && ~ismember(f{i}, {'RIE_mean', 'pearson_mean', ...
            'time_recon_median', 'iterations_mean', 'time_train'})
        p{end+1} = sprintf('%s=%.4g', f{i}, v); %#ok<AGROW>
    end
end
s = strjoin(p, ', ');
end
