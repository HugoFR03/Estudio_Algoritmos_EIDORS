function tables = generate_results(res, ds, ph, Timg, tuning, grid, cfg, out_dir)
%GENERATE_RESULTS  Figuras y tablas del estudio (Experimentos 1, 2 y 3).
%
%   Figuras (.fig editables + .png) en out_dir/figures:
%     visual_<escenario>           comparación visual: filas = niveles de ruido
%                                  (realización 1), columnas = real + 6 algoritmos
%     rie_by_scenario_N0           RIE por algoritmo y escenario, sin ruido
%     rie_vs_noise / rie_vs_noise_by_scenario
%     pearson_by_noise             correlación de Pearson por algoritmo
%     time_reconstruction          tiempo por reconstrucción (diagrama de cajas)
%     time_setup                   preparación, precomputación y entrenamiento
%     time_total                   proceso completo de cada algoritmo, de la
%                                  construcción del modelo a una reconstrucción
%     rie_per_iteration            RIE de cada iterado de TV y GN
%     tuning_curves                RIE de validación frente a los hiperparámetros
%   Tablas (.csv; las principales también en .tex con booktabs) en out_dir/tables.

fdir = fullfile(out_dir, 'figures');
tdir = fullfile(out_dir, 'tables');
if ~exist(tdir, 'dir'), mkdir(tdir); end
nA = numel(res);
labels = {res.label};
nS = numel(ph);
snr = cfg.noise.snr_db;
nL = numel(snr);
scen = [ds.scenario];
nidx = [ds.noise_idx];
rlz = [ds.realization];
N = numel(ds);
RIE = vertcat(res.RIE);
PR  = vertcat(res.pearson);
TM  = vertcat(res.time);
IT  = vertcat(res.iterations);
CV  = vertcat(res.converged);
if isfield(cfg.noise, 'pct')
    noise_names = arrayfun(@noise_name_pct, cfg.noise.pct, 'UniformOutput', false);
else
    noise_names = arrayfun(@noise_name, snr, 'UniformOutput', false);
end
mk = {'o', 's', 'd', '^', 'v', 'p'};

%% ---------------- A. Comparaciones visuales: un escenario por figura,
%   una fila por nivel de ruido (realización 1), mismos límites de color
for s = 1:nS
    lim = max(max(abs(Timg(:, :, s)), [], 'omitnan'));
    C = cell(nL, nA + 1);
    C{1, 1} = Timg(:, :, s);                     % la distribución real, una sola vez
    for n = 1:nL
        k = find(scen == s & nidx == n & rlz == 1, 1);
        if isempty(k), continue, end
        for a = 1:nA
            C{n, a + 1} = res(a).images(:, :, k);
        end
    end
    row_labels = noise_names;
    row_labels{1} = 'Sin ruido';
    plot_comparison_grid(C, row_labels, [{'Real'}, labels], [-lim lim], grid, ...
        sprintf('%s (realización 1 de cada nivel de ruido)', ph(s).label), ...
        fullfile(fdir, sprintf('visual_%s', ph(s).name)), cfg);
end

%% ---------------- B1. RIE por algoritmo y escenario (sin ruido)
if ~any(isinf(snr))
    error('EIT:results:N0', 'El Experimento 1 necesita un nivel sin ruido (SNR = Inf) en cfg.noise.snr_db.');
end
k0 = arrayfun(@(s) find(scen == s & isinf([ds.snr_db]), 1), 1:nS);
M = RIE(:, k0)';                                  % nS x nA
f = new_fig(cfg, [2 2 20 9]);
bar(M);
set(gca, 'XTickLabel', {ph.label}, 'XTickLabelRotation', 15);
ylabel('RIE'); grid_on();
legend(labels, 'Location', 'eastoutside');
title('Experimento 1: RIE sin ruido por escenario');
save_figure(f, fullfile(fdir, 'rie_by_scenario_N0'), cfg);

%% ---------------- B2. RIE frente al nivel de ruido
[muR, sdR] = agg_by(RIE, scen, nidx, nS, nL);     % nA x nS x nL
f = new_fig(cfg, [2 2 22 16]);
tl = tiledlayout(f, 2, ceil(nS / 2), 'TileSpacing', 'compact');
for s = 1:nS
    ax = nexttile(tl); hold(ax, 'on');
    for a = 1:nA
        errorbar(ax, (1:nL) + (a - 3.5) * 0.04, squeeze(muR(a, s, :)), squeeze(sdR(a, s, :)), ...
            ['-' mk{mod(a - 1, numel(mk)) + 1}], 'MarkerSize', 4, 'LineWidth', 1);
    end
    set(ax, 'XTick', 1:nL, 'XTickLabel', noise_names); xlim(ax, [0.6 nL + 0.4]);
    ylabel(ax, 'RIE'); title(ax, ph(s).label, 'FontWeight', 'normal');
    set(ax, 'XGrid', 'on', 'YGrid', 'on');   % aquí 'grid' es la cuadrícula de evaluación, no la función
end
lg = legend(ax, labels); lg.Layout.Tile = 'east';
xl = 'Nivel de ruido (RMS relativo a \DeltaV)';
title(tl, 'Experimento 2: RIE (media \pm desviación típica) frente al nivel de ruido');
xlabel(tl, xl);
save_figure(f, fullfile(fdir, 'rie_vs_noise_by_scenario'), cfg);

[muRp, sdRp] = agg_pooled(RIE, nidx, nL);         % nA x nL
f = new_fig(cfg, [2 2 16 10]); hold on
for a = 1:nA
    errorbar((1:nL) + (a - 3.5) * 0.04, muRp(a, :), sdRp(a, :), ...
        ['-' mk{mod(a - 1, numel(mk)) + 1}], 'MarkerSize', 5, 'LineWidth', 1.2);
end
set(gca, 'XTick', 1:nL, 'XTickLabel', noise_names); xlim([0.6 nL + 0.4]);
ylabel('RIE'); xlabel(xl); grid_on(); legend(labels, 'Location', 'eastoutside');
title('Experimento 2: RIE medio de los cuatro escenarios frente al nivel de ruido');
save_figure(f, fullfile(fdir, 'rie_vs_noise'), cfg);

%% ---------------- B3. Correlación de Pearson por algoritmo
[muP, sdP] = agg_pooled(PR, nidx, nL);
f = new_fig(cfg, [2 2 20 10]);
b = bar(muP');                                    % grupos = niveles de ruido
hold on
for a = 1:nA
    try
        errorbar(b(a).XEndPoints, muP(a, :), sdP(a, :), 'k.', 'LineWidth', 0.8);
    catch
    end
end
set(gca, 'XTickLabel', noise_names); ylim([-1 1]); ylabel('Correlación de Pearson');
grid_on(); legend(labels, 'Location', 'eastoutside');
title('Correlación de Pearson por algoritmo y nivel de ruido (media \pm desviación típica)');
save_figure(f, fullfile(fdir, 'pearson_by_noise'), cfg);

%% ---------------- B4. Tiempo de reconstrucción
f = new_fig(cfg, [2 2 16 10]);
g = categorical(repelem(labels(:), N, 1), labels);
tv = reshape(TM', [], 1);
boxchart(g, tv);
set(gca, 'YScale', 'log'); ylabel('Tiempo por reconstrucción (s)'); grid_on();
title(sprintf('Experimento 3: tiempo de reconstrucción (%d reconstrucciones por algoritmo)', N));
save_figure(f, fullfile(fdir, 'time_reconstruction'), cfg);

%% ---------------- B5. Preparación (izquierda) y ajuste de hiperparámetros (derecha)
%   El ajuste va en un panel propio y en escala logarítmica: es varios
%   órdenes de magnitud mayor y depende del tamaño de la rejilla de búsqueda.
[t_tune, n_cand, n_rtune] = tuning_cost(tuning, {res.name});
S = [[res.time_build]', [res.time_precompute]', [res.time_train]'];
f = new_fig(cfg, [2 2 28 10]);
tl = tiledlayout(f, 1, 2, 'TileSpacing', 'compact');
ax = nexttile(tl);
bar(ax, S, 'stacked');
set(ax, 'XTickLabel', labels); ylabel(ax, 'Tiempo (s)');
set(ax, 'XGrid', 'on', 'YGrid', 'on', 'Box', 'on');
legend(ax, {'Configuración del modelo', 'Precomputación (1.ª llamada)', 'Entrenamiento (GREIT)'}, ...
    'Location', 'northwest');
title(ax, 'Preparación con hiperparámetros fijos', 'FontWeight', 'normal');
ax = nexttile(tl);
bar(ax, t_tune, 'FaceColor', [0.45 0.45 0.45]);
set(ax, 'XTickLabel', labels, 'YScale', 'log'); ylabel(ax, 'Tiempo (s)');
set(ax, 'XGrid', 'on', 'YGrid', 'on', 'Box', 'on');
for a = 1:nA
    if isfinite(t_tune(a))
        text(ax, a, t_tune(a) * 1.15, sprintf('%d comb.\n%d rec.', n_cand(a), n_rtune(a)), ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 7);
    end
end
if any(t_tune > 0), ylim(ax, [min(t_tune(t_tune > 0)) / 2, max(t_tune) * 6]); end
title(ax, 'Ajuste de hiperparámetros en validación (una sola vez)', 'FontWeight', 'normal');
title(tl, 'Experimento 3: coste de preparación');
save_figure(f, fullfile(fdir, 'time_setup'), cfg);
%% ---------------- B6. Error de imagen por iteración (TV y GN)
%   Media del RIE de cada iterado sobre las realizaciones. Si una ejecución
%   se detiene antes (convergencia), su imagen ya no cambia: se prolonga con
%   su último valor hasta la iteración máxima de ese algoritmo.
iterA = find(arrayfun(@(r) isfield(r, 'rie_iter') && any(~cellfun(@isempty, r.rie_iter)), res));
RI = table();
if ~isempty(iterA)
    f = new_fig(cfg, [2 2 6 * nS + 4, 6 * numel(iterA) + 1.5]);
    tl = tiledlayout(f, numel(iterA), nS, 'TileSpacing', 'compact');
    cols = lines(nL);
    for j = 1:numel(iterA)
        a = iterA(j);
        kmax = max(cellfun(@numel, res(a).rie_iter));
        for s = 1:nS
            ax = nexttile(tl); hold(ax, 'on');
            for n = 1:nL
                sel = find(scen == s & nidx == n);
                M_it = NaN(numel(sel), kmax);
                for q = 1:numel(sel)
                    ri = res(a).rie_iter{sel(q)};
                    if isempty(ri), continue, end
                    M_it(q, :) = [ri, repmat(ri(end), 1, kmax - numel(ri))];
                end
                mu = mean(M_it, 1, 'omitnan'); sd = std(M_it, 0, 1, 'omitnan');
                plot(ax, 1:kmax, mu, '-o', 'Color', cols(n, :), 'MarkerSize', 3, ...
                    'MarkerFaceColor', cols(n, :), 'LineWidth', 1.1, 'DisplayName', noise_names{n});
                RI = [RI; table(repmat(string(res(a).name), kmax, 1), repmat(string(ph(s).name), kmax, 1), ...
                    repmat(string(noise_names{n}), kmax, 1), (1:kmax)', mu(:), sd(:), ...
                    'VariableNames', {'algoritmo', 'escenario', 'ruido', 'iteracion', 'RIE_media', 'RIE_std'})]; %#ok<AGROW>
            end
            set(ax, 'XGrid', 'on', 'YGrid', 'on', 'Box', 'on');
            xlim(ax, [0.5 kmax + 0.5]);
            if j == 1, title(ax, ph(s).label, 'FontWeight', 'normal', 'FontSize', 9); end
            if s == 1, ylabel(ax, sprintf('%s: RIE', res(a).label)); end
            if j == numel(iterA), xlabel(ax, 'Iteración'); end
        end
    end
    lg = legend(ax, 'Location', 'eastoutside'); lg.Layout.Tile = 'east';
    title(tl, 'Error de imagen de cada iterado (media de las realizaciones)');
    save_figure(f, fullfile(fdir, 'rie_per_iteration'), cfg);
    writetable(RI, fullfile(tdir, 'rie_per_iteration.csv'));
end

%% ---------------- B8. Proceso completo de cada algoritmo
%   modelo FEM + configuración + entrenamiento + precomputación + una
%   reconstrucción media. No incluye el ajuste de hiperparámetros, que es
%   una fase de diseño y no se repite al usar el algoritmo.
t_model = zeros(nA, 1);
for a = 1:nA
    if isfield(res, 'time_model') && ~isempty(res(a).time_model), t_model(a) = res(a).time_model; end
end
TT = [t_model, [res.time_build]', [res.time_train]', [res.time_precompute]', mean(TM, 2, 'omitnan')];
tot = sum(TT, 2);
[~, o] = sort(tot, 'ascend');
f = new_fig(cfg, [2 2 18 10]);
bh = barh(1:nA, TT(o, :), 'stacked');
set(gca, 'YTick', 1:nA, 'YTickLabel', labels(o));
hold on
for q = 1:nA
    text(tot(o(q)), q, sprintf('  %.3g s', tot(o(q))), 'VerticalAlignment', 'middle', 'FontSize', 8);
end
xlim([0 max(tot) * 1.18]);
xlabel('Tiempo (s)'); grid_on();
legend(bh, {'Modelo FEM', 'Configuración', 'Entrenamiento (GREIT)', 'Precomputación', ...
    'Reconstrucción (media)'}, 'Location', 'southeast');
title('Experimento 3: tiempo hasta la primera imagen (hiperparámetros ya ajustados)');
save_figure(f, fullfile(fdir, 'time_total'), cfg);

%% ---------------- B7. Curvas de ajuste en validación
plot_tuning_curves(tuning, cfg, fullfile(fdir, 'tuning_curves'));

%% ---------------- Tablas
% Larga: una fila por conjunto y algoritmo
nR = N * nA;
L = table('Size', [nR 15], 'VariableTypes', {'string', 'string', 'string', 'double', 'double', ...
    'double', 'double', 'double', 'string', 'double', 'double', 'double', 'double', 'logical', 'string'}, ...
    'VariableNames', {'algorithm', 'scenario', 'noise', 'snr_db', 'realization', 'seed_v0', ...
    'seed_v1', 'snr_eff_dv_db', 'hp', 'RIE', 'pearson', 'time_s', 'iterations', 'converged', 'stop_reason'});
i = 0;
for a = 1:nA
    hs = hp_string(res(a).hp);
    for k = 1:N
        i = i + 1;
        L(i, :) = {string(res(a).name), string(ds(k).scenario_name), string(ds(k).noise_label), ds(k).snr_db, ...
            ds(k).realization, ds(k).seed_v0, ds(k).seed_v1, ds(k).snr_eff_dv_db, string(hs), ...
            RIE(a, k), PR(a, k), TM(a, k), IT(a, k), CV(a, k) == 1, string(res(a).stop_reason{k})};
    end
end
writetable(L, fullfile(tdir, 'metrics_long.csv'));

% Experimento 1: sin ruido, escenario x algoritmo
E1 = table();
for s = 1:nS
    for a = 1:nA
        E1 = [E1; table(string(ph(s).label), string(res(a).label), RIE(a, k0(s)), PR(a, k0(s)), ...
            'VariableNames', {'escenario', 'algoritmo', 'RIE', 'pearson'})]; %#ok<AGROW>
    end
end
writetable(E1, fullfile(tdir, 'exp1_noiseless.csv'));
W = array2table(M, 'VariableNames', matlab.lang.makeValidName({res.name}));
W = [table(string({ph.label})', 'VariableNames', {'escenario'}), W];
write_latex_table(W, fullfile(tdir, 'exp1_rie_noiseless.tex'), ...
    'RIE sin ruido por escenario y algoritmo (Experimento 1).', 'tab:exp1_rie', '%.3f');
WP = array2table(PR(:, k0)', 'VariableNames', matlab.lang.makeValidName({res.name}));
WP = [table(string({ph.label})', 'VariableNames', {'escenario'}), WP];
write_latex_table(WP, fullfile(tdir, 'exp1_pearson_noiseless.tex'), ...
    'Correlación de Pearson sin ruido por escenario y algoritmo (Experimento 1).', 'tab:exp1_pearson', '%.3f');

% Experimento 2: por escenario y nivel de ruido
[muP2, sdP2] = agg_by(PR, scen, nidx, nS, nL);
snr_eff = arrayfun(@(s, n) mean([ds(scen == s & nidx == n).snr_eff_dv_db]), ...
    repmat((1:nS)', 1, nL), repmat(1:nL, nS, 1));
E2 = table();
for s = 1:nS
    for n = 1:nL
        nk = nnz(scen == s & nidx == n);
        for a = 1:nA
            E2 = [E2; table(string(ph(s).label), string(noise_names{n}), snr(n), snr_eff(s, n), ...
                string(res(a).label), nk, muR(a, s, n), sdR(a, s, n), muP2(a, s, n), sdP2(a, s, n), ...
                'VariableNames', {'escenario', 'ruido', 'snr_db', 'snr_eff_dv_db', 'algoritmo', ...
                'n', 'RIE_media', 'RIE_std', 'pearson_media', 'pearson_std'})]; %#ok<AGROW>
        end
    end
end
writetable(E2, fullfile(tdir, 'exp2_noise_by_scenario.csv'));
c = cell(nL, nA);
for n = 1:nL
    for a = 1:nA
        c{n, a} = sprintf('%.3f $\\pm$ %.3f', muRp(a, n), sdRp(a, n));
    end
end
W2 = cell2table(c, 'VariableNames', matlab.lang.makeValidName({res.name}));
W2 = [table(string(noise_names(:)), 'VariableNames', {'ruido'}), W2];
writetable(W2, fullfile(tdir, 'exp2_rie_pooled.csv'));
write_latex_table(W2, fullfile(tdir, 'exp2_rie_pooled.tex'), ...
    'RIE medio $\pm$ desviación típica de los cuatro escenarios por nivel de ruido (Experimento 2).', ...
    'tab:exp2_rie', '%.3f');

% Experimento 3: tiempos e iteraciones
E3 = table();
for a = 1:nA
    t = TM(a, :); t = t(isfinite(t));
    E3 = [E3; table(string(res(a).label), t_model(a), res(a).time_build, res(a).time_precompute, ...
        res(a).time_train, mean(t), median(t), std(t), pct(t, 25), pct(t, 75), tot(a), ...
        mean(IT(a, :), 'omitnan'), conv_rate(CV(a, :)), t_tune(a), n_cand(a), n_rtune(a), ...
        t_tune(a) / n_rtune(a), ...
        'VariableNames', {'algoritmo', 't_modelo_s', 't_config_s', 't_precomp_s', 't_entren_s', 't_rec_media_s', ...
        't_rec_mediana_s', 't_rec_std_s', 't_rec_p25_s', 't_rec_p75_s', 't_total_s', 'iter_media', 'tasa_conv', ...
        't_ajuste_s', 'n_combinaciones', 'n_rec_ajuste', 't_ajuste_por_rec_s'})]; %#ok<AGROW>
end
writetable(E3, fullfile(tdir, 'exp3_timing.csv'));
write_latex_table(E3(:, [1 3 4 5 7 8 11 12]), fullfile(tdir, 'exp3_timing.tex'), ...
    'Coste computacional por algoritmo (Experimento 3). Tiempos en segundos.', 'tab:exp3_time', '%.3g');

% Hiperparámetros seleccionados
HP = table(string({res.label})', string(arrayfun(@(r) hp_string(r.hp), res, 'UniformOutput', false))', ...
    [res.scale]', [res.normalize]', 'VariableNames', {'algoritmo', 'hiperparametros', 'escala', 'normalizacion'});
writetable(HP, fullfile(tdir, 'hyperparameters.csv'));
write_latex_table(HP(:, 1:2), fullfile(tdir, 'hyperparameters.tex'), ...
    'Hiperparámetros seleccionados en el conjunto de validación.', 'tab:hp', '%.3g');

% Curvas de ajuste
fn = fieldnames(tuning);
for i = 1:numel(fn)
    if isstruct(tuning.(fn{i})) && isfield(tuning.(fn{i}), 'table')
        writetable(tuning.(fn{i}).table, fullfile(tdir, sprintf('tuning_%s.csv', fn{i})));
    end
end

tables = struct('long', L, 'exp1', E1, 'exp2', E2, 'exp2_pooled', W2, 'exp3', E3, 'hp', HP);
log_msg('Figuras en %s y tablas en %s', fdir, tdir);
end

% ======================================================================
function f = new_fig(cfg, pos)
f = figure('Visible', cfg.fig.visible, 'Color', 'w', 'Units', 'centimeters', 'Position', pos);
end

function grid_on()
grid on; box on;
end

function s = noise_name(x)
if isinf(x), s = 'Sin ruido'; else, s = sprintf('%g dB', x); end
end

function s = noise_name_pct(p)
if p == 0, s = 'Sin ruido'; else, s = sprintf('%g %%', p); end
end

function [mu, sd] = agg_by(X, scen, nidx, nS, nL)
nA = size(X, 1);
mu = NaN(nA, nS, nL); sd = NaN(nA, nS, nL);
for s = 1:nS
    for n = 1:nL
        sel = scen == s & nidx == n;
        if ~any(sel), continue, end
        mu(:, s, n) = mean(X(:, sel), 2, 'omitnan');
        sd(:, s, n) = std(X(:, sel), 0, 2, 'omitnan');
    end
end
end

function [mu, sd] = agg_pooled(X, nidx, nL)
nA = size(X, 1);
mu = NaN(nA, nL); sd = NaN(nA, nL);
for n = 1:nL
    sel = nidx == n;
    mu(:, n) = mean(X(:, sel), 2, 'omitnan');
    sd(:, n) = std(X(:, sel), 0, 2, 'omitnan');
end
end

function [t, nc, nr] = tuning_cost(tuning, names)
% Tiempo total de ajuste, combinaciones evaluadas y reconstrucciones realizadas.
% Resultados antiguos sin recuento: combinaciones = filas de la tabla de ajuste.
n = numel(names);
t = NaN(n, 1); nc = NaN(n, 1); nr = NaN(n, 1);
for a = 1:n
    if ~isfield(tuning, names{a}), continue, end
    q = tuning.(names{a});
    if isfield(q, 'time'), t(a) = q.time; end
    if isfield(q, 'n_candidates'), nc(a) = q.n_candidates;
    elseif isfield(q, 'table'), nc(a) = height(q.table);
    end
    if isfield(q, 'n_recon'), nr(a) = q.n_recon; end
end
end

function r = conv_rate(c)
% Fracción de ejecuciones convergidas; NaN para métodos no iterativos.
c = c(~isnan(c));
if isempty(c), r = NaN; else, r = mean(c == 1); end
end

function v = pct(x, p)
% Percentil con interpolación lineal (sin Statistics Toolbox).
x = sort(x(:));
if isempty(x), v = NaN; return, end
pos = 1 + (numel(x) - 1) * p / 100;
lo = floor(pos); hi = ceil(pos);
v = x(lo) + (pos - lo) * (x(hi) - x(lo));
end
