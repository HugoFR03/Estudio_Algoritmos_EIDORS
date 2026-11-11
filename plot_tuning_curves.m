function plot_tuning_curves(tuning, cfg, fbase)
%PLOT_TUNING_CURVES  RIE medio de validación frente a los hiperparámetros.
%   Un panel por algoritmo: el primer hiperparámetro en abscisas (escala
%   logarítmica si es positivo y abarca más de una década) y una curva por
%   valor del segundo. El punto elegido se marca con un círculo relleno.

metric = {'RIE_mean', 'pearson_mean', 'time_recon_median', 'iterations_mean', ...
          'time_train', 'scale', 'transform_k', 'transform_s'};
algs = cfg.algorithms(isfield(tuning, cfg.algorithms));
f = figure('Visible', cfg.fig.visible, 'Color', 'w', 'Units', 'centimeters', 'Position', [2 2 24 14]);
tl = tiledlayout(f, 2, ceil(numel(algs) / 2), 'TileSpacing', 'compact');
for i = 1:numel(algs)
    T = tuning.(algs{i}).table;
    p = setdiff(T.Properties.VariableNames, metric, 'stable');
    ax = nexttile(tl); hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    if isempty(p)
        bar(ax, T.RIE_mean);
    else
        x = T.(p{1});
        if numel(p) >= 2, g = T.(p{2}); else, g = ones(height(T), 1); end
        ug = unique(g);
        for j = 1:numel(ug)
            sel = g == ug(j);
            [xs, o] = sort(x(sel));
            y = T.RIE_mean(sel); y = y(o);
            if numel(p) >= 2
                plot(ax, xs, y, '-o', 'MarkerSize', 3, 'DisplayName', sprintf('%s=%g', p{2}, ug(j)));
            else
                plot(ax, xs, y, '-o', 'MarkerSize', 3, 'HandleVisibility', 'off');
            end
        end
        [~, ib] = min(T.RIE_mean);
        plot(ax, x(ib), T.RIE_mean(ib), 'ko', 'MarkerFaceColor', 'k', 'HandleVisibility', 'off');
        if all(x > 0) && max(x) / min(x) > 10, set(ax, 'XScale', 'log'); end
        xlabel(ax, strrep(p{1}, '_', '\_'));
        if numel(p) >= 2, legend(ax, 'Location', 'best', 'FontSize', 6); end
    end
    ylabel(ax, 'RIE medio (validación)');
    title(ax, algs{i}, 'FontWeight', 'normal');
end
title(tl, 'Ajuste de hiperparámetros en el conjunto de validación');
save_figure(f, fbase, cfg);
end
