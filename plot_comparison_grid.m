function plot_comparison_grid(imgs, row_labels, col_titles, clim, grid, title_str, fbase, cfg)
%PLOT_COMPARISON_GRID  Rejilla de imágenes con idéntico tamaño, escala
%   espacial, límites de color y mapa de color.
%
%   imgs        cell R x C; cada celda es una imagen npix x npix (NaN fuera
%               de la región válida) o [] para dejar la casilla vacía
%   row_labels  cell 1 x R ({} o '' para no rotular filas)
%   col_titles  cell 1 x C (solo se escriben en la primera fila)
%   clim        límites de color comunes [min max]

[R, C] = size(imgs);
if ~(all(isfinite(clim)) && clim(2) > clim(1))
    clim = [-1 1];
end
cmap = eit_colormap(cfg.fig.colormap);
f = figure('Visible', cfg.fig.visible, 'Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 3.4 * C + 3.5, 3.3 * R + 1.6]);
tl = tiledlayout(f, R, C, 'TileSpacing', 'tight', 'Padding', 'compact');
th = linspace(0, 2 * pi, 200);
last_ax = [];
for r = 1:R
    for c = 1:C
        ax = nexttile(tl);
        im = imgs{r, c};
        has_label = c == 1 && ~isempty(row_labels) && ~isempty(row_labels{r});
        if isempty(im)
            % Casilla vacía: mismo marco que las demás para que los rótulos
            % de fila queden alineados
            xlim(ax, [-1.03 1.03]); ylim(ax, [-1.03 1.03]); axis(ax, 'square');
            axis(ax, 'off');
            if has_label
                ylabel(ax, row_labels{r}, 'FontSize', 10, 'FontWeight', 'bold', 'Color', 'k');
                ax.YLabel.Visible = 'on';
            end
            continue
        end
        h = imagesc(ax, grid.centers, grid.centers, im, clim);
        set(h, 'AlphaData', ~isnan(im));
        axis(ax, 'image'); axis(ax, 'xy');
        set(ax, 'XTick', [], 'YTick', [], 'Color', 'w', 'Box', 'off', ...
            'XColor', 'none', 'YColor', 'none');
        hold(ax, 'on');
        % Contorno sobre el borde de la región evaluada (radio de la máscara)
        rc = cfg.eval.mask_radius;
        plot(ax, rc * cos(th), rc * sin(th), 'Color', [0.25 0.25 0.25], 'LineWidth', 0.8);
        xlim(ax, [-1.03 1.03]); ylim(ax, [-1.03 1.03]);
        colormap(ax, cmap);
        if r == 1 || isempty(imgs{r - 1, c})
            title(ax, col_titles{c}, 'FontSize', 9, 'FontWeight', 'normal');
        end
        if has_label
            % El eje Y es invisible ('YColor' none), pero el rótulo se colorea aparte
            ylabel(ax, row_labels{r}, 'FontSize', 10, 'FontWeight', 'bold');
            ax.YLabel.Color = 'k';
        end
        last_ax = ax;
    end
end
if ~isempty(last_ax)
    cb = colorbar(last_ax);
    cb.Layout.Tile = 'east';
    cb.Label.String = '\Delta\sigma (S/m)';
end
title(tl, title_str, 'FontSize', 10);
save_figure(f, fbase, cfg);
end
