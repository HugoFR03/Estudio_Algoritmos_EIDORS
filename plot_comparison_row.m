function plot_comparison_row(imgs, titles, clim, grid, title_str, fbase, cfg)
%PLOT_COMPARISON_ROW  Fila de imágenes comparables (usada en la prueba piloto).
%   imgs: npix x npix x K. Es un caso particular de plot_comparison_grid.
K = size(imgs, 3);
C = cell(1, K);
for k = 1:K
    C{k} = imgs(:, :, k);
end
plot_comparison_grid(C, {}, titles, clim, grid, title_str, fbase, cfg);
end
