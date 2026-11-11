function [T, Timg] = truth_on_grid(fmdl, ph, grid)
%TRUTH_ON_GRID  Distribuciones reales (dsigma) llevadas a la cuadrícula común.
%   T     n_valid x n_fantomas (solo píxeles de la máscara)
%   Timg  npix x npix x n_fantomas (NaN fuera de la máscara)
T = zeros(grid.n_valid, numel(ph));
Timg = NaN(grid.npix, grid.npix, numel(ph));
for s = 1:numel(ph)
    img = sample_on_grid(fmdl, ph(s).dsigma, grid);
    v = img(grid.mask);
    if any(isnan(v))
        error('EIT:truth:mask', ['La máscara de evaluación sale de la malla directa en %s. ', ...
            'Reduzca cfg.eval.mask_radius.'], ph(s).name);
    end
    T(:, s) = v;
    img(~grid.mask) = NaN;
    Timg(:, :, s) = img;
end
end
