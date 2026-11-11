function grid = make_eval_grid(cfg)
%MAKE_EVAL_GRID  Cuadrícula común de evaluación (npix x npix) en [-1, 1]^2.
%
%   Los centros de píxel se sitúan en el punto medio de cada celda. La
%   máscara selecciona los píxeles cuyo centro está a menos de
%   cfg.eval.mask_radius del centro: solo esos intervienen en las métricas.
%   Filas ordenadas por y creciente (representar con axis xy).

n = cfg.eval.npix;
e = linspace(-1, 1, n + 1);
c = (e(1:end-1) + e(2:end)) / 2;
[X, Y] = meshgrid(c, c);
grid.npix = n;
grid.X = X;
grid.Y = Y;
grid.centers = c;
grid.mask = (X.^2 + Y.^2) <= cfg.eval.mask_radius^2;
grid.n_valid = nnz(grid.mask);
grid.pixel_size = 2 / n;
end
