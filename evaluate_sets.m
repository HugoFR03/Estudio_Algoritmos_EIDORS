function E = evaluate_sets(A, ds, T, grid, keep_history, track_iter)
%EVALUATE_SETS  Reconstruye y evalúa un algoritmo sobre una lista de conjuntos.
%
%   E = evaluate_sets(A, ds, T, grid) devuelve
%     E.R        n_valid x N, reconstrucciones en la cuadrícula común
%     E.RIE, E.pearson (1 x N)
%     E.time, E.iterations, E.converged (1 x N), E.stop_reason (1 x N cell)
%     E.history  (cell) iterados, si keep_history
%     E.rie_iter (cell) RIE de cada iterado (métodos iterativos), si track_iter
%     E.grid_method  procedimiento de muestreo empleado
%   T es n_valid x n_fantomas; ds(k).scenario indexa la columna de T.

if nargin < 5, keep_history = false; end
if nargin < 6, track_iter = false; end
keep_history = keep_history || track_iter;
N = numel(ds);
E.R = NaN(grid.n_valid, N);
E.RIE = NaN(1, N); E.pearson = NaN(1, N);
E.time = NaN(1, N); E.iterations = NaN(1, N); E.converged = NaN(1, N);
E.stop_reason = cell(1, N);
E.history = cell(1, N);
E.delta = cell(1, N);
E.rie_iter = cell(1, N);
E.fmdl_rec = [];
E.grid_method = '';
for k = 1:N
    [delta, fr, info] = reconstruct_one(A, ds(k).vh, ds(k).vi, keep_history);
    [img, E.grid_method] = sample_on_grid(fr, delta, grid);
    r = img(grid.mask);
    m = compute_metrics(T(:, ds(k).scenario), r);
    E.R(:, k) = r;
    E.RIE(k) = m.RIE;
    E.pearson(k) = m.pearson;
    E.time(k) = info.time;
    E.iterations(k) = info.iterations;
    E.converged(k) = info.converged;
    E.stop_reason{k} = info.stop_reason;
    E.delta{k} = delta;
    if keep_history, E.history{k} = info.history; end
    if track_iter && ~isempty(info.history)
        % Error de imagen de cada iterado, con el mismo muestreo y la misma
        % transformación (offset, escala) que la reconstrucción final.
        H = (info.history - A.offset) * A.scale;
        ri = NaN(1, size(H, 2));
        for j = 1:size(H, 2)
            im_j = sample_on_grid(fr, H(:, j), grid);
            mj = compute_metrics(T(:, ds(k).scenario), im_j(grid.mask));
            ri(j) = mj.RIE;
        end
        E.rie_iter{k} = ri;
    end
    if k == 1, E.fmdl_rec = fr; end
end
end
