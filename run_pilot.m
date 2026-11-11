function [pilot, grid] = run_pilot(models, sim, ph, grid, cfg, out_dir)
%RUN_PILOT  Prueba piloto: inclusión central, sin ruido, seis algoritmos,
%   una reconstrucción por algoritmo, con hiperparámetros provisionales
%   (cfg.pilot.hp). Comprueba que:
%     1. los seis algoritmos reciben exactamente los mismos datos;
%     2. todos se ejecutan y devuelven imágenes finitas y no nulas;
%     3. las imágenes se llevan a la cuadrícula común con dimensiones
%        compatibles y la geometría (centro, escala) es coherente;
%     4. una reconstrucción de datos nulos (vi = vh) da ~0, es decir, la
%        salida es una VARIACIÓN de conductividad;
%     5. las métricas pueden calcularse y los tiempos se registran.
%   Además identifica el significado del umbral de TSVD y prueba si
%   inv_solve_gn admite reconstrucción diferencial en esta versión.
%
%   Si algún píxel de la máscara queda fuera de algún modelo de
%   reconstrucción, la máscara común se reduce (igual para todos) y se
%   informa: grid se devuelve actualizado.

s = cfg.pilot.scenario;
vh = sim.V0;
vi = sim.V1(:, s);
pilot.scenario = ph(s).name;
pilot.data_checksum = [sum(vh), sum(vi), numel(vh)];
[Tp, Timg] = truth_on_grid(models.fwd, ph(s), grid);
contrast_sign = sign(mean(Tp(Tp ~= 0)));

nA = numel(cfg.algorithms);
imgs = NaN(grid.npix, grid.npix, nA);
C = cell(1, nA);
for a = 1:nA
    name = cfg.algorithms{a};
    r = struct('algorithm', name, 'ok', false, 'message', '', 'grid_method', '', ...
        'zero_probe_max', NaN, 'offset', 0, 'n_values', NaN, 'nan_in_mask', NaN, ...
        'sign_ok', false, 'RIE', NaN, 'pearson', NaN, 'time', NaN, 'iterations', NaN, ...
        'converged', NaN, 'geometry_err_px', NaN, 'data_checksum', [NaN NaN NaN]);
    try
        A = make_algorithm(name, cfg.pilot.hp.(name), models, cfg);
        % --- Sonda de datos nulos
        d0 = reconstruct_one(A, vh, vh);
        if strcmpi(name, 'GN') && strcmpi(cfg.gn.implementation, 'inv_solve_gn') ...
                && abs(mean(d0) - cfg.model.sigma_bkg) < 0.1 * cfg.model.sigma_bkg
            A.offset = cfg.model.sigma_bkg;    % inv_solve_gn devuelve conductividad absoluta
            d0 = reconstruct_one(A, vh, vh);
        end
        r.offset = A.offset;
        r.zero_probe_max = max(abs(d0));
        % --- Reconstrucción piloto
        r.data_checksum = [sum(vh), sum(vi), numel(vh)];
        [delta, fr, info] = reconstruct_one(A, vh, vi);
        r.n_values = numel(delta);
        [img, r.grid_method] = sample_on_grid(fr, delta, grid);
        if ~isequal(size(img), [grid.npix grid.npix])
            error('EIT:pilot:size', 'Imagen de %dx%d.', size(img, 1), size(img, 2));
        end
        r.nan_in_mask = nnz(isnan(img(grid.mask)));
        r.geometry_err_px = grid_geometry_error(fr, grid);
        v = img(grid.mask);
        okv = ~isnan(v);
        m = compute_metrics(Tp(okv), v(okv));
        r.RIE = m.RIE; r.pearson = m.pearson;
        r.time = info.time; r.iterations = info.iterations; r.converged = info.converged;
        r.sign_ok = sign(mean(v(Tp ~= 0 & okv))) == contrast_sign;
        imgs(:, :, a) = img;
        r.ok = all(isfinite(delta)) && norm(delta) > 0 && isfinite(r.RIE) && r.time > 0 ...
            && r.zero_probe_max <= 1e-3 * max(abs(delta)) ...
            && r.geometry_err_px < 1.5;
        if ~r.ok
            r.message = 'Comprobación fallida (ver campos)';
        end
    catch err
        r.message = err.message;
    end
    C{a} = r;
    log_msg('  Piloto %-6s ok=%d RIE=%.3f r=%.3f t=%.4f s signo=%d sonda0=%.1e geom=%.2f px %s', ...
        name, r.ok, r.RIE, r.pearson, r.time, r.sign_ok, r.zero_probe_max, r.geometry_err_px, r.message);
end
pilot.checks = [C{:}];
cs = reshape([pilot.checks.data_checksum], 3, [])';
cs = cs(all(isfinite(cs), 2), :);
pilot.same_input = ~isempty(cs) && size(unique(cs, 'rows'), 1) == 1;

% --- Máscara común (región válida en todos los modelos)
valid = grid.mask & all(~isnan(imgs(:, :, [pilot.checks.ok])), 3);
pilot.mask_reduced_by = grid.n_valid - nnz(valid);
if pilot.mask_reduced_by > 0
    log_msg('  Máscara común reducida en %d píxeles (fuera de algún modelo de reconstrucción).', ...
        pilot.mask_reduced_by);
    grid.mask = valid;
    grid.n_valid = nnz(valid);
end

% --- Interpretación del umbral de TSVD
pilot.tsvd = [];
if cfg.tsvd.verify_interpretation && ismember('TSVD', cfg.algorithms)
    try
        pilot.tsvd = verify_tsvd_interpretation(models, cfg);
    catch err
        log_msg('  Verificación de TSVD no completada: %s', err.message);
    end
end

% --- ¿Admite inv_solve_gn reconstrucción diferencial? (informativo)
pilot.inv_solve_gn_probe = probe_inv_solve_gn(models, cfg, vh, vi);

pilot.passed = all([pilot.checks.ok]) && pilot.same_input;
pilot.images = imgs;
pilot.truth = Timg;

% --- Figura del piloto
if nargin >= 6 && ~isempty(out_dir)
    lim = max(abs(Timg(:)), [], 'omitnan');
    plot_comparison_row(cat(3, Timg(:, :, 1), imgs), [{'Real'}, cfg.alg_labels], ...
        [-lim lim], grid, sprintf('Prueba piloto: %s, sin ruido (hiperparámetros provisionales)', ...
        ph(s).label), fullfile(out_dir, 'figures', 'pilot'), cfg);
end
log_msg('Piloto %s (mismos datos: %d).', ternary(pilot.passed, 'SUPERADO', 'NO SUPERADO'), pilot.same_input);
end

% ======================================================================
function e = grid_geometry_error(fmdl, grid)
% Error medio (en píxeles) entre el radio normalizado de cada elemento y
% el radio del píxel donde cae: comprueba centro y escala del muestreo.
xy = fmdl.nodes(:, 1:2);
lo = min(xy); hi = max(xy);
ctr = (lo + hi) / 2; hw = max(hi - lo) / 2;
rad = sqrt(sum(((elem_centroids(fmdl) - ctr) / hw).^2, 2));
img = sample_on_grid(fmdl, rad, grid);
d = abs(img - sqrt(grid.X.^2 + grid.Y.^2));
d = d(grid.mask & ~isnan(img));
e = mean(d) / grid.pixel_size;
end

function P = probe_inv_solve_gn(models, cfg, vh, vi)
P = struct('available', exist('inv_solve_gn', 'file') > 0, 'runs', false, ...
    'output_is_absolute', NaN, 'iterations_readable', false, 'message', '');
if ~P.available, P.message = 'inv_solve_gn no disponible'; return, end
c2 = cfg; c2.gn.implementation = 'inv_solve_gn';
try
    A = make_algorithm('GN', cfg.pilot.hp.GN, models, c2);
    img0 = inv_solve(A.imdl, vh, vh);
    img1 = inv_solve(A.imdl, vh, vi);
    P.runs = all(isfinite(img1.elem_data(:)));
    P.output_is_absolute = abs(mean(img0.elem_data(:, end)) - cfg.model.sigma_bkg) < 0.1;
    [~, ~, info] = reconstruct_one(A, vh, vi);
    P.iterations_readable = ~isnan(info.iterations);
    P.message = 'ejecución correcta';
catch err
    P.message = err.message;
end
log_msg('  Sonda inv_solve_gn (diferencial): ejecuta=%d, salida absoluta=%d, iteraciones legibles=%d (%s)', ...
    P.runs, P.output_is_absolute, P.iterations_readable, P.message);
end

function out = ternary(c, a, b)
if c, out = a; else, out = b; end
end
