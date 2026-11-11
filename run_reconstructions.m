function res = run_reconstructions(algs, ds, T, grid, cfg)
%RUN_RECONSTRUCTIONS  Evaluación definitiva: todos los conjuntos de prueba
%   con los seis algoritmos y los hiperparámetros ya fijados.
%
%   Para cada algoritmo:
%     1. se vacía la caché de EIDORS y se cronometra la primera llamada
%        (incluye la precomputación de matrices y jacobianos);
%     2. se hacen cfg.timing.n_warmup reconstrucciones de calentamiento,
%        no registradas;
%     3. se reconstruyen los N conjuntos, cronometrando solo la llamada al
%        solucionador.
%   Precomputación = max(0, primera llamada - mediana de las siguientes).

N = numel(ds);
nA = numel(algs);
C = cell(1, nA);
for a = 1:nA
    A = algs{a};
    log_msg('Reconstruyendo %d conjuntos con %s ...', N, A.name);
    clear_eidors_cache();
    t = tic;
    reconstruct_one(A, ds(1).vh, ds(1).vi);
    t_first = toc(t);
    for w = 1:cfg.timing.n_warmup
        reconstruct_one(A, ds(1).vh, ds(1).vi);
    end
    % En los métodos iterativos se registra además el RIE de cada iterado
    % (fuera del intervalo cronometrado de cada reconstrucción).
    track = ismember(A.name, cfg.iter_tracking.algorithms);
    E = evaluate_sets(A, ds, T, grid, false, track);
    if ~isfield(A, 'time_model'), A.time_model = 0; end
    imgs = NaN(grid.npix, grid.npix, N);
    for k = 1:N
        im = NaN(grid.npix);
        im(grid.mask) = E.R(:, k);
        imgs(:, :, k) = im;
    end
    C{a} = struct('name', A.name, 'label', A.label, 'hp', A.hp, ...
        'RIE', E.RIE, 'pearson', E.pearson, 'time', E.time, ...
        'iterations', E.iterations, 'converged', E.converged, ...
        'stop_reason', {E.stop_reason}, 'images', imgs, 'delta', {E.delta}, ...
        'rie_iter', {E.rie_iter}, 'time_model', A.time_model, ...
        'grid_method', E.grid_method, ...
        'time_build', A.time_build, 'time_train', A.time_train, ...
        'time_setup_total', A.time_setup_total, 'time_first_call', t_first, ...
        'time_precompute', max(0, t_first - median(E.time, 'omitnan')), ...
        'scale', A.scale, 'offset', A.offset, 'normalize', A.normalize);
    log_msg('  %s: RIE medio %.4f, Pearson medio %.4f, t mediano %.4f s', A.name, ...
        mean(E.RIE, 'omitnan'), mean(E.pearson, 'omitnan'), median(E.time, 'omitnan'));
end
res = [C{:}];
end
