function algs = build_inverse_models(models, hp_final, cfg, pilot)
%BUILD_INVERSE_MODELS  Configura los seis algoritmos con los hiperparámetros
%   definitivos (fijados en validación) y mide su tiempo de preparación.
%
%   algs{a} es la estructura de make_algorithm, con además:
%     .label            nombre para tablas y figuras
%     .time_setup_total tiempo total de configuración (incluye entrenamiento)
%   El entrenamiento de GREIT se repite aquí con la misma semilla que en el
%   ajuste, de modo que su tiempo se mide sin caché.

algs = cell(1, numel(cfg.algorithms));
for a = 1:numel(cfg.algorithms)
    name = cfg.algorithms{a};
    clear_eidors_cache();
    t = tic;
    A = make_algorithm(name, hp_final.(name), models, cfg);
    A.time_setup_total = toc(t);
    if nargin >= 4 && ~isempty(pilot)
        chk = pilot.checks(strcmp({pilot.checks.algorithm}, name));
        if ~isempty(chk), A.offset = chk.offset; end
    end
    A.label = cfg.alg_labels{a};
    % Construcción del modelo FEM del que parte cada método (común, se mide
    % una vez en build_models): BP usa su rejilla precalculada (incluida en
    % time_build); GREIT necesita además la extrusión a 3D.
    switch name
        case 'BP',    A.time_model = 0;
        case 'GREIT', A.time_model = models.time.inv_build + models.time.greit_build;
        otherwise,    A.time_model = models.time.inv_build;
    end
    algs{a} = A;
    log_msg('  %-6s preparado en %.3f s (entrenamiento %.3f s)', name, A.time_setup_total, A.time_train);
end
end
