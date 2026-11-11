function A = make_algorithm(name, hp, models, cfg)
%MAKE_ALGORITHM  Configura uno de los seis algoritmos con unos hiperparámetros.
%
%   A = make_algorithm(name, hp, models, cfg) devuelve una estructura con
%   todo lo necesario para reconstruir con reconstruct_one:
%     A.name, A.kind ('eidors' | 'gn_explicit'), A.imdl, A.hp
%     A.perm, A.sgn   reordenación de medidas (solo BP; [] en el resto)
%     A.scale         calibración global de amplitud (1 si no procede)
%     A.offset        valor a restar a la salida (solo inv_solve_gn, si
%                     devuelve conductividad absoluta)
%     A.time_build    tiempo de configuración del modelo (s)
%     A.time_train    tiempo de entrenamiento (GREIT; 0 en el resto)
%
%   Todos reciben los MISMOS datos (vh, vi) y reconstruyen la variación de
%   conductividad respecto al fondo homogéneo.

A = struct('name', name, 'kind', 'eidors', 'imdl', [], 'gn', [], 'hp', hp, ...
           'perm', [], 'sgn', [], 'scale', 1, 'offset', 0, ...
           'time_build', 0, 'time_train', 0, 'normalize', 0, 'term_tol', NaN);
if isfield(hp, 'scale') && ~isempty(hp.scale)
    A.scale = hp.scale;
end
base = models.inv;
t = tic;

switch upper(name)
    case 'BP'
        % Sheffield Mk I: matriz de retroproyección incluida en EIDORS.
        try
            imdl = mk_common_gridmdl('backproj');
        catch err
            error('EIT:alg:BP', ['mk_common_gridmdl(''backproj'') no está disponible en esta ', ...
                'versión de EIDORS (%s). No se sustituye por otro método.'], err.message);
        end
        imdl.reconst_type = 'difference';
        imdl = set_normalization(imdl, hp.normalize);
        A.normalize = hp.normalize;
        if isfield(hp, 'transform') && ~isempty(hp.transform) && ~isequal(hp.transform(:)', [0 1])
            [A.perm, A.sgn] = bp_orientation_permutation(models.stim, cfg.model.n_elec, ...
                hp.transform(1), hp.transform(2));
        end
        A.imdl = imdl;

    case 'NOSER'
        % Gauss-Newton linealizado de un paso con prior NOSER:
        %   dsigma = (J'J + lambda^2 * diag(J'J)^p)^-1 J' dV
        imdl = strip_priors(base);
        imdl.solve = @inv_solve_diff_GN_one_step;
        imdl.RtR_prior = @prior_noser;
        imdl.prior_noser.exponent = hp.exponent;
        imdl.hyperparameter.value = hp.lambda;
        A.imdl = imdl;

    case 'TSVD'
        % Pseudo-inversa truncada del jacobiano (sin término de Tikhonov).
        imdl = strip_priors(base);
        imdl.solve = @inv_solve_TSVD;
        imdl.hyperparameter.value = hp.threshold;
        A.imdl = imdl;

    case 'TV'
        % Variación total, método primal-dual de punto interior (Borsic).
        imdl = strip_priors(base);
        imdl.solve = @inv_solve_TV_pdipm;
        imdl.R_prior = @prior_TV;
        imdl.hyperparameter.value = hp.lambda;
        imdl.parameters.max_iterations = hp.max_iter;
        imdl.parameters.term_tolerance = cfg.tv.term_tolerance;
        imdl.parameters.keep_iterations = 1;   % permite contar iteraciones y medir convergencia
        A.term_tol = cfg.tv.term_tolerance;
        A.imdl = imdl;

    case 'GN'
        imdl = strip_priors(base);
        imdl.RtR_prior = @prior_tikhonov;
        imdl.hyperparameter.value = hp.lambda;
        switch lower(cfg.gn.implementation)
            case 'explicit'
                % Gauss-Newton iterativo sobre el modelo diferencial no lineal,
                % construido con fwd_solve, calc_jacobian y calc_RtR_prior.
                A.kind = 'gn_explicit';
                G.fmdl = base.fwd_model;
                G.sigma0 = cfg.model.sigma_bkg * ones(size(G.fmdl.elems, 1), 1);
                G.RtR = calc_RtR_prior(imdl);
                G.lambda = hp.lambda;
                G.max_iter = hp.max_iter;
                G.tol_step = cfg.gn.tol_step;
                G.tol_cost = cfg.gn.tol_cost;
                G.max_backtrack = cfg.gn.max_backtrack;
                G.min_sigma = cfg.gn.min_sigma;
                v0 = fwd_solve(mk_image(G.fmdl, G.sigma0));
                G.v0 = v0.meas(:);
                A.gn = G;
            case 'inv_solve_gn'
                imdl.solve = @inv_solve_gn;
                imdl.inv_solve_gn.max_iterations = hp.max_iter;
                imdl.inv_solve_gn.verbose = 0;
                imdl.inv_solve_gn.return_working_variables = 1;
                imdl.inv_solve_gn.elem_working = 'conductivity';
                imdl.inv_solve_gn.elem_output  = 'conductivity';
                % La salida puede ser absoluta: se calibra con una sonda de
                % datos nulos en run_pilot (A.offset).
            otherwise
                error('EIT:alg:GN', 'Implementación GN desconocida: %s', cfg.gn.implementation);
        end
        A.imdl = imdl;

    case 'GREIT'
        % Matriz de reconstrucción entrenada con objetivos simulados en la
        % malla GRUESA extruida a 3D (mk_GREIT_model no admite modelos 2D);
        % nunca en la fina que genera los datos de prueba.
        if isempty(models.greit_fwd)
            error('EIT:alg:GREIT3d', 'Falta el modelo 3D de entrenamiento (models.greit_fwd).');
        end
        fm = set_normalization(struct('fwd_model', models.greit_fwd), cfg.greit.normalize);
        opt = struct();
        opt.target_plane = cfg.greit.extrude_height / 2;
        opt.imgsz = [cfg.eval.npix cfg.eval.npix];
        if isfield(hp, 'noise_figure_used') && ~isempty(hp.noise_figure_used) && isfinite(hp.noise_figure_used)
            hp.noise_figure = hp.noise_figure_used;  % reproduce exactamente el modelo elegido en el ajuste
        end
        opt.target_size = hp.target_size;
        opt.Nsim = cfg.greit.Nsim;
        if ~isempty(cfg.greit.distr), opt.distr = cfg.greit.distr; end
        A.time_build = toc(t);
        % La búsqueda interna del peso que da el Noise Figure pedido puede
        % quedarse a ~1e-4 del objetivo y rechazarlo ("Cannot find an
        % accurate enough match"). En ese caso se reintenta con NF
        % ligeramente perturbados; el NF realmente usado queda en A.hp.
        nf0 = hp.noise_figure;
        factors = [1, cfg.greit.nf_retry_factors];
        imdl = [];
        msgs = {};
        state = rng;
        for f = factors
            opt.noise_figure = nf0 * f;
            clear_eidors_cache();                    % entrenamiento sin caché previa
            rng(cfg.greit.seed, 'twister');
            t_train = tic;
            try
                imdl = mk_GREIT_model(fm.fwd_model, cfg.greit.radius, [], opt);
                A.time_train = toc(t_train);
                break
            catch err
                msgs{end+1} = sprintf('NF=%.5g: %s', opt.noise_figure, err.message); %#ok<AGROW>
                if ~contains(err.message, 'Noise Figure', 'IgnoreCase', true)
                    break                            % otro tipo de error: no se reintenta
                end
            end
        end
        rng(state);
        if isempty(imdl)
            error('EIT:alg:GREIT', 'mk_GREIT_model ha fallado:\n  %s', strjoin(msgs, '\n  '));
        end
        if opt.noise_figure ~= nf0
            log_msg('  GREIT: NF=%.4g no convergió; se usa NF=%.5g', nf0, opt.noise_figure);
        end
        A.hp.noise_figure_used = opt.noise_figure;
        A.normalize = cfg.greit.normalize;
        A.imdl = imdl;
        return

    otherwise
        error('EIT:alg:name', 'Algoritmo desconocido: %s', name);
end
A.time_build = toc(t);
end

% ======================================================================
function imdl = strip_priors(imdl)
% Elimina el prior y la regularización por defecto de mk_common_model
% para que cada algoritmo declare explícitamente los suyos.
f = {'RtR_prior', 'R_prior', 'prior_use_fwd_not_rec'};
for i = 1:numel(f)
    if isfield(imdl, f{i}), imdl = rmfield(imdl, f{i}); end
end
if isfield(imdl, 'hyperparameter')
    imdl.hyperparameter = struct('value', []);
end
end
