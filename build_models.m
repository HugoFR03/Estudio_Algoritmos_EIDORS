function models = build_models(cfg)
%BUILD_MODELS  Construye los modelos FEM directo (fino) e inverso (grueso).
%
%   models = build_models(cfg) devuelve:
%     models.fwd      fwd_model fino, usado SOLO para simular las medidas
%     models.inv      inv_model base (malla gruesa) del que derivan NOSER,
%                     TSVD, TV, GN y GREIT
%     models.stim     patrón de estimulación/medida común (adyacente)
%     models.n_meas   número de medidas por fotograma
%     models.R        radio del dominio
%     models.time     tiempos de construcción (s)
%
%   Ambos modelos comparten geometría nominal, electrodos, estimulación y
%   protocolo de medida; difieren en la discretización (evita el inverse
%   crime). Se verifican radio, posición angular de los electrodos y número
%   de medidas.

n = cfg.model.n_elec;
t_all = tic;

% --- Patrón de estimulación común
stim = mk_stim_patterns(n, 1, cfg.model.stim.inj, cfg.model.stim.meas, ...
                        cfg.model.stim.options, cfg.model.stim.amplitude);

% --- Modelo directo (fino)
t = tic;
try
    imdl_f = mk_common_model(cfg.model.fwd_mesh, n);
catch err
    error('EIT:models:fwd', ['mk_common_model(''%s'',%d) ha fallado (%s). ', ...
        'Revise cfg.model.fwd_mesh en config_experiment.m.'], cfg.model.fwd_mesh, n, err.message);
end
fmdl = prepare_fmdl(imdl_f.fwd_model, stim);
models.time.fwd_build = toc(t);

% --- Modelo inverso (grueso)
t = tic;
try
    imdl_i = mk_common_model(cfg.model.inv_mesh, n);
catch err
    error('EIT:models:inv', ['mk_common_model(''%s'',%d) ha fallado (%s). ', ...
        'Revise cfg.model.inv_mesh en config_experiment.m.'], cfg.model.inv_mesh, n, err.message);
end
imdl_i.fwd_model = prepare_fmdl(imdl_i.fwd_model, stim);
imdl_i.reconst_type = 'difference';
imdl_i.jacobian_bkgnd.value = cfg.model.sigma_bkg;
imdl_i = set_normalization(imdl_i, 0);
models.time.inv_build = toc(t);

% --- Verificaciones geométricas
Rf = domain_radius(fmdl);
Ri = domain_radius(imdl_i.fwd_model);
if abs(Rf - Ri) / Rf > cfg.model.radius_tol
    error('EIT:models:radius', 'Radios distintos: directo %.4f, inverso %.4f.', Rf, Ri);
end
if numel(fmdl.electrode) ~= n || numel(imdl_i.fwd_model.electrode) ~= n
    error('EIT:models:nelec', 'Los modelos no tienen %d electrodos.', n);
end
af = electrode_angles(fmdl);
ai = electrode_angles(imdl_i.fwd_model);
dang = abs(angle(exp(1i * (af - ai)))) * 180 / pi;
if max(dang) > cfg.model.elec_angle_tol_deg
    error('EIT:models:elec', ['La posición angular de los electrodos difiere hasta %.2f grados ', ...
        'entre el modelo directo y el inverso.'], max(dang));
end

% --- Número de medidas (resolviendo el caso homogéneo en ambos modelos)
vf = fwd_solve(mk_image(fmdl, cfg.model.sigma_bkg));
vi = fwd_solve(mk_image(imdl_i.fwd_model, cfg.model.sigma_bkg));
if numel(vf.meas) ~= numel(vi.meas)
    error('EIT:models:nmeas', 'Número de medidas distinto: directo %d, inverso %d.', ...
        numel(vf.meas), numel(vi.meas));
end
expected = n * (n - 3);
if any(strcmp(cfg.model.stim.options, 'no_meas_current')) && numel(vf.meas) ~= expected
    warning('EIT:models:nmeasExp', 'Se esperaban %d medidas y hay %d.', expected, numel(vf.meas));
end

% --- Modelo 3D equivalente para el entrenamiento de GREIT
models.greit_fwd = [];
models.greit_equiv = [];
if ismember('GREIT', cfg.algorithms)
    [models.greit_fwd, models.greit_equiv] = build_greit_model(imdl_i.fwd_model, cfg);
    models.time.greit_build = models.greit_equiv.time_extrude;
else
    models.time.greit_build = 0;
end

models.fwd    = fmdl;
models.inv    = imdl_i;
models.stim   = stim;
models.n_meas = numel(vf.meas);
models.R      = Rf;
models.mesh.fwd = struct('name', cfg.model.fwd_mesh, 'n_elems', size(fmdl.elems, 1), ...
                         'n_nodes', size(fmdl.nodes, 1));
models.mesh.inv = struct('name', cfg.model.inv_mesh, 'n_elems', size(imdl_i.fwd_model.elems, 1), ...
                         'n_nodes', size(imdl_i.fwd_model.nodes, 1));
models.elec_angle_max_diff_deg = max(dang);
models.time.total = toc(t_all);

log_msg('Modelo directo %s: %d elementos | modelo inverso %s: %d elementos | %d medidas', ...
    cfg.model.fwd_mesh, models.mesh.fwd.n_elems, cfg.model.inv_mesh, ...
    models.mesh.inv.n_elems, models.n_meas);
log_msg('Radio %.4f; desfase angular máximo entre electrodos %.3f grados', Rf, max(dang));
end

% ======================================================================
function [f3, eq] = build_greit_model(f2, cfg)
% Extrusión de la malla gruesa y verificación de su equivalencia con el 2D:
%   err_homog = ||h*V3 - V2|| / ||V2||          (medio homogéneo)
%   err_diff  = ||h*dV3 - dV2|| / ||dV2||       (inclusión central de prueba)
h = cfg.greit.extrude_height;
t_ex = tic;
f3 = extrude_fmdl(f2, h);
eq.time_extrude = toc(t_ex);           % solo la extrusión; la verificación no cuenta
s0 = cfg.model.sigma_bkg;
V2 = fwd_solve(mk_image(f2, s0));  V2 = V2.meas(:);
try
    V3 = fwd_solve(mk_image(f3, s0));  V3 = V3.meas(:);
catch err
    error('EIT:models:greit3dSolve', ['EIDORS no ha podido resolver el modelo 3D extruido ', ...
        '(electrodos lineales de 2 nodos): %s'], err.message);
end
if numel(V3) ~= numel(V2)
    error('EIT:models:greit3dMeas', 'El modelo 3D da %d medidas y el 2D %d.', numel(V3), numel(V2));
end
c2 = elem_centroids(f2);
in2 = sum(c2.^2, 2) <= 0.2^2;
s2 = s0 * ones(size(f2.elems, 1), 1); s2(in2) = 2 * s0;
s3 = repmat(s2, 3, 1);                % los 3 tetraedros de cada prisma, en el orden de extrude_fmdl
dV2 = fwd_solve(mk_image(f2, s2)); dV2 = dV2.meas(:) - V2;
dV3 = fwd_solve(mk_image(f3, s3)); dV3 = dV3.meas(:) - V3;
eq.height = h;
eq.err_homog = norm(h * V3 - V2) / norm(V2);
eq.err_diff  = norm(h * dV3 - dV2) / norm(dV2);
eq.n_elems = size(f3.elems, 1);
log_msg('Modelo 3D para GREIT: losa de %.3g con %d tetraedros; equivalencia con 2D: %.2e (homogéneo), %.2e (diferencial)', ...
    h, eq.n_elems, eq.err_homog, eq.err_diff);
if eq.err_homog > cfg.greit.equiv_tol || eq.err_diff > cfg.greit.equiv_tol
    error('EIT:models:greit3d', ['El modelo 3D extruido no reproduce el 2D (errores %.3g y %.3g, ', ...
        'tolerancia %.3g). Revise el tratamiento de electrodos puntuales en 3D.'], ...
        eq.err_homog, eq.err_diff, cfg.greit.equiv_tol);
end
end

function fmdl = prepare_fmdl(fmdl, stim)
% Impone el patrón común. Con 'no_meas_current' el patrón ya excluye los
% electrodos de corriente, así que se elimina cualquier meas_select previo.
fmdl.stimulation = stim;
if isfield(fmdl, 'meas_select')
    fmdl = rmfield(fmdl, 'meas_select');
end
end

function R = domain_radius(fmdl)
xy = fmdl.nodes(:, 1:2);
c  = (max(xy) + min(xy)) / 2;
R  = max(sqrt(sum((xy - c).^2, 2)));
end

function a = electrode_angles(fmdl)
xy = fmdl.nodes(:, 1:2);
c  = (max(xy) + min(xy)) / 2;
a  = zeros(numel(fmdl.electrode), 1);
for k = 1:numel(fmdl.electrode)
    p = mean(xy(fmdl.electrode(k).nodes, :), 1) - c;
    a(k) = atan2(p(2), p(1));
end
end
