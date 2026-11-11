function V = verify_tsvd_interpretation(models, cfg, thresholds)
%VERIFY_TSVD_INTERPRETATION  Identifica qué significa el hiperparámetro de
%   inv_solve_TSVD en la versión instalada de EIDORS.
%
%   Para cada umbral se obtiene la matriz de reconstrucción efectiva de
%   EIDORS (RM_e) reconstruyendo las columnas de la identidad, y se compara
%   con la pseudo-inversa truncada explícita del mismo jacobiano,
%       RM_k = V_k * diag(1./s_k) * U_k',   k = 1..rank(J).
%   El k que mejor la reproduce indica cuántos valores singulares conserva
%   EIDORS; se contrastan cuatro hipótesis:
%       'relativo %': s_i >= (umbral/100) * s_1   (documentada en EIDORS 3.12)
%       'relativo'  : s_i >= umbral * s_1
%       'absoluto'  : s_i >= umbral
%       'energía'   : s_i^2 >= umbral * s_1^2
%   Si ningún k reproduce RM_e (error relativo > 1e-3), EIDORS no aplica
%   una TSVD pura (por ejemplo, incluye ponderaciones) y se avisa.

if nargin < 3, thresholds = cfg.hp.TSVD.threshold([1 round(end/2) end]); end
base = make_algorithm('TSVD', struct('threshold', thresholds(1)), models, cfg);
fmdl = base.imdl.fwd_model;
nm = models.n_meas;
J = calc_jacobian(mk_image(fmdl, cfg.model.sigma_bkg));
[U, S, W] = svd(full(J), 'econ');
s = diag(S);
r = nnz(s > max(size(J)) * eps(s(1)));
V.singular_values = s;
V.rank = r;
V.rows = struct('threshold', {}, 'k_best', {}, 'rel_error', {}, 'k_percent', {}, ...
    'k_relative', {}, 'k_absolute', {}, 'k_energy', {}, 'match', {});
for i = 1:numel(thresholds)
    A = make_algorithm('TSVD', struct('threshold', thresholds(i)), models, cfg);
    try
        img = inv_solve(A.imdl, zeros(nm, 1), eye(nm));   % un fotograma de referencia
    catch
        img = inv_solve(A.imdl, zeros(nm, nm), eye(nm));  % referencia por fotograma
    end
    RMe = img.elem_data;
    if ~isequal(size(RMe), [size(J, 2) nm])
        error('EIT:tsvd:dim', 'La matriz efectiva de TSVD tiene %dx%d.', size(RMe, 1), size(RMe, 2));
    end
    nrm = norm(RMe, 'fro');
    err = zeros(r, 1);
    RMk = zeros(size(RMe));
    for k = 1:r
        RMk = RMk + W(:, k) * (U(:, k)' / s(k));
        err(k) = norm(RMe - RMk, 'fro') / max(nrm, eps);
    end
    [e, kb] = min(err);
    kp = nnz(s(1:r) >= thresholds(i) / 100 * s(1));
    kr = nnz(s(1:r) >= thresholds(i) * s(1));
    ka = nnz(s(1:r) >= thresholds(i));
    ke = nnz(s(1:r).^2 >= thresholds(i) * s(1)^2);
    m = {};
    if kb == kp, m{end+1} = 'relativo %'; end %#ok<AGROW>
    if kb == kr, m{end+1} = 'relativo'; end %#ok<AGROW>
    if kb == ka, m{end+1} = 'absoluto'; end %#ok<AGROW>
    if kb == ke, m{end+1} = 'energía'; end %#ok<AGROW>
    V.rows(i) = struct('threshold', thresholds(i), 'k_best', kb, 'rel_error', e, 'k_percent', kp, ...
        'k_relative', kr, 'k_absolute', ka, 'k_energy', ke, 'match', strjoin(m, '/'));
    log_msg('  TSVD umbral %.3g: EIDORS conserva k=%d valores singulares (error %.2e); hipótesis: %s', ...
        thresholds(i), kb, e, ifempty(strjoin(m, '/'), 'ninguna'));
end
V.pure_tsvd = all([V.rows.rel_error] < 1e-3);
if ~V.pure_tsvd
    warning('EIT:tsvd:interp', ['inv_solve_TSVD no coincide con una pseudo-inversa truncada ', ...
        'del jacobiano (error > 1e-3). Revise su implementación antes de interpretar el umbral.']);
end
end

function s = ifempty(s, d)
if isempty(s), s = d; end
end
