function [delta, info] = gn_tikhonov_diff(G, vh, vi, keep_history)
%GN_TIKHONOV_DIFF  Gauss-Newton iterativo con regularización de Tikhonov
%   para EIT diferencial, construido con las rutinas de EIDORS
%   (fwd_solve, calc_jacobian, calc_RtR_prior).
%
%   Problema (malla gruesa, sigma0 = fondo homogéneo):
%     min_d  Phi(d) = || dV - [F(sigma0 + d) - F(sigma0)] ||^2 + lambda^2 * d' R d
%   con dV = vi - vh y R = prior_tikhonov.
%
%   Iteración k (d_0 = 0):
%     J_k = dF/dsigma en sigma0 + d_k                 (calc_jacobian)
%     r_k = dV - [F(sigma0 + d_k) - F(sigma0)]        (fwd_solve)
%     (J_k'J_k + lambda^2 R) p_k = J_k' r_k - lambda^2 R d_k
%     d_{k+1} = d_k + a_k p_k,  a_k = 1, 1/2, 1/4, ... hasta que Phi descienda
%                                (y sigma0 + d se mantenga >= min_sigma)
%   Criterios de parada (por orden):
%     'step_tolerance' : ||a_k p_k|| / ||d_{k+1}|| < tol_step
%     'cost_tolerance' : (Phi_k - Phi_{k+1}) / Phi_k < tol_cost
%     'line_search'    : ningún paso reduce Phi (no convergido)
%     'max_iterations' : se alcanza max_iter (no convergido)
%
%   La primera iteración coincide con la reconstrucción Tikhonov de un paso;
%   las siguientes corrigen la no linealidad del problema directo.

if nargin < 4, keep_history = false; end
t = tic;
dV = vi(:) - vh(:);
n = numel(G.sigma0);
delta = zeros(n, 1);
lam2 = G.lambda^2;
img = mk_image(G.fmdl, G.sigma0);
r = dV;                         % residuo en d = 0
cost = r' * r;
info.cost = cost;
info.alpha = [];
info.rel_step = [];
info.history = [];
if keep_history, info.history = zeros(n, 0); end
reason = 'max_iterations';
converged = false;
k = 0;
while k < G.max_iter
    k = k + 1;
    img.elem_data = G.sigma0 + delta;
    J = calc_jacobian(img);
    if size(J, 1) ~= numel(dV) || size(J, 2) ~= n
        error('EIT:gn:dim', 'Jacobiano de %dx%d incompatible con %d medidas y %d elementos.', ...
            size(J, 1), size(J, 2), numel(dV), n);
    end
    H = J' * J + lam2 * G.RtR;
    g = J' * r - lam2 * (G.RtR * delta);
    p = H \ g;
    a = 1; accepted = false;
    for b = 0:G.max_backtrack
        cand = delta + a * p;
        if min(G.sigma0 + cand) > G.min_sigma
            img.elem_data = G.sigma0 + cand;
            vc = fwd_solve(img);
            rc = dV - (vc.meas(:) - G.v0);
            cost_c = rc' * rc + lam2 * (cand' * (G.RtR * cand));
            if cost_c <= cost
                accepted = true; break
            end
        end
        a = a / 2;
    end
    if ~accepted
        k = k - 1;                      % la iteración no se completa
        reason = 'line_search';
        break
    end
    rel_step = norm(a * p) / max(norm(cand), eps);
    rel_cost = (cost - cost_c) / max(cost, eps);
    delta = cand; r = rc; cost = cost_c;
    info.cost(end+1) = cost;
    info.alpha(end+1) = a;
    info.rel_step(end+1) = rel_step;
    if keep_history, info.history(:, end+1) = delta; end
    if rel_step < G.tol_step
        reason = 'step_tolerance'; converged = true; break
    end
    if rel_cost < G.tol_cost
        reason = 'cost_tolerance'; converged = true; break
    end
end
info.iterations = k;
info.converged = converged;
info.stop_reason = reason;
if isempty(info.rel_step), info.rel_change_last = NaN; else, info.rel_change_last = info.rel_step(end); end
info.time = toc(t);
end
