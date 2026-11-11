function [delta, fmdl_rec, info] = reconstruct_one(A, vh, vi, keep_history)
%RECONSTRUCT_ONE  Una reconstrucción diferencial con el algoritmo A.
%
%   [delta, fmdl_rec, info] = reconstruct_one(A, vh, vi)
%     delta     variación de conductividad reconstruida (por elemento o
%               por píxel, según el modelo de reconstrucción)
%     fmdl_rec  modelo sobre el que está definida delta
%     info      .time (s, solo la llamada al solucionador), .iterations,
%               .converged, .stop_reason, .rel_change_last
%
%   La reordenación de medidas de BP (si existe) se aplica por igual a vh y
%   vi; es una renumeración de electrodos, no un preprocesado de los datos.

if nargin < 4, keep_history = false; end
vh = vh(:); vi = vi(:);
if ~isempty(A.perm)
    vh = A.sgn .* vh(A.perm);
    vi = A.sgn .* vi(A.perm);
end
info = struct('time', NaN, 'iterations', NaN, 'converged', NaN, ...
              'stop_reason', 'n/a', 'rel_change_last', NaN, 'history', []);

switch A.kind
    case 'gn_explicit'
        [delta, gi] = gn_tikhonov_diff(A.gn, vh, vi, keep_history);
        fmdl_rec = A.gn.fmdl;
        info.time = gi.time;
        info.iterations = gi.iterations;
        info.converged = gi.converged;
        info.stop_reason = gi.stop_reason;
        info.rel_change_last = gi.rel_change_last;
        info.history = gi.history;

    case 'eidors'
        t = tic;
        img = inv_solve(A.imdl, vh, vi);
        info.time = toc(t);
        ed = img.elem_data;
        fmdl_rec = img.fwd_model;
        if strcmpi(A.name, 'TV')
            % keep_iterations = 1: una columna por iteración
            info.iterations = size(ed, 2);
            if size(ed, 2) > 1
                info.rel_change_last = norm(ed(:, end) - ed(:, end-1)) / max(norm(ed(:, end)), eps);
                info.converged = size(ed, 2) < A.hp.max_iter ...
                    || info.rel_change_last < A.term_tol;
                if size(ed, 2) < A.hp.max_iter
                    info.stop_reason = 'tolerance';
                else
                    info.stop_reason = 'max_iterations';
                end
            else
                info.stop_reason = 'unknown (keep_iterations no soportado)';
            end
            if keep_history, info.history = ed; end
        elseif strcmpi(A.name, 'GN')
            info = gn_core_info(img, info, A);
        end
        delta = ed(:, end);

    otherwise
        error('EIT:rec:kind', 'Tipo de algoritmo desconocido: %s', A.kind);
end
delta = (delta - A.offset) * A.scale;
if any(~isfinite(delta))
    warning('EIT:rec:nan', '%s ha devuelto valores no finitos.', A.name);
end
end

function info = gn_core_info(img, info, A)
% Extracción defensiva del número de iteraciones de inv_solve_gn: los
% nombres de los campos dependen de la versión de EIDORS.
info.iterations = NaN;
if isfield(img, 'info')
    f = {'iters', 'iter', 'iterations', 'n_iter'};
    for i = 1:numel(f)
        if isfield(img.info, f{i}), info.iterations = img.info.(f{i})(end); break, end
    end
    if isnan(info.iterations) && isfield(img.info, 'res_data')
        info.iterations = numel(img.info.res_data) - 1;
    end
end
if ~isnan(info.iterations)
    info.converged = info.iterations < A.hp.max_iter;
    info.stop_reason = 'inv_solve_gn';
end
end
