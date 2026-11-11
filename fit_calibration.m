function [s, rie] = fit_calibration(T, R, scen)
%FIT_CALIBRATION  Escalar global de amplitud que minimiza el RIE medio.
%
%   [s, rie] = fit_calibration(T, R, scen)
%     T     n_valid x n_fantomas, distribuciones reales
%     R     n_valid x N, reconstrucciones sin calibrar
%     scen  1 x N, columna de T asociada a cada reconstrucción
%   Minimiza mean_k ||t_k - s r_k|| / ||t_k||, función convexa en s, con
%   fminbnd. Se aplica solo a métodos cuya salida no está en unidades
%   físicas (BP de Sheffield, GREIT); es un único valor por algoritmo,
%   fijado en validación, no una normalización imagen a imagen. Admite
%   s < 0 (por ejemplo, si la salida es un cambio de resistividad).

ok = all(isfinite(R), 1);
if ~any(ok)
    s = NaN; rie = Inf; return
end
T = T(:, scen(ok)); R = R(:, ok);
w = 1 ./ sum(T.^2, 1);
s0 = sum(w .* sum(T .* R, 1)) / max(sum(w .* sum(R.^2, 1)), eps);
f = @(x) mean(sqrt(sum((T - x * R).^2, 1)) ./ sqrt(sum(T.^2, 1)));
a = 5 * abs(s0) + 1e-12;
[s, rie] = fminbnd(f, -a, a, optimset('TolX', 1e-6 * a));
if f(s0) < rie
    s = s0; rie = f(s0);
end
end
