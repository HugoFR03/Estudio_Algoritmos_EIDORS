function m = compute_metrics(t, r)
%COMPUTE_METRICS  Métricas comunes sobre los píxeles válidos del dominio.
%
%   m = compute_metrics(t, r)
%     t   variación real de conductividad (vector de píxeles válidos)
%     r   variación reconstruida (mismos píxeles, mismo orden)
%   Sin normalización individual de ninguna de las dos imágenes.
%
%   m.RIE      = ||t - r||_2 / ||t||_2           (0 = reconstrucción exacta)
%   m.pearson  = correlación de Pearson entre t y r, en [-1, 1]
%   m.pearson_valid = false si r o t no tienen variabilidad o r no es
%                finita (en ese caso m.pearson = NaN)

t = t(:); r = r(:);
if numel(t) ~= numel(r)
    error('EIT:metrics:dim', 'Imágenes de distinto tamaño (%d y %d píxeles).', numel(t), numel(r));
end
if norm(t) == 0
    error('EIT:metrics:zero', 'La distribución real no tiene variación: el RIE no está definido.');
end
m.RIE = NaN;
m.pearson = NaN;
m.pearson_valid = false;
if any(~isfinite(r))
    return
end
m.RIE = norm(t - r) / norm(t);
if std(t) > 0 && std(r) > 0
    C = corrcoef(t, r);
    m.pearson = C(1, 2);
    m.pearson_valid = true;
end
end
