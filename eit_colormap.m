function c = eit_colormap(name, n)
%EIT_COLORMAP  Mapas de color para imágenes de variación de conductividad.
%
%   c = eit_colormap('eit', n)   (por defecto)
%       Pensado para límites simétricos [-L, L]: el centro (dsigma = 0, el
%       fondo) es azul; los aumentos de conductividad avanzan hacia colores
%       cálidos (cian, verde, amarillo, rojo) y las disminuciones hacia azul
%       muy oscuro. Da más contraste que un fondo blanco al proyectar las
%       imágenes y sigue distinguiendo el signo del cambio.
%   c = eit_colormap('diverging', n)
%       Azul-blanco-rojo, blanco = sin cambio.

if nargin < 1 || isempty(name), name = 'eit'; end
if nargin < 2, n = 256; end
switch lower(name)
    case 'eit'
        h = floor(n / 2);
        % Mitad inferior: -L (casi negro azulado) -> 0 (azul de fondo)
        neg = interp_stops([0 0.5 1], [0.02 0.02 0.10; 0.05 0.09 0.42; 0.12 0.26 0.78], h);
        % Mitad superior: 0 (azul) -> cian -> verde -> amarillo -> rojo oscuro (+L)
        pos = interp_stops([0 0.22 0.45 0.70 0.88 1], ...
            [0.12 0.26 0.78; 0.00 0.62 0.95; 0.35 0.82 0.45; 1.00 0.85 0.10; ...
             0.95 0.40 0.05; 0.70 0.05 0.05], n - h + 1);
        c = [neg; pos(2:end, :)];
    case 'diverging'
        c = diverging_cmap(n);
    otherwise
        error('EIT:cmap:name', 'Mapa de color desconocido: %s', name);
end
end

function c = interp_stops(t, rgb, n)
x = linspace(0, 1, n)';
c = [interp1(t, rgb(:, 1), x), interp1(t, rgb(:, 2), x), interp1(t, rgb(:, 3), x)];
c = min(max(c, 0), 1);
end
