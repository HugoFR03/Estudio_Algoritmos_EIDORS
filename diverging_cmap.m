function c = diverging_cmap(n)
%DIVERGING_CMAP  Mapa de color divergente azul-blanco-rojo (blanco = 0).
%   Azul: disminución de conductividad; rojo: aumento.
if nargin < 1, n = 256; end
h = floor(n / 2);
blue = [0.13 0.40 0.67]; red = [0.70 0.09 0.17]; w = [1 1 1];
t = linspace(0, 1, h)';
c_neg = blue + t .* (w - blue);
c_pos = w + t .* (red - w);
c = [c_neg; c_pos(2:end, :)];
while size(c, 1) < n, c = [c; red]; end %#ok<AGROW>
end
