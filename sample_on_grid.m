function [img, method] = sample_on_grid(fmdl, elem_data, grid)
%SAMPLE_ON_GRID  Lleva una imagen por elementos a la cuadrícula común.
%
%   [img, method] = sample_on_grid(fmdl, elem_data, grid) devuelve una
%   matriz npix x npix con el valor del elemento que contiene el centro de
%   cada píxel (interpolación constante a trozos), y NaN fuera de la malla.
%   Es el MISMO procedimiento para la distribución real (malla fina) y
%   para todas las reconstrucciones (malla gruesa, rejilla de Sheffield o
%   rejilla de GREIT).
%
%   Las coordenadas de cada modelo se normalizan con su caja envolvente
%   (centro y semilado), de modo que el dominio circular ocupa [-1, 1]^2
%   en todos ellos, cualquiera que sea la unidad interna del modelo.
%
%   Si el número de valores no coincide con el de elementos y el modelo
%   tiene coarse2fine (modelos de reconstrucción por píxeles), se aplica
%   esa correspondencia. Como último recurso se usa calc_slices de EIDORS.

ed = elem_data(:, end);
nE = size(fmdl.elems, 1);
method = 'pointLocation';
if numel(ed) ~= nE
    if isfield(fmdl, 'coarse2fine') && size(fmdl.coarse2fine, 2) == numel(ed) ...
            && size(fmdl.coarse2fine, 1) == nE
        ed = fmdl.coarse2fine * ed;
        method = 'coarse2fine+pointLocation';
    else
        [img, method] = fallback_slices(fmdl, ed, grid);
        return
    end
end
try
    xy = fmdl.nodes(:, 1:2);
    lo = min(xy); hi = max(xy);
    ctr = (lo + hi) / 2; hw = max(hi - lo) / 2;
    xq = ctr(1) + grid.X(:) * hw;
    yq = ctr(2) + grid.Y(:) * hw;
    if size(fmdl.elems, 2) ~= 3
        error('EIT:grid:elem', 'Elementos no triangulares.');
    end
    % Algunos modelos de EIDORS tienen nodos que no pertenecen a ningún
    % elemento; no afecta al muestreo, así que se silencia ese aviso.
    ws = warning('off', 'MATLAB:triangulation:PtsNotInTriWarnId');
    TR = triangulation(double(fmdl.elems), xy);
    warning(ws);
    id = pointLocation(TR, xq, yq);
    v = NaN(numel(xq), 1);
    ok = ~isnan(id);
    v(ok) = ed(id(ok));
    img = reshape(v, grid.npix, grid.npix);
catch err
    warning('EIT:grid:fallback', 'Muestreo directo fallido (%s); se usa calc_slices.', err.message);
    [img, method] = fallback_slices(fmdl, ed, grid);
end
end

function [img, method] = fallback_slices(fmdl, ed, grid)
im = mk_image(fmdl, 0);
im.elem_data = ed;
im.calc_colours.npoints = grid.npix;
s = calc_slices(im);
if ~isequal(size(s), [grid.npix grid.npix])
    error('EIT:grid:slices', 'calc_slices devolvió %dx%d en lugar de %dx%d.', ...
        size(s, 1), size(s, 2), grid.npix, grid.npix);
end
img = flipud(s);   % calc_slices: primera fila = y máxima
method = 'calc_slices';
end
