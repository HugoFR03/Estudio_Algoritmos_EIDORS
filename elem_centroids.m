function c = elem_centroids(fmdl)
%ELEM_CENTROIDS  Centroides (x, y) de los elementos de un fwd_model 2D.
xy = fmdl.nodes(:, 1:2);
el = fmdl.elems;
c = zeros(size(el, 1), 2);
for j = 1:size(el, 2)
    c = c + xy(el(:, j), :);
end
c = c / size(el, 2);
end
