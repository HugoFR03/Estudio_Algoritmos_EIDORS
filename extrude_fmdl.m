function f3 = extrude_fmdl(f2, h)
%EXTRUDE_FMDL  Extruye un fwd_model 2D a una losa 3D de una capa de espesor h.
%
%   mk_GREIT_model solo admite modelos 3D (lee la coordenada z de los
%   electrodos para fijar el plano de reconstrucción). Para entrenar GREIT
%   con una física equivalente a la del estudio 2D se extruye la MISMA
%   malla gruesa del modelo inverso:
%     - cada triángulo pasa a ser un prisma, dividido en 3 tetraedros con la
%       regla del índice global ordenado (garantiza una malla conforme);
%     - cada electrodo puntual 2D pasa a ser una línea vertical que ocupa
%       todo el espesor (nodo inferior y superior);
%     - mismas posiciones de electrodos, estimulación e impedancia de contacto.
%   Con caras superior e inferior aisladas, el potencial no depende de z y
%   las tensiones cumplen V_3D = V_2D / h. build_models verifica esta
%   equivalencia numéricamente antes de usar el modelo.

if size(f2.nodes, 2) ~= 2 || size(f2.elems, 2) ~= 3
    error('EIT:extrude:dim', 'Se esperaba un modelo 2D de triángulos.');
end
xy = f2.nodes;
n = size(xy, 1);
T = sort(double(f2.elems), 2);                 % i < j < k en cada triángulo
i = T(:, 1); j = T(:, 2); k = T(:, 3);
E = [i, j, k, k + n;
     i, j, j + n, k + n;
     i, i + n, j + n, k + n];
N = [xy, zeros(n, 1); xy, h * ones(n, 1)];

% Orientación positiva de todos los tetraedros
a = N(E(:, 2), :) - N(E(:, 1), :);
b = N(E(:, 3), :) - N(E(:, 1), :);
c = N(E(:, 4), :) - N(E(:, 1), :);
v = sum(a .* cross(b, c, 2), 2);
if any(abs(v) < eps)
    error('EIT:extrude:degenerate', 'La extrusión ha generado tetraedros degenerados.');
end
neg = v < 0;
E(neg, [3 4]) = E(neg, [4 3]);

f3 = f2;
if isfield(f2, 'name'), f3.name = [char(string(f2.name)) '_extruido']; else, f3.name = 'extruido'; end
f3.nodes = N;
f3.elems = E;
f3.boundary = find_boundary(E);
for e = 1:numel(f2.electrode)
    nd = f2.electrode(e).nodes(:)';
    f3.electrode(e).nodes = [nd, nd + n];
end
end
