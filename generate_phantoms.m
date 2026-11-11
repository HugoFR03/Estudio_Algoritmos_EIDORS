function ph = generate_phantoms(fmdl, scenarios, sigma_bkg)
%GENERATE_PHANTOMS  Distribuciones de conductividad sobre la malla directa.
%
%   ph = generate_phantoms(fmdl, scenarios, sigma_bkg) devuelve, para cada
%   escenario, la conductividad por elemento (ph(s).sigma), su variación
%   respecto al fondo (ph(s).dsigma) y la descripción geométrica.
%
%   Un elemento pertenece a una inclusión si su centroide cae dentro del
%   círculo. Las coordenadas de los escenarios están normalizadas al radio
%   del dominio. Si dos inclusiones se solapan prevalece la última.

validate_scenarios(scenarios);
xy = fmdl.nodes(:, 1:2);
ctr = (max(xy) + min(xy)) / 2;
R = max(sqrt(sum((xy - ctr).^2, 2)));
cen = elem_centroids(fmdl) - ctr;
nE = size(fmdl.elems, 1);

C = cell(1, numel(scenarios));
for s = 1:numel(scenarios)
    sig = sigma_bkg * ones(nE, 1);
    inc = scenarios(s).inclusions;
    n_in = zeros(1, numel(inc));
    for k = 1:numel(inc)
        c = inc(k).center(:)' * R;
        r = inc(k).radius * R;
        inside = sum((cen - c).^2, 2) <= r^2;
        n_in(k) = nnz(inside);
        if n_in(k) == 0
            error('EIT:phantom:empty', ...
                'La inclusión %d del escenario %s no contiene ningún elemento de la malla.', ...
                k, scenarios(s).name);
        end
        sig(inside) = inc(k).sigma;
    end
    C{s} = struct('name', scenarios(s).name, 'label', scenarios(s).label, ...
        'inclusions', inc, 'sigma', sig, 'dsigma', sig - sigma_bkg, ...
        'n_elems_inclusion', n_in);
end
ph = [C{:}];
end

% ======================================================================
function validate_scenarios(sc)
for s = 1:numel(sc)
    inc = sc(s).inclusions;
    if isempty(inc)
        error('EIT:phantom:noinc', 'El escenario %s no tiene inclusiones.', sc(s).name);
    end
    for k = 1:numel(inc)
        c = inc(k).center;
        if numel(c) ~= 2 || inc(k).radius <= 0 || inc(k).sigma <= 0
            error('EIT:phantom:def', 'Inclusión %d de %s mal definida.', k, sc(s).name);
        end
        if norm(c) + inc(k).radius > 1
            error('EIT:phantom:out', 'La inclusión %d de %s sale del dominio.', k, sc(s).name);
        end
    end
end
end
