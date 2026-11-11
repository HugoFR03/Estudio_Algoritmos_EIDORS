function sc = make_validation_scenarios(cfg)
%MAKE_VALIDATION_SCENARIOS  Fantomas aleatorios del conjunto de validación.
%
%   Genera cfg.val.n_phantoms escenarios con 1..cfg.val.max_inclusions
%   inclusiones de radio, posición y conductividad aleatorios, con una
%   semilla propia (cfg.val.seed). Se rechaza cualquier escenario que
%   reproduzca uno de prueba (mismo número de inclusiones y centros a menos
%   de 0,1), de modo que el ajuste nunca ve los casos de evaluación.

rs = RandStream('mt19937ar', 'Seed', cfg.val.seed);
test = cfg.scenarios;
C = cell(1, cfg.val.n_phantoms);
p = 0; tries = 0;
while p < cfg.val.n_phantoms
    tries = tries + 1;
    if tries > 10000
        error('EIT:val:gen', 'No se han podido generar fantomas de validación válidos.');
    end
    nin = randi(rs, cfg.val.max_inclusions);
    inc = struct('center', {}, 'radius', {}, 'sigma', {});
    ok = true;
    for k = 1:nin
        placed = false;
        for attempt = 1:200
            r = cfg.val.radius_range(1) + diff(cfg.val.radius_range) * rand(rs);
            rmax = cfg.val.max_extent - r;
            rho = rmax * sqrt(rand(rs));
            th = 2 * pi * rand(rs);
            c = rho * [cos(th) sin(th)];
            if all(arrayfun(@(q) norm(q.center - c) >= q.radius + r + cfg.val.min_gap, inc))
                placed = true; break
            end
        end
        if ~placed, ok = false; break, end
        sig = cfg.val.sigma_choices(randi(rs, numel(cfg.val.sigma_choices)));
        inc(end+1) = struct('center', c, 'radius', r, 'sigma', sig); %#ok<AGROW>
    end
    if ~ok || replicates_test(inc, test), continue, end
    p = p + 1;
    C{p} = struct('name', sprintf('V%02d', p), ...
        'label', sprintf('Validación %d (%d inclusiones)', p, nin), 'inclusions', inc);
end
sc = [C{:}];
end

function tf = replicates_test(inc, test)
tf = false;
for s = 1:numel(test)
    ti = test(s).inclusions;
    if numel(ti) ~= numel(inc), continue, end
    d = zeros(numel(inc), 1);
    for k = 1:numel(inc)
        d(k) = min(arrayfun(@(q) norm(q.center - inc(k).center), ti));
    end
    if all(d < 0.1), tf = true; return, end
end
end
