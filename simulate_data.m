function sim = simulate_data(fmdl, ph, sigma_bkg)
%SIMULATE_DATA  Problema directo sobre la malla fina: medidas limpias.
%
%   sim.V0  tensiones del medio homogéneo de referencia (n_meas x 1)
%   sim.V1  tensiones con inclusiones, una columna por fantoma
%   sim.dV  datos diferenciales limpios, dV = V1 - V0
%
%   Convención única de signo y orden: el de fwd_solve con el patrón de
%   estimulación de fmdl. Ningún algoritmo recibe otra representación.

t = tic;
vh = fwd_solve(mk_image(fmdl, sigma_bkg));
V0 = vh.meas(:);
V1 = zeros(numel(V0), numel(ph));
for s = 1:numel(ph)
    if numel(ph(s).sigma) ~= size(fmdl.elems, 1)
        error('EIT:sim:dim', 'El fantoma %s no coincide con la malla directa.', ph(s).name);
    end
    vi = fwd_solve(mk_image(fmdl, ph(s).sigma));
    V1(:, s) = vi.meas(:);
end
if any(~isfinite(V0)) || any(~isfinite(V1(:)))
    error('EIT:sim:nan', 'El problema directo ha devuelto valores no finitos.');
end
sim.V0 = V0;
sim.V1 = V1;
sim.dV = V1 - V0;
sim.names = {ph.name};
sim.time  = toc(t);
for s = 1:numel(ph)
    log_msg('  %-20s RMS(dV)/RMS(V0) = %.4f (%.1f dB)', ph(s).name, ...
        vrms(sim.dV(:, s)) / vrms(V0), 20 * log10(vrms(sim.dV(:, s)) / vrms(V0)));
end
end
