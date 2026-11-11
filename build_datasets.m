function ds = build_datasets(sim, snr_db, labels, n_real, seed_base, reference)
%BUILD_DATASETS  Conjuntos de medidas (vh, vi) para todos los algoritmos.
%
%   ds = build_datasets(sim, snr_db, labels, n_real, seed_base, reference)
%   genera, para cada fantoma y nivel de ruido, n_real realizaciones (una
%   sola si el nivel es "sin ruido"). El ruido se añade por separado a V0 y
%   a V1 ANTES de calcular la diferencia, con semillas deterministas:
%       seed_V0 = seed_base + 100000*s + 1000*n + 2*r
%       seed_V1 = seed_V0 + 1
%   Los conjuntos se generan una sola vez y se entregan idénticos a los
%   seis algoritmos.
%
%   reference = 'acquisition': sigma por adquisición = RMS(V_limpia)/10^(SNR/20)
%   reference = 'difference' : sigma = RMS(dV_limpio)/(sqrt(2)*10^(SNR/20))

nS = size(sim.V1, 2);
C = {};
for s = 1:nS
    V1c = sim.V1(:, s);
    dvc = V1c - sim.V0;
    for n = 1:numel(snr_db)
        if isinf(snr_db(n)), nr = 1; else, nr = n_real; end
        for r = 1:nr
            seed0 = seed_base + 100000 * s + 1000 * n + 2 * r;
            seed1 = seed0 + 1;
            switch lower(reference)
                case 'acquisition'
                    ref0 = vrms(sim.V0); ref1 = vrms(V1c);
                case 'difference'
                    ref0 = vrms(dvc) / sqrt(2); ref1 = ref0;
                otherwise
                    error('EIT:noise:ref', 'Referencia de ruido desconocida: %s', reference);
            end
            [vh, s0] = add_measurement_noise(sim.V0, snr_db(n), seed0, ref0);
            [vi, s1] = add_measurement_noise(V1c, snr_db(n), seed1, ref1);
            dV = vi - vh;
            nd = dV - dvc;
            if any(nd), snr_eff = 20 * log10(vrms(dvc) / vrms(nd)); else, snr_eff = Inf; end
            C{end+1} = struct('scenario', s, 'scenario_name', sim.names{s}, ...
                'noise_idx', n, 'noise_label', labels{n}, 'snr_db', snr_db(n), ...
                'realization', r, 'seed_v0', seed0, 'seed_v1', seed1, ...
                'sigma_v0', s0, 'sigma_v1', s1, 'snr_eff_dv_db', snr_eff, ...
                'vh', vh, 'vi', vi); %#ok<AGROW>
        end
    end
end
ds = [C{:}];
end
