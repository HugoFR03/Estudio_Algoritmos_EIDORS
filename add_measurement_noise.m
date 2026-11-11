function [Vn, sigma_n] = add_measurement_noise(V, snr_db, seed, ref_rms)
%ADD_MEASUREMENT_NOISE  Ruido gaussiano aditivo, blanco y reproducible.
%
%   [Vn, sigma_n] = add_measurement_noise(V, snr_db, seed, ref_rms)
%     V        adquisición limpia (vector columna)
%     snr_db   SNR en dB (Inf = sin ruido)
%     seed     semilla del generador (se usa un RandStream local, de modo
%              que el estado global de rng no se altera)
%     ref_rms  RMS de referencia; por defecto RMS(V), es decir,
%              sigma_n = RMS(V_limpia) / 10^(SNR/20)

if nargin < 4 || isempty(ref_rms)
    ref_rms = vrms(V);
end
if ~isscalar(snr_db) || isnan(snr_db)
    error('EIT:noise:snr', 'SNR no válido.');
end
if isinf(snr_db)
    Vn = V; sigma_n = 0;
    return
end
rs = RandStream('mt19937ar', 'Seed', seed);
sigma_n = ref_rms / 10^(snr_db / 20);
Vn = V + sigma_n * randn(rs, size(V));
end
