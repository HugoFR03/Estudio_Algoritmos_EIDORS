function [perm, sgn] = bp_orientation_permutation(stim, n_elec, k, s)
%BP_ORIENTATION_PERMUTATION  Reordenación de medidas equivalente a renumerar
%   los electrodos mediante e -> mod(s*(e-1) + k, n) + 1  (s = +1 rotación,
%   s = -1 reflexión seguida de rotación).
%
%   La matriz de Sheffield de mk_common_gridmdl('backproj') tiene su propia
%   numeración de electrodos. Para alinearla con el modelo de simulación se
%   prueban las 2*n transformaciones del grupo diédrico; cada una equivale a
%   reordenar (y, en su caso, cambiar de signo) las medidas:
%       v_transformada(i) = sgn(i) * v(perm(i))
%   La transformación correcta se elige en validación (ver
%   tune_hyperparameters). No se modifica ninguna otra cosa del algoritmo.

if ~ismember(s, [1 -1])
    error('EIT:bp:sign', 's debe ser +1 o -1.');
end
T = mod(s * ((1:n_elec) - 1) + k, n_elec) + 1;   % T(e): electrodo transformado

% Tabla de medidas: [patrón de corriente, patrón de medida]
K = zeros(0, 2 * n_elec);
for i = 1:numel(stim)
    c = full(stim(i).stim_pattern(:))';
    M = full(stim(i).meas_pattern);
    K = [K; repmat(c, size(M, 1), 1), M]; %#ok<AGROW>
end
nm = size(K, 1);
perm = zeros(nm, 1);
sgn  = zeros(nm, 1);
tol = 1e-9 * max(1, max(abs(K(:))));
for i = 1:nm
    c = zeros(1, n_elec); m = zeros(1, n_elec);
    c(T) = K(i, 1:n_elec);
    m(T) = K(i, n_elec+1:end);
    for sc = [1 -1]
        for sm = [1 -1]
            target = [sc * c, sm * m];
            j = find(all(abs(K - target) < tol, 2), 1);
            if ~isempty(j)
                perm(i) = j; sgn(i) = sc * sm;
                break
            end
        end
        if perm(i), break, end
    end
    if ~perm(i)
        error('EIT:bp:perm', ['El protocolo de medida no es invariante ante la transformación ', ...
            '(k=%d, s=%d): la medida %d no tiene equivalente.'], k, s, i);
    end
end
if numel(unique(perm)) ~= nm
    error('EIT:bp:perm2', 'La reordenación obtenida no es una permutación.');
end
end
