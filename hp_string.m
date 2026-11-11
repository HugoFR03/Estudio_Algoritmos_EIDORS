function s = hp_string(hp)
%HP_STRING  Texto compacto con los hiperparámetros de un algoritmo.
skip = {'RIE_mean', 'pearson_mean', 'time_recon_median', 'iterations_mean', ...
        'time_train', 'transform_k', 'transform_s', 'scale'};
f = fieldnames(hp);
p = {};
for i = 1:numel(f)
    if ismember(f{i}, skip), continue, end
    v = hp.(f{i});
    if isnumeric(v) || islogical(v)
        p{end+1} = sprintf('%s=%s', f{i}, mat2str(v, 4)); %#ok<AGROW>
    end
end
s = strjoin(p, ', ');
end
