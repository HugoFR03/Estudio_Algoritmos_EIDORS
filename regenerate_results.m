function tables = regenerate_results(out_dir)
%REGENERATE_RESULTS  Rehace figuras y tablas a partir de un results.mat ya
%   guardado, sin repetir el ajuste ni las reconstrucciones.
%
%   tables = regenerate_results('results/20261108_101500')
S = load(fullfile(out_dir, 'results.mat'));
tables = generate_results(S.res, S.ds_test, S.ph_test, S.Timg_test, S.tuning, ...
    S.grid, S.cfg, out_dir);
save(fullfile(out_dir, 'results.mat'), 'tables', '-append');
log_msg('Figuras y tablas regeneradas en %s', out_dir);
end
