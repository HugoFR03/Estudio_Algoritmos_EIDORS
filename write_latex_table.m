function write_latex_table(T, fname, caption, label, numfmt)
%WRITE_LATEX_TABLE  Exporta una tabla de MATLAB a LaTeX (booktabs).
%   Los números se escriben con numfmt; el texto se escapa (_ % & #),
%   salvo las expresiones matemáticas entre $...$ que ya traiga.
if nargin < 5, numfmt = '%.3f'; end
fid = fopen(fname, 'w', 'n', 'UTF-8');
if fid < 0, error('EIT:latex:open', 'No se puede escribir %s', fname); end
c = onCleanup(@() fclose(fid));
nc = width(T);
fprintf(fid, '%% Generado automáticamente por generate_results.m\n');
fprintf(fid, '\\begin{table}[htbp]\n\\centering\n\\small\n');
fprintf(fid, '\\caption{%s}\n\\label{%s}\n', caption, label);
fprintf(fid, '\\begin{tabular}{l%s}\n\\toprule\n', repmat('c', 1, nc - 1));
fprintf(fid, '%s \\\\\n\\midrule\n', strjoin(cellfun(@esc, T.Properties.VariableNames, ...
    'UniformOutput', false), ' & '));
for i = 1:height(T)
    cells = cell(1, nc);
    for j = 1:nc
        v = T{i, j};
        if iscell(v), v = v{1}; end
        if isnumeric(v) || islogical(v)
            if isnan(v), cells{j} = '--';
            elseif isinf(v), cells{j} = '$\infty$';
            else, cells{j} = sprintf(numfmt, v);
            end
        else
            cells{j} = esc(char(v));
        end
    end
    fprintf(fid, '%s \\\\\n', strjoin(cells, ' & '));
end
fprintf(fid, '\\bottomrule\n\\end{tabular}\n\\end{table}\n');
end

function s = esc(s)
s = char(s);
s = strrep(s, '_', '\_');
s = strrep(s, '%', '\%');
s = strrep(s, '&', '\&');
s = strrep(s, '#', '\#');
end
