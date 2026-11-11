function save_figure(f, fbase, cfg)
%SAVE_FIGURE  Guarda una figura editable (.fig) y su exportación PNG.
d = fileparts(fbase);
if ~isempty(d) && ~exist(d, 'dir'), mkdir(d); end
% Las figuras se crean ocultas; al abrir el .fig se mostrarán.
set(f, 'CreateFcn', 'set(gcbo, ''Visible'', ''on'')');
savefig(f, [fbase '.fig']);
try
    exportgraphics(f, [fbase '.png'], 'Resolution', cfg.fig.dpi);
catch
    print(f, [fbase '.png'], '-dpng', sprintf('-r%d', cfg.fig.dpi));
end
close(f);
end
