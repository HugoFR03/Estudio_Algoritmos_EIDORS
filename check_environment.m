function env = check_environment(cfg)
%CHECK_ENVIRONMENT  Comprueba MATLAB, EIDORS y las funciones que usa el estudio.
%
%   env = check_environment(cfg) registra las versiones instaladas, la
%   disponibilidad y la ruta de cada función de EIDORS requerida y la
%   cabecera de su ayuda (para documentar la interfaz de la versión
%   instalada). Si falta una función imprescindible para alguno de los
%   seis algoritmos, se detiene con un error que indica cuál y por qué:
%   el estudio NO sustituye un algoritmo por otro.

env = struct();
env.timestamp      = datestr(now, 'yyyy-mm-dd HH:MM:SS'); %#ok<TNOW1,DATST>
env.matlab_version = version;
env.matlab_release = ['R' version('-release')];
env.computer       = computer;

% --- EIDORS en el path
if ~exist('eidors_obj', 'file') && ~isempty(cfg.eidors_startup)
    if ~exist(cfg.eidors_startup, 'file')
        error('EIT:env:startup', 'No existe el fichero de arranque de EIDORS: %s', cfg.eidors_startup);
    end
    run(cfg.eidors_startup);
end
if ~exist('eidors_obj', 'file')
    error('EIT:env:noEIDORS', ['EIDORS no está en el path de MATLAB. Ejecute eidors/startup.m ', ...
        'o indique su ruta en cfg.eidors_startup.']);
end
try
    env.eidors_version = eidors_obj('eidors_version');
catch
    env.eidors_version = 'desconocida';
end
try
    eidors_msg('log_level', cfg.eidors_log_level);
catch
    % versiones antiguas: se ignora
end

% --- Funciones a comprobar
fn = {'mk_common_model', 'mk_stim_patterns', 'mk_image', 'fwd_solve', 'inv_solve', ...
      'calc_jacobian', 'calc_RtR_prior', 'eidors_cache', 'mdl_normalize', ...
      'mk_common_gridmdl', 'prior_noser', 'inv_solve_diff_GN_one_step', ...
      'inv_solve_TSVD', 'calc_TSVD_RM', 'inv_solve_TV_pdipm', 'prior_TV', ...
      'mk_GREIT_model', 'find_boundary', 'inv_solve_core', 'inv_solve_gn', 'prior_tikhonov'};
env.functions = struct();
for i = 1:numel(fn)
    f = fn{i};
    info.available = exist(f, 'file') > 0;
    info.path      = which(f);
    info.help_head = '';
    if info.available
        try
            h = help(f);
            lines = regexp(h, '\n', 'split');
            info.help_head = strjoin(lines(1:min(25, numel(lines))), newline);
        catch
        end
    end
    env.functions.(f) = info;
end

% --- Requisitos por algoritmo
req.common = {'mk_common_model', 'mk_stim_patterns', 'mk_image', 'fwd_solve', ...
              'inv_solve', 'calc_jacobian', 'eidors_cache'};
req.BP     = {'mk_common_gridmdl'};
req.NOSER  = {'prior_noser', 'inv_solve_diff_GN_one_step'};
req.TSVD   = {'inv_solve_TSVD'};
req.TV     = {'inv_solve_TV_pdipm', 'prior_TV'};
switch lower(cfg.gn.implementation)
    case 'explicit'
        req.GN = {'calc_RtR_prior', 'prior_tikhonov'};
    case 'inv_solve_gn'
        req.GN = {'inv_solve_gn', 'inv_solve_core', 'prior_tikhonov'};
    otherwise
        error('EIT:env:gnImpl', 'cfg.gn.implementation desconocida: %s', cfg.gn.implementation);
end
req.GREIT  = {'mk_GREIT_model', 'find_boundary'};
env.requirements = req;

missing = {};
groups = [{'common'}, cfg.algorithms];
for g = 1:numel(groups)
    names = req.(groups{g});
    for i = 1:numel(names)
        if ~env.functions.(names{i}).available
            missing{end+1} = sprintf('%s (necesaria para %s)', names{i}, groups{g}); %#ok<AGROW>
        end
    end
end
env.missing = missing;

% Funciones opcionales: se informa, no se detiene
opt = {'calc_TSVD_RM', 'inv_solve_core', 'inv_solve_gn', 'mdl_normalize'};
env.optional_missing = opt(~cellfun(@(f) env.functions.(f).available, opt));

% --- Informe
log_msg('MATLAB %s (%s), EIDORS %s', env.matlab_release, env.matlab_version, env.eidors_version);
for i = 1:numel(fn)
    if env.functions.(fn{i}).available, st = 'OK'; else, st = 'NO DISPONIBLE'; end
    log_msg('  %-28s %s', fn{i}, st);
end
if ~isempty(env.optional_missing)
    log_msg('Aviso: funciones opcionales ausentes: %s', strjoin(env.optional_missing, ', '));
end
if ~isempty(missing)
    error('EIT:env:missing', ['Faltan funciones de EIDORS imprescindibles:\n  %s\n', ...
        'El estudio no sustituye algoritmos. Revise la instalación o la versión de EIDORS.'], ...
        strjoin(missing, '\n  '));
end
end
