function log_msg(fmt, varargin)
%LOG_MSG  Mensaje con marca de tiempo (se recoge también en el diario log.txt).
fprintf('[%s] %s\n', datestr(now, 'HH:MM:SS'), sprintf(fmt, varargin{:})); %#ok<TNOW1,DATST>
end
