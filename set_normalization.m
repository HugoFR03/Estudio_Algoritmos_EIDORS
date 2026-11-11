function imdl = set_normalization(imdl, flag)
%SET_NORMALIZATION  Fija si los datos diferenciales se normalizan, (vi-vh)./vh,
%   o se usan en bruto, vi-vh. EIDORS consulta el indicador en el fwd_model
%   a través de mdl_normalize; en versiones sin esa función se escribe el
%   campo clásico normalize_measurements.
try
    imdl.fwd_model = mdl_normalize(imdl.fwd_model, flag);
catch
    imdl.fwd_model.normalize_measurements = flag;
end
end
