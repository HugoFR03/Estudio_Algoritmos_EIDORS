function clear_eidors_cache()
%CLEAR_EIDORS_CACHE  Vacía la caché de EIDORS (matrices y jacobianos precalculados).
%   Se usa antes de medir tiempos de precomputación y de entrenamiento.
try
    eidors_cache('clear');
catch
    warning('EIT:cache', 'No se ha podido vaciar la caché de EIDORS.');
end
end
