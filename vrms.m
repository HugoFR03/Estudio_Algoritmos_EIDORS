function r = vrms(x)
%VRMS  Valor eficaz de un vector (sin depender de toolboxes).
r = sqrt(mean(x(:).^2));
end
