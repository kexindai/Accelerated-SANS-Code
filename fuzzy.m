% Fuzzy sphere model
function I_fuzzy = fuzzy(q, param)
scale = param(1);
R = param(2);
sig = param(3);
sld_solvent = 3;
sld_sphere = param(4);
b = param(5);
V = 4/3 * pi* (R)^3;
%R = R + sig;
Phi = 3 * V * (sin(q.*R) - q.*R.* cos(q.* R))./(q.*R).^3;
A = Phi.*exp(-0.5*(sig*q).^2);
I_fuzzy = (scale/V * (sld_sphere - sld_solvent).^2 .* A.^2)*10^(-4) + b;
end