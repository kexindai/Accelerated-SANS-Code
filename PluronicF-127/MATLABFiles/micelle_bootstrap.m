function [saveparam_spheremicelle, I_fit_micelle, residual, ci] = micelle_bootstrap(M_data, I_rep_mat)
saveparam_spheremicelle = zeros(1,18);  
q = M_data(22:end-2,1);
I = I_rep_mat(22:end-2,1);
dI = M_data(22:end-2,3);
% Fit to spherical micelle with polydispersity
ndensity = 10;
Rc = 50;
Rg = 30;
d = 1;
bkg = 0.01;
sld_core = 3;
sld_shell = 4;
polydisp = 0.18;
% initial guess
x0 = [ndensity, Rc, Rg, d, bkg, sld_core, sld_shell, polydisp];
lb = [0, 20, 10, 0.5, -inf,  0, 0, 0];
ub = [inf, inf, inf, 2, inf, inf, inf, 1];
numparams = 8;
options=optimoptions('lsqnonlin','MaxFunEvals',inf,'MaxIter',inf,...
         'TolFun',1e-12,'TolX',1e-12);
% A = [0, -1, 1, 0, 0, 0, 0];
% b = 0; 
% Aeq = [];
% beq = [];
% [x_micelle,resnorm,residual,~,~,~,jacobian] = ...
%           lsqnonlin(@(x) LSfunc(@(q,x) micelle(q, x),q,I,dI,x),x0,lb,ub, A,b, Aeq,beq, [], options);
[x_micelle,resnorm,residual,~,~,~,jacobian] = ...
        lsqnonlin(@(x) LSfunc(@(q,x) micelle_SZ(q, x),q,I,dI,x),x0,lb,ub, options);
ci = nlparci(x_micelle,residual,'jacobian',jacobian);
RSS_N = sum((residual.*I).^2)/length(residual);
AIC_micelle = length(residual)*log(RSS_N) + 2*numparams;
BIC_micelle = length(residual)*log(RSS_N) + numparams * log(length(residual));
[I_fit_micelle] = micelle_SZ(q,x_micelle);
 saveparam_spheremicelle(1 ,1:8) = x_micelle';
 saveparam_spheremicelle(1,9:16) = (x_micelle' - ci(:,1))';
 saveparam_spheremicelle(1,17) = AIC_micelle;
 saveparam_spheremicelle(1,18) = BIC_micelle;
end

function LS = LSfunc(func,q,I_expt,dI_expt,params)
%Obtain fit I data from correlation peak model (CorrLength)
I_fit = func(q,params); % any function in this code!!! :D :D :D 
%Calculate the residual
diff = I_fit - I_expt;
%Create the LS objective function
LS = diff./I_expt;
end

%% Model: Polymer micelle model with dispersity

function [I_micelle_SZ] = micelle_SZ(q, x)
ndensity = x(1);
Rc_mean = x(2);       % Mean core radius
Rg = x(3);
d = x(4);
b = x(5);
sld_core = x(6);
sld_corona = x(7);
sld_solvent = 6.33;
v_core = 6283;
v_corona = 14667;
pdisp = x(8);         % Polydispersity index (sigma / Rc_mean)

% Schulz-Zimm shape parameter
z = (1 / pdisp^2) - 1;

% Schulz-Zimm distribution over Rc
n_points = 50;
Rc_values = linspace((Rc_mean) * 0.1, (Rc_mean) * 3.0, n_points);  % Reasonable range
dRc = Rc_values(2) - Rc_values(1);

% Schulz-Zimm weights
D = ((z+1).^(z+1) ./ gamma(z+1)) .* (Rc_values / Rc_mean).^z ...
    .* exp(-(z+1) * Rc_values / Rc_mean) ./ Rc_mean;
weights = D * dRc;
weights = weights / sum(weights);  % Normalize

I_q = zeros(size(q));

for i = 1:n_points
    Rc = Rc_values(i);
    weight = weights(i);

    n_agg = 4/3 * pi * Rc^3 / v_core;
    beta_s = v_core * (sld_core - sld_solvent);
    beta_c = v_corona * (sld_corona - sld_solvent);

    Phi = 3 * (sin(q.*Rc) - q.*Rc.*cos(q.*Rc)) ./ (q.*Rc).^3;
    Fs = Phi.^2;

    xq = (q .* Rg).^2;
    Fc = 2 * (exp(-xq) - 1 + xq) ./ xq.^2;
    psi = (1 - exp(-xq)) ./ xq;

    Ssc = Phi .* psi .* (sin(q .* (Rc + d*Rg))) ./ (q .* (Rc + d*Rg));
    Scc = psi.^2 .* ((sin(q .* (Rc + d*Rg))) ./ (q .* (Rc + d*Rg))).^2;

    term1 = n_agg.^2 * beta_s.^2 .* Fs;
    term2 = n_agg * beta_c.^2 .* Fc;
    term3 = 2 * n_agg.^2 * beta_s * beta_c .* Ssc;
    term4 = n_agg * (n_agg - 1) * beta_c.^2 .* Scc;

    P_micelle = term1 + term2 + term3 + term4;

    I_q = I_q + weight * P_micelle;
end

I_micelle_SZ = I_q * ndensity * 1e-13 + b;
end