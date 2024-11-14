function [saveparam_spheremicelle, I_fit, residual, ci] = micelle_bootstrap(M_data, I_rep_mat)
q = M_data(58:end,1);
I = I_rep_mat(58:end);
dI = M_data(58:end,3);
% Fit to spherical micelle
ndensity = 100;
Rc = 30;
Rg = 10;
d = 1;
A = 1;
b = 0;
n_agg = 2;
sld_core = 0.344;
sld_shell = 0.638;
x0 = [ndensity; Rc; Rg; d;  A; b; sld_core; sld_shell];
lb = [0; 0; 5; 1; 0; 0; 0;0];
ub = [inf ; inf; inf; 2; inf; inf;inf;inf];
options=optimoptions('lsqnonlin','MaxFunEvals',inf,'MaxIter',inf,...
         'TolFun',1e-12,'TolX',1e-12);
[x_micelle,~,residual,~,~,~,jacobian] = ...
        lsqnonlin(@(x) LSfunc(@(q,x) micelle(q, x),q,I,dI,x),x0,lb,ub, options);
ci = nlparci(x_micelle,residual,'jacobian',jacobian);
RSS_N = sum((residual.*dI).^2)/length(residual);
AIC_micelle = length(residual)*log(RSS_N) + 2*7;
BIC_micelle = length(residual)*log(RSS_N) + 7 * log(length(residual));
saveparam_spheremicelle(1,1:8) = x_micelle';
saveparam_spheremicelle(1,9:16) = (x_micelle - ci(:,1))';
saveparam_spheremicelle(1,17) = AIC_micelle;
saveparam_spheremicelle(1,18) = BIC_micelle;
I_fit = micelle(q,x_micelle);
end

function LS = LSfunc(func,q,I_expt,dI_expt,params)
%Obtain fit I data from correlation peak model (CorrLength)
I_fit = func(q,params); % any function in this code!!! :D :D :D 
%Calculate the residual
diff = I_fit - I_expt;
%Create the LS objective function
LS = diff./dI_expt;
end
%% Model: Polymer micelle model

function I_micelle = micelle(q,x)
ndensity = x(1);
Rc = x(2);
Rg = x(3);
d = x(4);
A = x(5);
b = x(6);

% SLD and v for F127
sld_core = x(7);
sld_corona = x(8);
sld_solvent = 6.33;
v_core = 6283;
v_corona = 14667;
n_agg = 4/3 * pi *Rc^3/v_corona;
%n_agg = 4/3 * pi* Rc^3/v_core;
% Example testing
% sld_core =  0.34;
% sld_corona = 0.8;
% v_core = 62624;
% v_corona = 61940;
% Pederson paper
% sld_core =  6.3 +1;
% sld_corona = 6.3 +1 ;
% sld_solvent = 6.3 ;
% v_corona = 4000;
% v_core = 4000;
beta_s = v_core * (sld_core - sld_solvent);
beta_c = v_corona * (sld_corona - sld_solvent);
Phi = 3 * (sin(q.*Rc) - q.*Rc.* cos(q.* Rc))./(q.*Rc).^3;
Fs = Phi.^2;
x = (q.*Rg).^2;
Fc = 2 * (exp(-x) - 1 + x)./x.^2;
psi = (1 - exp(-x))./(x);
Ssc = Phi.*psi.* (sin(q.*(Rc + d*Rg)))./(q*(Rc+d*Rg));
Scc = psi.^2.* ((sin(q.*(Rc + d*Rg)))./(q*(Rc+d*Rg))).^2;
term1 = n_agg.^2 * beta_s.^2* Fs;
term2 = n_agg * beta_c.^2*Fc;
term3 = 2 * n_agg.^2 * beta_s * beta_c * Ssc;
term4 = n_agg * (n_agg - 1) * beta_c.^2 * Scc;
P_micelle = term1 + term2+ term3 + term4;
I_micelle = (A*P_micelle + b)*ndensity*10^(-13);
end