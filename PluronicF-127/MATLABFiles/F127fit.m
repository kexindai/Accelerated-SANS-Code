clc
clear
close all
k = 1;
saveparam_spheremicelle = zeros(k,18);  
saveparam_sphere = zeros(k,10);  
saveparam_fuzzy = zeros(k,7);  
%% load the data
for i = 1:k
clearvars -except i saveparam_spheremicelle saveparam_sphere saveparam_fuzzy
F127 = importdata(sprintf('F127_%d_merged_1.txt', i-1));

q = F127(22:end-2,1);
I = F127(22:end-2,2);
dI = F127(22:end-2,3);

% Fit to spherical micelle

ndensity = 10;
Rc = 50;
Rg = 30;
d = 1;
bkg = 0.004;
sld_core = 3;
sld_shell = 4;
polydisp = 0.18;
numparams = 8;
%% this is for fitting just micelle model
% x0 = [ndensity, Rc, Rg, d, bkg, sld_core, sld_shell];
% lb = [0, 20, 10, 0.5, -inf,  0, 0];
% ub = [inf, inf, inf, 2, inf, inf, inf];

% fitting micelle model with polydispersity
x0 = [ndensity, Rc, Rg, d, bkg, sld_core, sld_shell, polydisp];
lb = [0, 20, 10, 0.5, -inf,  0, 0, 0];
ub = [inf, inf, inf,1.5, inf, inf, inf, 1];
figure('Visible','off')
I_micelle = micelle_SZ(q,x0);
errorbar(q*10, I, dI, 's','Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
hold on
loglog(q*10, I_micelle, '-','LineWidth',3,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
xticks([10.^(-1) 10.^(0) 10.^1])
xticklabels({'10^{-1}','1','10^{1}'})
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
xlabel('q [nm^-^1]','FontWeight','bold');
ylabel('I(q) [cm^-^1]','FontWeight','bold');
%%
options=optimoptions('lsqnonlin','MaxFunEvals',inf,'MaxIter',inf,...
         'TolFun',1e-12,'TolX',1e-12);
% A = [0, -1, 1, 0, 0, 0, 0, 0];
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

figure('Visible','on')
%ax.Box = 'on';
errorbar(q*10, I, dI, 's','LineWidth',2,'Color','k','MarkerEdgeColor','k', 'MarkerSize',12)
hold on
plot(q*10, I_fit_micelle, 'LineWidth',2)

set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);

xlabel('q [nm^-^1]','FontWeight','bold');
ylabel('I(q) [cm^-^1]','FontWeight','bold');
xticks([10.^(-1) 10.^(0) 10.^1])
xticklabels({'10^{-1}','1','10^{1}'})
xlim([0.046 5.33])
ylim([10^(-3) 10^2])
%% saving parameters
 saveparam_spheremicelle(i ,1:8) = x_micelle';
 saveparam_spheremicelle(i,9:16) = (x_micelle' - ci(:,1))';
 saveparam_spheremicelle(i,17) = AIC_micelle;
 saveparam_spheremicelle(i,18) = BIC_micelle;
 %writematrix(saveparam_spheremicelle, 'F127_newmicelle_parameters_1.csv')
 saveas(gcf,sprintf('F127_%d_fit_spheremicelle_0p5.png',i-1))

%% Bias Testing
%% this section tests for the correlation between errors in Pluronics F-127 micelles
% Work flow
% 1. Take the model and the fitted parameters f(q, theta) where theta is
% the parmaters
% 2. Fix theta and add noise to the model. We will first test Gaussian
% noise (the assumption for MC bootstrapping)
% 3. We will create 100 replicates and fit the model with broad peak model
% 4. If nonlinear least squares is a biased estimator (errors are correlated), 
% then your theta’ will have expectation something other than theta

% the noise is gaussian centered at 0 but porportional to the dI
% the noise is gaussian centered at 0 but porportional to the dI
Num_rep = 100;
noise = dI.*randn(length(q), Num_rep);
I_rep_data = repmat(I_fit_micelle, [1, Num_rep]) + noise;
for i = 1:Num_rep
[x_micelle_test(i,:),resnorm,residual,~,~,~,jacobian] = ...
        lsqnonlin(@(x) LSfunc(@(q,x) micelle_SZ(q, x),q, I_rep_data(:,i),dI,x),x0,lb,ub, options);
end
%% calculate bias and plot histogram
E = mean(x_micelle_test);
bias = (E - x_micelle)./x_micelle * 100; 
figure()
histogram(x_micelle_test(:,1),'LineWidth',2,'FaceColor',[0 0.4470 0.7410]);
xlabel('N_{density} values')
ylabel('Frequency')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
saveas(gcf, 'F-127 Ndensity distribution.png')
figure()
histogram(x_micelle_test(:,2),'LineWidth',2);
xlabel('R values')
ylabel('Frequency')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
saveas(gcf, 'F-127 R distribution.png')
figure()
histogram(x_micelle_test(:,3),'LineWidth',2);
xlabel('R_{g} values')
ylabel('Frequency')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
saveas(gcf, 'F-127 Rg distribution.png')
figure()
histogram(x_micelle_test(:,4),'LineWidth',2);
xlabel('d values')
ylabel('Frequency')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
saveas(gcf, 'F-127 d distribution.png')
figure()
histogram(x_micelle_test(:,5),'LineWidth',2);
xlabel('b values')
ylabel('Frequency')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
saveas(gcf, 'F-127 b distribution.png')

figure()
histogram(x_micelle_test(:,8),'LineWidth',2);
xlabel('Polydispersity values')
ylabel('Frequency')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
saveas(gcf, 'F-127 Polydispersity distribution.png')


%% Fit to sphere model
% scale = 1;
% R = 50;
% sld_sphere = 5;
% b = 0.01;
% x0 = [scale; R; sld_sphere; b];
% lb = [0;0;0;0];
% ub = [inf ; inf; inf; inf];
% options=optimoptions('lsqnonlin','MaxFunEvals',inf,'MaxIter',inf,...
%          'TolFun',1e-12,'TolX',1e-12);
% [x_sphere,resnorm,residual,~,~,~,jacobian] = ...
%         lsqnonlin(@(param) LSfunc(@(q,param) sphere(q, param),q,I,dI,param),x0,lb,ub, options);
% ci = nlparci(x_sphere,residual,'jacobian',jacobian);
% q_fit = logspace(-3,1,200);
% I_fit = sphere(q_fit,x_sphere);
% figure('visible', 'off')
% errorbar(q*10, I, dI, 'o','LineWidth',2,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
% hold on
% plot(q_fit*10,I_fit', 'Color','r','LineWidth',2)
% xlabel('q [nm^-^1]','FontWeight','bold');
% ylabel('I(q) [cm^-^1]','FontWeight','bold');
% xticks([10.^(-1) 10.^(0) 10.^1])
% xticklabels({'10^{-1}','1','10^{1}'})
% xlim([0.046 5.33])
% ylim([10^(-3) 10^2])
% set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
% set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
% saveas(gcf,sprintf('F127_%d_fit_sphere.png',i-1))
% RSS_N = sum((residual.*I).^2)/length(residual);
% AIC_sphere = length(residual)*log(RSS_N) + 2*4;
% BIC_sphere = length(residual)*log(RSS_N) + 4 * log(length(residual));
% saveparam_sphere(i ,1:4) = x_sphere';
% saveparam_sphere(i,5:8) = (x_sphere - ci(:,1))';
% saveparam_sphere(i,9) = AIC_sphere;
% saveparam_sphere(i,10) = BIC_sphere;
% 
% %% Fit to fuzzy sphere
% scale = 1;
% R = 50;
% sig = 1;
% sld_sphere = 5;
% b = 0.004;
% x0 = [scale; R; sig; sld_sphere; b];
% lb = [0;30;0;0;0];
% ub = [inf ; inf; inf; inf;inf];
% options=optimoptions('lsqnonlin','MaxFunEvals',inf,'MaxIter',inf,...
%          'TolFun',1e-12,'TolX',1e-12);
% [x_fuzzy,resnorm,residual,~,~,~,jacobian] = ...
%         lsqnonlin(@(param) LSfunc(@(q,param) fuzzy(q, param),q,I,dI,param),x0,lb,ub, options);
% ci = nlparci(x_fuzzy,residual,'jacobian',jacobian);
% q_fit = logspace(-3,1,200);
% I_fit = fuzzy(q_fit,x_fuzzy);
% figure('visible', 'off')
% errorbar(q*10, I, dI, 'o','LineWidth',2,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
% hold on
% plot(q_fit*10,I_fit', 'Color','r','LineWidth',2)
% xlabel('q [nm^-^1]','FontWeight','bold');
% ylabel('I(q) [cm^-^1]','FontWeight','bold');
% xlim([0.13 4.5])
% ylim([0 10])
% xticks([10.^(-1) 10.^(0) 10.^1])
% xticklabels({'10^{-1}','1','10'})
% xlim([0.046 5.33])
% ylim([10^(-3) 10^2])
% set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
% set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
% saveas(gcf,sprintf('F127_%d_fit_fuzzy_sphere.png',i-1))
% RSS_N = sum((residual.*I).^2)/length(residual);
% AIC_fuzzy = length(residual)*log(RSS_N) + 2*5;
% BIC_fuzzy = length(residual)*log(RSS_N) + 5 * log(length(residual));
% saveparam_fuzzy(i ,1:5) = x_fuzzy';
% saveparam_fuzzy(i,6:10) = (x_fuzzy - ci(:,1))';
% saveparam_fuzzy(i,11) = AIC_fuzzy;
% saveparam_fuzzy(i,12) = BIC_fuzzy;

end

%% save the parameters

% T_spheremicelle = table(saveparam_spheremicelle);
% T_spheremicelle_output = splitvars(T_spheremicelle,'saveparam_spheremicelle','NewVariableNames',{'ndensity', 'Rc', 'Rg', 'd', 'b', 'sld_core', 'sld_shell','polydisp',...
%      'ndensity_Err', 'Rc_Err', 'Rg_Err', 'd_Err','A_Err', 'b_Err', 'sld_core_Err', 'sld_shell_Err','AIC','BIC'});
% writetable( T_spheremicelle_output, 'F127_spheremicelle_parameters_0p5.csv')
% csvwrite('F127_sphere_parameters_0p01.csv',saveparam_sphere)
% csvwrite('F127_fuzzy_parameters_0p01.csv',saveparam_fuzzy)

%% compute the average and std of key parameters
% data = readtable('F127_spheremicelle_parameters_0p5.csv');
% data = table2array(data);
% Rc_avg = mean(data(:,2));
% Rc_std = std(data(:,2));
% Rg_avg = mean(data(:,3));
% Rg_std = std(data(:,3));
%% LS function to minimize

function LS = LSfunc(func,q,I_expt,dI_expt,params)
%Obtain fit I data from correlation peak model (CorrLength)
I_fit = func(q,params); % any function in this code!!! :D :D :D 
%Calculate the residual
diff = I_fit - I_expt;
%Create the LS objective function
LS = diff./(I_expt);

end

%% Model: Polymer micelle model

function [I_micelle] = micelle(q,x)
ndensity = x(1);
Rc = x(2);
Rg = x(3);
d = x(4);
b = x(5);
% SLD and v for F127
sld_core = x(6);
sld_corona = x(7);
sld_solvent = 6.33;
v_core = 6283;
v_corona = 14667;
n_agg = 4/3 * pi* Rc^3/v_core;
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
I_micelle = (P_micelle) *ndensity*10^(-13) + b;
end

%% Micelle model with polydispersity Gaussian
function [I_micelle_Gauss] = micelle_Gauss(q, x)
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
pdisp = x(8);         % Polydispersity (sigma / Rc)

% Discretize Rc over Gaussian distribution
n_points = 40;  % Integration resolution
sigma = pdisp * Rc_mean;
Rc_values = linspace(Rc_mean - 3*sigma, Rc_mean + 3*sigma, n_points);
weights = exp(-(Rc_values - Rc_mean).^2 ./ (2 * sigma^2));
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

I_micelle_Gauss = I_q * ndensity * 1e-13 + b;
end

%% lognormal distribution
function [I_micelle_p] = micelle_p(q, x)
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
pdisp = x(8);         % Polydispersity (σ / Rc_mean)

% Log-normal distribution parameters
sigma_ln = sqrt(log(1 + pdisp^2));
mu_ln = log(Rc_mean) - 0.5 * sigma_ln^2;

% Discretize Rc over log-normal distribution
n_points = 50;
Rc_values = linspace(Rc_mean * 0.3, Rc_mean * 2.5, n_points);  % Sample around mean
dRc = Rc_values(2) - Rc_values(1);  % For proper integration weights
lognorm_weights = (1 ./ (Rc_values * sigma_ln * sqrt(2*pi))) .* ...
    exp(-(log(Rc_values) - mu_ln).^2 / (2 * sigma_ln^2));
lognorm_weights = lognorm_weights / sum(lognorm_weights);  % Normalize

I_q = zeros(size(q));

for i = 1:n_points
    Rc = Rc_values(i);
    weight = lognorm_weights(i);

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

I_micelle_p = I_q * ndensity * 1e-13 + b;
end

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

%% Core-shell ellipsoid micelle model
function I_ellipsoid = coreshell_ellipsoid(q, param)
R = param(1);
eps = param(2);
d = param(3);
Rg = param(4);
A = param(5);
b = param(6);
ndensity = param(7);
sld_core = param(8);
sld_corona = param(9);
% SLD and v input
% sld_core =  (6.3 +1);
% sld_corona = 6.3 ;
% sld_solvent = 6.3 ;
% v_core = 4000;
% v_corona = 4000;
% F127 parameter
sld_solvent = 6.33;
v_core = 6283;
v_corona = 14667;
n_agg = v_core/v_corona;
beta_s = v_core * (sld_core - sld_solvent);
beta_c = v_corona * (sld_corona - sld_solvent);
x = (q.*Rg).^2;
Fc = 2 * (exp(-x) - 1 + x)./x.^2;
tspan = [0 pi/2];
options = odeset('RelTol',1e-8,'AbsTol',1e-8);
for i = 1:length(q)
psi(i) = (1 - exp(-q(i).*Rg))./(q(i).*Rg);
[~,y] = ode45(@(alpha,int_Fs_ellipsoid) Fsfun_ellipsoid(q(i), alpha, R, eps), tspan, 0, options);
Fs(i) = y(end,1);
Ssc(i) = psi(i).*integral(@(alpha) Ssc_fun_ellipsoid(q (i), alpha, R, Rg, d,eps), 0, pi/2);
Scc(i) = psi(i).^2*integral(@(alpha) Scc_fun_ellipsoid(q(i), R, eps,alpha, d, Rg), 0, pi/2);
end

term1 = n_agg.^2 * beta_s.^2* Fs;
term2 = n_agg * beta_c.^2.*Fc';
term3 = 2 * n_agg.^2 * beta_s * beta_c.* Ssc;
term4 = n_agg* (n_agg - 1) * beta_c.^2 * Scc;
P_ellipsoid = term1 + term2+ term3 + term4;
I_ellipsoid = (A*P_ellipsoid + b)'*ndensity*10^(-13);
end

%% Core-shell cylindrical micelle model
function I_cylinder = coreshell_cylinder(q, param)
R = param(1);
Rg = param(2);
L = param(3);
d = param(4);
A = param(5);
b = param(6);
ndensity = param(7);
sld_core = param(8);
sld_corona = param(9);
% SLD and v for F127
sld_solvent = 6.33;
v_core = 6283;
v_corona = 14667;
% For example testing
% sld_core =  6.3 +1;
% sld_corona = 6.3+1  ;
% sld_solvent = 6.3 ;
% v_core = 4000;
% v_corona = 4000;
n_agg = pi*R^2*L/v_corona;
beta_s = v_core * (sld_core - sld_solvent);
beta_c = v_corona * (sld_corona - sld_solvent);
for i = 1:length(q)
x(i) = (q(i)*Rg).^2;
Fc(i) = 2 * (exp(-x(i)) - 1 + x(i))/(x(i)^2);
psi(i) = (1 - exp(-q(i)*Rg))./(q(i)*Rg);
Fs(i) = integral(@(alpha) Fsfun_cylinder(q(i), alpha, R, L), 0, pi/2);
Ssc(i) = psi(i)*integral(@(alpha) Ssc_fun_cylinder(q(i), alpha, R, Rg, L, d), 0, pi/2);
Scc(i) = psi(i)^2*integral(@(alpha) Scc_fun_cylinder(q(i), alpha, R, Rg, L, d), 0, pi/2);
end
% When L and R are large the scale of term 2 dominates the others
term1 = n_agg.^2 * beta_s.^2* Fs;
term2 = n_agg * beta_c.^2 * Fc;
term3 = 2 * n_agg.^2 * beta_s * beta_c * Ssc;
term4 = n_agg * (n_agg - 1) * beta_c.^2 * Scc;
P_cylinder = term1 + term2+ term3 + term4;
%I_cylinder = (term1 + term2 + term3 + term4);
% loglog(q,term1)
% hold on
% loglog(q,term2)
% hold on
% loglog(q,term3)
% hold on
% loglog(q,term4)
% hold on
% legend('term1','term2','term3','term4')
I_cylinder = (A*P_cylinder + b)' * ndensity * 10^(-13);
end

%% Worm like micelle model
function I_worm = coreshell_worm(q, param)
R = param(1);
Rg = param(2);
L = param(3);
d = param(4);
A = param(5);
b = param(6);
ndensity = param(7);
sld_core = param(8);
sld_corona = param(9);
% SLD and v for F127

sld_solvent = 6.33;
v_core = 6283;
v_corona = 14667;
% For example testing
% sld_core =  6.3 +1;
% sld_corona = 6.3 +1;
% sld_solvent = 6.3 ;
% v_core = 4000;
% v_corona = 4000;
n_agg = pi*R^2*L/v_corona;
beta_s = v_core * (sld_core - sld_solvent);
beta_c = v_corona * (sld_corona - sld_solvent);
x = (q.*Rg).^2;
Fc = 2 * (exp(-x) - 1 + x)./x.^2;
for i = 1:length(q)
psi(i) = (1 - exp(-q(i)*Rg))./(q(i)*Rg);
Fcs(i) = (2*besselj(1,q(i)*R)/(q(i)*R))^2;
Si(i) = integral(@(t) Fun_Si(t), 0, q(i)*L);
FL(i) = 2*Si(i)/(q(i)*L) - 4*(sin((q(i)*L)/2))^2./(q(i)*L)^2;
Fs(i) = Fcs(i)*FL(i);
Ssc(i) = psi(i)*2*besselj(1, q(i)*R)./(q(i)*R) * besselj(0, q (i)*(R+d*Rg)) * FL(i);
Scc(i) = psi(i)^2* (besselj(0, q(i)*(R+d*Rg)))^2 * FL(i);
end
term1 = n_agg.^2 * beta_s.^2* Fs;
term2 = n_agg * beta_c.^2.*Fc';
term3 = 2 * n_agg.^2 * beta_s * beta_c.* Ssc;
term4 = n_agg* (n_agg - 1) * beta_c.^2 * Scc;
P_worm = term1 + term2+ term3 + term4;
% loglog(q,term1)
% hold on
% loglog(q,term2)
% hold on
% loglog(q,term3)
% hold on
% loglog(q,term4)
% hold on
% legend('term1','term2','term3','term4')
I_worm = (A*P_worm + b)' * ndensity * 10^(-13);
end

function I_disk = coreshell_disk (q,param)
R = param(1);
Rg = param(2);
n_agg = param(3);
L = param(4);
d = param(5);
A = param(6);
b = param(7);
ndensity = param(8);

% SLD and v for F127
% sld_core = 0.344;
% sld_corona = 0.638;
% sld_solvent = 6.33;
% v_core = 6283;
% v_corona = 14667;
% For example testing
sld_core =  6.3 +1;
sld_corona = 6.3 +1;
sld_solvent = 6.3 ;
v_core = 4000;
v_corona = 4000;
beta_s = v_core * (sld_core - sld_solvent);
beta_c = v_corona * (sld_corona - sld_solvent);
for i = 1:length(q)
x(i) = (q(i)*Rg).^2;
Fc(i) = 2 * (exp(-x(i)) - 1 + x(i))/(x(i)^2);  
psi(i) = (1 - exp(-q(i)*Rg))./(q(i)*Rg);
Fcs (i) = (sin(q(i)*L/2)/(q(i)*L/2))^2;
FR(i) = 2 / (q(i)*R)^2 * (1 - (besselj(1, 2* q(i) * R)/(q(i)*R)));
Fs(i) = Fcs(i) * FR(i); 
Ssc(i) = psi(i) * sin(q(i)*L/2)/(q(i)*L/2) * cos(q(i) * (L/2 + d*Rg)) *FR(i);
Scc(i) = psi(i)^2 * (cos(q(i)*(L/2 + d*Rg)))^2 * FR(i);
end
% term 1 and term 4 are different
term1 = n_agg.^2 * beta_s.^2* Fs;
term2 = n_agg * beta_c.^2.*Fc;
term3 = 2 * n_agg.^2 * beta_s * beta_c.* Ssc;
term4 = n_agg* (n_agg - 1) * beta_c.^2 * Scc;
I_disk = term1 + term2+ term3 + term4;
% figure()
% loglog(q,term1)
% hold on
% loglog(q,term2)
% hold on
% loglog(q,term3)
% hold on
% loglog(q,term4)
% hold on
% legend('term1','term2','term3','term4')
end
% %% Core Shell sphere model
% function I = coreshell_sphere(q, A,b, rc, thickness)
% sld_core = 0.325;
% sld_shell = 0.547; 
% sld_solvent = 6.33;
% rs = rc + thickness;
% % Calculate the contribution of the core
% v_core = 4/3 * pi * rc.^3;
% core_qr = q.*rc;
% bess_core = 3 * (sin(core_qr) - core_qr* cos(core_qr))./(core_qr).^3;
% core_contrast = sld_core - sld_shell;
% F_core = v_core.*bess_core.*core_contrast;
% % calculate the contribution of the shell
% v_shell = 4/3 * pi * rs.^3;
% shell_qr = q.*rs;
% bess_shell = 3 * (sin(shell_qr) - shell_qr* cos(shell_qr))./(shell_qr).^3;
% shell_contrast = sld_shell - sld_solvent;
% F_shell = v_shell.*bess_shell.*shell_contrast;
% F = F_core + F_shell;
% I = A*F.^2./v_shell + b;
% end

function int_Fs_ellipsoid = Fsfun_ellipsoid(q, alpha, R, eps)
r = R *(sin(alpha).^2 + eps^2 * cos(alpha).^2)^(1/2);
phi = besselJ(q*r);
int_Fs_ellipsoid = phi.^2 * sin(alpha);
end

function int_Ssc = Ssc_fun_ellipsoid(q, alpha, R, Rg, d, eps)
r = R * (sin(alpha).^2 + eps.^2 * cos(alpha).^2).^(1/2);
phi = 3 * (sin(q.*r) - q.*r.* cos(q.* r))./(q.*r).^3;
int_Ssc = phi.*(sin(q*(r + d*Rg))./(q*(r + d*Rg))).* sin(alpha);
end

function int_Scc = Scc_fun_ellipsoid(q, R, eps,alpha, d, Rg)
r = R * (sin(alpha).^2 + eps.^2 * cos(alpha).^2).^(1/2);
int_Scc = (sin(q*(r + d*Rg))./(q*(r + d*Rg))).^2.* sin(alpha);
end

function J = besselJ(x)
J = 3 * (sin(x)./x - cos(x))./(x.^2);
end

function int_Fs_cylinder = Fsfun_cylinder(q, alpha, R, L)
Psi = 2*besselj(1,q*R*sin(alpha)).*sin(q*L*cos(alpha)/2)./ (q*R*sin(alpha).*q*L.*cos(alpha)./2);
int_Fs_cylinder = (Psi).^2 .* sin(alpha);
end

function int_Ssc_cylinder = Ssc_fun_cylinder(q, alpha, R, Rg, L, d)
Psi =  2*besselj(1,q*R*sin(alpha)).*sin(q*L*cos(alpha)/2)./ (q*R*sin(alpha).*q*L.*cos(alpha)./2);
xi_val = Xi_cylinder(q, R + d*Rg, L + 2*d*Rg, alpha);
int_Ssc_cylinder = Psi.* xi_val.* sin(alpha);
end

function int_Scc_cylinder = Scc_fun_cylinder(q, alpha, R, Rg, L, d)
xi_val = Xi_cylinder(q, R + d*Rg, L + 2*d*Rg, alpha);
int_Scc_cylinder = xi_val.^2.* sin(alpha);
end

function xi_val = Xi_cylinder(q, R, L, alpha)
xi_val = (R*2*besselj(1,q*R*sin(alpha)).*cos(q*L*cos(alpha)/2)./((R+L) * q*R.*sin(alpha)) ...
    + L*besselj(0,q*R*sin(alpha)).*sin(q*L*cos(alpha)./2)./((R+L)*q*L*cos(alpha)./2));
end

function Si_int = Fun_Si(t)
Si_int = sin(t)./t;
end

function I_fuzzy = fuzzy(q, param)
scale = param(1);
R = param(2);
sig = param(3);
sld_solvent = 6.33;
sld_sphere = param(4);
b = param(5);
V = 4/3 * pi* (R)^3;
%R = R + sig;
Phi = 3 * V * (sin(q.*R) - q.*R.* cos(q.* R))./(q.*R).^3;
A = Phi.*exp(-0.5*(sig*q).^2);
I_fuzzy = (scale/V * (sld_sphere - sld_solvent).^2 .* A.^2)*10^(-4) + b;
end

% Sphere model
function I_sphere = sphere(q, param)
scale = param(1);
R = param(2);
V = 4/3 * pi* R^3;
sld_solvent = 6.33;
sld_sphere = param(3);
b = param(4);
Phi = 3 * V * (sin(q.*R) - q.*R.* cos(q.* R))./(q.*R).^3;
rho = sld_sphere - sld_solvent;
I_sphere =  scale/V * rho.^2* Phi.^2 * 10^(-4) + b;
end

