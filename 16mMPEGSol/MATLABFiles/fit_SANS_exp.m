close all
clc
clear

%% Inputs
Nsample = 1;
paramsave = zeros(Nsample,16);
for k = 1:Nsample
filename = sprintf('16mMSol_0_merged.txt', k-1);
%N_input = 1000; % number of counts for the input file starting point
% a_sim = linspace(0.1,1,10); % number of counts to simulate
% a_sim(end+1) = 0.001;
% a_sim = logspace(-3,0,10);
%a_sim = 1;
a_sim = [2.5 5 10 25 50 100];
N_rep = 100;
saveout = 'y';
phi = 0.0685;

%% Generate simulated files
% Simulated files are named in this fashion: strcat(name,'_Nsim_',num2str(N_sim(i)),ext)
simulate_sigma(filename,a_sim);

%% Bootstrap from the simulated files
for i = 1:length(a_sim)
    [~,name,ext] = fileparts(filename);
    sim_file = strcat(name,ext);
    [I_rep,M_data] = bootstrap_SANS(sim_file, N_rep);
    if i == 1
        I_rep_mat = zeros(size(I_rep,1),size(I_rep,2),length(a_sim));
    end
    I_rep_mat(:,:,i) = I_rep;
end

%% Perform the fit 16 mM Gel
I_debye = zeros(size(I_rep_mat));
A = zeros(N_rep,length(a_sim)); Rg = zeros(N_rep,length(a_sim)); 
NSLDpfit = zeros(N_rep,length(a_sim)); bkg = zeros(N_rep,length(a_sim));
std_Rg = zeros(N_rep,length(a_sim)); std_NSLD = zeros(N_rep,length(a_sim)); 
std_bkg = zeros(N_rep,length(a_sim));
for i = 1:N_rep
    for j = 1:length(a_sim)
        sim_file = strcat(name, ext);
        M_data = importdata(sim_file);
        [A(i,j),Rg(i,j),NSLDpfit(i,j),bkg(i,j),std_Rg(i,j),std_NSLD(i,j),std_bkg(i,j)] = ...
            fit_GaussCoilfunc(M_data(:,1),I_rep_mat(:,i,j),M_data(:,3),phi);
        I_debye(:,i,j) = I_model(M_data(:,1),[Rg(i,j),A(i,j),bkg(i,j)]);
    end
end
%% Get the bootstrapped uncertainty
Rg_avg = mean(Rg); % [Angstrom]
NSLD_avg = mean(NSLDpfit); % [Angstrom^-2]
bkg_avg = mean(bkg); % [cm^-1]
Rg_bs_std = std(Rg);
NSLD_bs_std = std(NSLDpfit);
bkg_bs_std = std(bkg);

%% Get the average error from nlparci
std_Rg_avg = mean(std_Rg);
std_NSLD_avg = mean(std_NSLD);
std_bkg_avg = mean(std_bkg);

% %% Perform the fit-P4 to broad peak
% I_P4 = zeros(size(I_rep_mat));
% for i = 1:N_rep %N_rep
%      for j = 1:length(a_sim) %length(a_sim)
%          sim_file = strcat(name,ext);
%          M_data = importdata(sim_file);
%          [x_low, x_m, x_high, MSR, residual,ci] = correlation(M_data, I_rep_mat(:,i,j));
%          A (i,j)= x_low(1);
%          n_low (i,j) = x_low(2);
%          C (i,j) = x_m(3);
%          m (i,j) = x_m(4);
%          q0 (i,j) = x_m(5);
%          xi_corr(i,j) = x_m(6);
%          D(i,j) = x_high(1);
%          n_high (i,j) = x_high(2);
%          MSR_N(i,j) = MSR;
%      end
% end
% %% 
%  AIC(:,:,1) = length(residual)*log(MSR_N) + 2*6;
% %%
% A_mean = mean(A);
% n_low_mean = mean(n_low);
% C_mean = mean(C);
% m_mean = mean(m);
% q0_mean = mean(q0);
% xi_mean = mean(xi_corr);
% D_mean = mean(D);
% n_high_mean = mean(n_high);
% A_std = std(A);
% n_low_std = std(n_low);
% C_std = std(C);
% m_std = std(m);
% q0_std = std(q0);
% xi_std = std(xi_corr);
% D_std = std(D);
% n_high_std = std(n_high);
% %AIC_correlation_mean = mean(AIC_correlation);
% paramsave(k,1) = mean(A);
% paramsave(k,2) = mean(n_low);
% paramsave(k,3) = mean(C);
% paramsave(k,4) = mean(m);
% paramsave(k,5) = mean(q0);
% paramsave(k,6) = mean(xi_corr);
% paramsave(k,7) = mean(D);
% paramsave(k,8) = mean(n_high);
% paramsave(k,9) = std(A);
% paramsave(k,10) = std(n_low);
% paramsave(k,11) = std(C);
% paramsave(k,12) = std(m);
% paramsave(k,13) = std(q0);
% paramsave(k,14) = std(xi_corr);
% paramsave(k,15) = std(D);
% paramsave(k,16) = std(n_high);

% %% Fit to the glgel
% for i = 1:N_rep
%      for j = 1:length(a_sim)
%          sim_file = strcat(name,'_asim_',num2str(a_sim(j)),ext);
%          M_data = csvread(sim_file);
%          [Z, xi, MSR] = glgel(M_data, I_rep_mat(:,i,j));
%          Z_gl (i,j)= Z;
%          xi_gl(i,j) = xi;
%          MSR_glgel (i,j) = MSR/a_sim(j);
%      end
% end
%  AIC(:,:,2) = 2*MSR_glgel + 2 * 4;
%  
% %% Get the mean parameters from glgel
%  Z_gl_mean = mean (Z_gl);
%  xi_gl_mean = mean(xi_gl);
%  AIC_glgel_mean = mean (AIC_glgel);
%  
%  %% Fit to the fsgel
%  for i = 10
%      for j = 1
%          sim_file = strcat(name,ext);
%          M_data = csvread(sim_file);
%          [a1, a2, D, MSR] = fsgel(M_data, I_rep_mat(:,i,j));
%          a1 (i,j)= a1;
%          a2 (i,j) = a2;
%          D (i,j) = D;
%          MSR_fsgel (i,j) = MSR/a_sim(j);
%      end
%  end
% 
%  AIC(:,:,3) = 2*MSR_fsgel + 2 * 5;
%  
%  %% Get the mean parameters from fsgel
%  AIC_fsgel_mean = mean (AIC_fsgel);
%  
%  %% Fit to OZsqLorentz
%  for i = 1:N_Rep
%      for j = 1:length(a_sim)
%          sim_file = strcat(name,'_asim_',num2str(a_sim(j)),ext);
%          M_data = csvread(sim_file);
%          [xi_OZL, MSR] = OZsqLorentz(M_data, I_rep_mat);
%      end
%  end
%  
%   %% Fit to Debye Bueche
%  for i = 1
%      for j = 1
%          sim_file = strcat(name,'_asim_',num2str(a_sim(j)),ext);
%          M_data = csvread(sim_file);
%          [xi_DB, MSR] = DB(M_data, I_rep_mat);
%      end
%  end
%  
%  %% Fit to Ornstein-Zernike and pseudo-Voigt Lorentzian
%   for i = 1
%      for j = 1
%          sim_file = strcat(name,'_asim_',num2str(a_sim(j)),ext);
%          M_data = csvread(sim_file);
%          [xi_OZVL, MSR] = OZpVL(M_data, I_rep_mat);
%      end
%  end
%  %% Calculate the relative likelihood
%  % B_i=exp((AIC_min-AIC_i)/2); p_i = B_i/sum(B_i)
% 
% for i = 1:N_rep
%     for j = 1:length(a_sim)
%         for k = 1:3
%     B (i,j,k) = exp((min([AIC_correlation(i,j),AIC_fsgel(i,j),AIC_glgel(i,j)]) - AIC(i,j,k))/2);
%         end
%             corr_histogram(j) = sum(B(:,j,1));
%             glgel_histogram(j) = sum(B(:,j,2));
%             fsgel_histogram(j) = sum(B(:,j,3));
%     end
% end
% %%
% loglog(a_sim,corr_histogram,'ok','MarkerSize',8,'MarkerFaceColor','k')
% hold on
% loglog(a_sim,glgel_histogram,'or', 'MarkerSize',8,'MarkerFaceColor','r')
% hold on
% loglog(a_sim,fsgel_histogram,'ob', 'MarkerSize',8,'MarkerFaceColor','b')
% hold off
% set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2,'xscale','log');
% set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
% legend ('Broad Peak Model', 'Gaussian Lorentzian Gel', 'Fine-scale polymer gel')
% xlabel('Fraction of Counts')
% ylabel('Frequency of Best Fit')
% %%
% errorbar(a_sim,mean(MSR_correlation),std(MSR_correlation),'s','MarkerSize',8,'LineWidth',1,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k');
% hold on
% errorbar(a_sim,mean(MSR_glgel),std(MSR_glgel),'s','MarkerSize',8,'LineWidth',1,'Color','r','MarkerFaceColor','r','MarkerEdgeColor','r');
% hold on
% errorbar(a_sim,mean(MSR_fsgel),std(MSR_fsgel),'s','MarkerSize',8,'LineWidth',1,'Color','b','MarkerFaceColor','b','MarkerEdgeColor','b');
% set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2,'xscale','log','yscale','log');
% legend ('Broad Peak Model', 'Gaussian Lorentzian Gel', 'Fine-scale polymer gel')
% set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
% xlabel('Fraction of Counts')
% ylabel('Average Mean Square Residual')
% % loglog(a_sim,MSR_correlation,'ok','MarkerSize',8,'MarkerFaceColor','k')
% % hold on
% % loglog(a_sim,MSR_glgel,'or','MarkerSize',8,'MarkerFaceColor','r')
% % hold on
% % loglog(a_sim,MSR_fsgel,'ob','MarkerSize',8,'MarkerFaceColor','b')
% % hold off
% % legend ('Broad Peak Model', 'Gaussian Lorentzian Gel', 'Fine-scale polymer gel')
% % set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2,'xscale','log');
% % set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
% % legend ('Broad Peak Model', 'Gaussian Lorentzian Gel', 'Fine-scale polymer gel')
% % xlabel('Fraction of Counts')
% % ylabel('Average Mean Square Residual')

% 
%% Plot the fit
figure()

for i = length(a_sim)
    sim_file = strcat(name,'_asim_',num2str(a_sim(i)),ext);
    M_data = csvread(sim_file);
%     subplot(2,5,i)
%     box on;
    hold on;
    errorbar (M_data(:,1)*10,I_rep_mat(:,2,i),M_data(:,3),'o','LineWidth',3,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
    hold on
    loglog(M_data(:,1)*10,I_debye(:,2,i),'Color','r','LineWidth',2)
    xlabel("q [A^-^1]",'FontWeight','bold');
    ylabel('I(q) [cm^-^1]','FontWeight','bold');
    set(gca,'FontSize',20,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
    set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
    xlim([10^(-1) 4.3])
    ylim([0.0001 10]);
%     set(gca,'FontSize',14,'TickLength',[0.03 0.03],'LineWidth',2,'yscale','log','xscale','log');
%     set(gcf,'Color','w','units','normalized','outerposition',[0 0 1 1]);
    %saveas(gcf,'Iq_bootstrapped.png');
end
% 
%% Plot the parameters as a function of count
% % Rg
% % Calculate what 5% error means for Rg
errorbound = 0.05*(Rg_avg(1,length(Rg_avg))/10); % [nm]
figure()
box on;
hold on;
errorbar(a_sim,Rg_avg/10,Rg_bs_std/10,'s','MarkerSize',12,'LineWidth',2,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k');
plot(a_sim,(Rg_avg(1,length(Rg_avg))/10 - errorbound)*ones(length(a_sim)),'--','LineWidth',2,'Color','b');
plot(a_sim,(Rg_avg(1,length(Rg_avg))/10 + errorbound),'--','LineWidth',2,'Color','b');
%xlim([0.0009 1.1])
xlabel('Fraction of Counts','FontWeight','bold');
ylabel('R_g [nm]','FontWeight','bold');
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2,'xscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
%saveas(gcf,'Rg with count.png');
end
%%
%csvwrite('Bootstrap_16mMSol.csv',[A,Rg,NSLDpfit,bkg,std_Rg,std_NSLD,std_bkg])
%% Fitting functions
function [A,Rg,NSLDpfit,bkgfit,std_Rg,std_NSLD,std_bkg] = fit_GaussCoilfunc(q,I_expt,SD,phi)

%% Parameters for Fit
% Fitting parameters
% x = [Rg, A (prefactor), bkg]
% A = phi*d_SLD^2*N*v
[v,N] = mono_volume('PNIPAM',28050);
A0 = phi*2.68^2*N*v*1e8; %convert from 1/A to 1/cm because bkg is in 1/cm
% Rg is in A, which is good because q is in 1/A
lb = [10; 0; 0]; %lower bounds
ub = [80; 10; 10]; %upper bounds
x0 = [30; A0; 1]; %guesses
% Overall final Ifit(q) is in 1/cm, better scaling for fit
[x_out,ci] = Debyefit(q,I_expt,SD,lb,ub,x0);

%% Extract physical quantities
% Radius of gyration
Rg = x_out(1); %Angstrom
% Prefactor
A = x_out(2);
% Scattering length density of polymer
NSLDd = 3.334;
NSLDpfit = -sqrt((x_out(2)*1e-8)/(phi*N*v)) + NSLDd; %1/A^2
% Incoherent background scattering
bkgfit = x_out(3); %1/cm
% Confidence intervals for all fit parameters
ci_Rg = x_out(1)-ci(1,1);
ci_A = x_out(2)-ci(2,1);
ci_bkg = x_out(3)-ci(3,1);
% Standard deviation of each fit parameter
std_Rg = ci_Rg/1.96; %convert to standard deviation
std_A = ci_A/1.96; %convert to standard deviation
std_bkg = (ci_bkg/1.96); %convert to standard deviation
std_NSLD = abs((NSLDpfit*0.5*std_A)/x_out(2)); %the unit conversion canceled out
end


function [x_out,ci] = Debyefit(q,I_expt,SD,lb,ub,x0)
%% Set bounds
%Rg, A, bkg is the order
%Rg in Angstrom
%A in 1/cm
%bkg in 1/cm
% lb lower bounds
% ub upper bounds
% x0 guesses

%% Define the fitting problem
problem.options=optimoptions('lsqnonlin','MaxFunEvals',inf,'MaxIter',inf,...
    'TolFun',1e-12,'TolX',1e-12,'display','off');
problem.solver='lsqnonlin';
problem.objective=@(x)LSfunc(I_expt,SD,q,x);
problem.lb=lb;
problem.ub=ub;
problem.x0=x0;

%% Obtain the fit
% x_out = [Rg, A, bkg] in [Angstrom,1/cm,1/cm]
% ci = confidence intervals in [Rg, A, bkg]
[x_out,~,residuals,~,~,~,jacobians] = lsqnonlin(problem);
ci = nlparci(x_out,residuals,'jacobian',jacobians);
end

function LS = LSfunc(I_expt,SD,q,x)
%Obtain fit I data from GaussCoil (Gaussian chain form factor)
I_fit = I_model(q,x);
%Calculate the residual
diff = I_fit - I_expt;
%Create the LS objective function
LS = diff./SD;
end

function I_fit = I_model(q,x)
Rg = x(1); %radius of gyration
A = x(2); %prefactor
b = x(3); %incoherent background
P = PGaussCoil(q,Rg); %form factor
I_fit = A*P + b;
end

function P = PGaussCoil(q,Rg)
x = q.^2.*Rg^2;
P = (2.*(exp(-x)+x-1))./(x.^2);
end

%% MODEL: Correlation peak
% I(q) = A/q^n + C/(1+(|q-q0|xi)^m + B [no B, already taken care of]
% Guess in this order: [A, n, C, m, q0, xi]
function I_fit = CorrLength(q,params)
A = params(1); n = params(2); 
C = params(3); m = params(4);
q0 = params(5); xi = params(6);
I_fit = A./q.^n + C./(1+(abs(q-q0)*xi).^m);
end

%% MODEL: Correlation peak: low/high q clustering 
function I_fit = CorrLengthLowQ(q,params)
A = params(1); n = params(2);
I_fit = A./q.^n;
end

%% MODEL: Correlation peak: mid/high q solvation
function I_fit = CorrLengthMidQ(q,params)
C = params(1); m = params(2);
q0 = params(3); xi = params(4);
I_fit = C./(1+(abs(q-q0)*xi).^m);
end