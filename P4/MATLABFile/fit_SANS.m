close all
clc
clear

%% Inputs
k = 1; % change to 1 for bootstrapping
Nsample = 1;

% we need to change folder for every fraction of counts to get the
% bootstrapped parameters. Remember to change the filename when you do
% that.
filename = sprintf('P4_%d_merged.txt', k-1);
savefilename_1 = 'Bootstrap_P4_1_piecefit.csv';
savefilename_2 = 'Bootstrap_P4_1_average_piecefit.csv';
%N_input = 1000; % number of counts for the input file starting point
% a_sim = linspace(0.1,1,10); % number of counts to simulate
% a_sim(end+1) = 0.001;
% a_sim = logspace(-3,0,10);
%a_sim = [1 2.5 5 10 25 50 100];
a_sim = 1;
%a_sim = [0.001 0.0025 0.005 0.01 0.025 0.05 0.1 0.25 0.5 1];
N_rep = 100;
saveout = 'y';
phi = 0.05;

%% Generate simulated files
% Simulated files are named in this fashion: strcat(name,'_Nsim_',num2str(N_sim(i)),ext)
%simulate_sigma(filename,a_sim);

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
%% plot the bootstrap data file to compare with experimental file (original)
% plot only the data used for fitting
% qcut = [0.005 0.6]; 
% q = M_data(M_data(:,1)>qcut(1) & M_data(:,1)<qcut(2),1);
% p4_cut = I_rep_mat(M_data(:,1)>qcut(1) & M_data(:,1)<qcut(2),:);
% p4_cut_error = I_rep_mat(M_data(:,1)>qcut(1) & M_data(:,1)<qcut(2),3);
% figure()
% % this is plotting the bootstrap data
% plot(q, p4_cut(:,2:end),'s')
% hold on
% % this is plotting the original experimental file used for bootstrapping
% errorbar(q, p4_cut(:,1), p4_cut_error,'.','MarkerEdgeColor','k','MarkerSize',10,'Color','k')
% xlabel("q [nm^-^1]",'FontWeight','bold');
% ylabel('I(q) [cm^-^1]','FontWeight','bold');
% set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',1,'xscale','log');
% set(gcf,'Color','w','units','pixels','outerposition',[200 100 600 600]);


%% Perform the fit-P4 to broad peak
I_P4 = zeros(size(I_rep_mat));
for i = 1:N_rep %N_rep
     for j = 1:length(a_sim) %length(a_sim)
         sim_file = strcat(name,ext);
         M_data = importdata(sim_file);
         [params, I_out] = correlation(M_data, I_rep_mat(:,i,j));
         A (i,j)= params.x_low(1);
         n_low (i,j) = params.x_low(2);
         C (i,j) = params.x(3);
         m (i,j) = params.x(4);
         q0 (i,j) = params.x(5);
         xi_corr(i,j) = params.x(6);
         D(i,j) = params.x_high(1);
         n_high (i,j) = params.x_high(2);
         MSR_N(i,j) = params.MSR;
         I_fit(i,:) = I_out.I_fit;
         I_fit_lowq (i,:)= I_out.I_fit_lowq';
         I_fit_midq(i,:) = I_out.I_fit_midq';
         I_fit_highq (i,:)= I_out.I_fit_highq' ;
         bkg(i,:) = I_out.I_fit_bkg';
         I_sub_solv(i,:) = I_out.I_sub_solv ;
         dI_sub_solv(i,:) = I_out.dI_sub_solv;
         ci(i,:) = params.error_ci;
         residual(i,:) = params.residual;
     end
end
%% 
 AIC(:,:,1) = length(params.residual)*log(MSR_N) + 2*6;
 BIC(:,:,1) = length(params.residual)*log(MSR_N) + 6 * log(length(params.residual));
%% 
A_mean = mean(A);
n_low_mean = mean(n_low);
C_mean = mean(C);
m_mean = mean(m);
q0_mean = mean(q0);
xi_mean = mean(xi_corr);
D_mean = mean(D);
n_high_mean = mean(n_high);
A_std = std(A);
n_low_std = std(n_low);
C_std = std(C);
m_std = std(m);
q0_std = std(q0);
xi_std = std(xi_corr);
D_std = std(D);
n_high_std = std(n_high);
% x = [A, n, C, m, q0, xi]
if isfile(savefilename_1)
     disp('File exists! Please rename the filename')
else
     T = table(A, n_low, C, m, q0, xi_corr, D, n_high, AIC, BIC, ci);
     T_average = table(A_mean, n_low_mean, C_mean, m_mean, q0_mean, xi_mean, D_mean, n_high_mean, mean(AIC), mean(BIC), A_std, n_low_std, C_std, m_std, q0_std, xi_std, D_std, n_high_std);
     writetable(T, savefilename_1);
     writetable(T_average, savefilename_2);
end

%% Plot the 100 replicates fit on the data
figure()
for i = 1:length(a_sim)
    figure(1)
    hold on
    box on
    
    set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
    %set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log');
    set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
    errorbar(I_out.qfit*10,I_sub_solv(1,:),dI_sub_solv(1,:),'s','MarkerSize',10,'MarkerEdgeColor','k','LineWidth',1,'Color','k');
    plot(logspace(-3,0)*10,I_fit(1, :)+bkg(1, :),'Color','r','LineWidth',2);
    plot(logspace(-3,0)*10,bkg,'Color','b','LineWidth',2);
    xlabel("q [nm^-^1]",'FontWeight','bold');
    ylabel('I(q) [cm^-^1]','FontWeight','bold');
    xlim(I_out.qcut*10);
    ylim([0.01 10]);
    figure()
    hold on
    box on
    
    set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
    %set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log');
    set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
    plot(I_out.qfit*10,I_sub_solv(2:end,:),'s','MarkerSize',10);
    errorbar(I_out.qfit*10,I_sub_solv(1,:),dI_sub_solv(1,:),'.','MarkerSize',10,'MarkerEdgeColor','k','LineWidth',2,'Color','k');
    plot(logspace(-3,0)*10,I_fit(1, :)+bkg(1, :),'Color','r','LineWidth',2);
    
    
    plot(logspace(-3,0)*10,bkg,'Color','b','LineWidth',2);

    plot(logspace(-3,0)*10,I_fit(1, :)+bkg(1, :),'Color','r','LineWidth',2);

    plot(logspace(-3,0)*10,I_fit(2:end, :)+bkg(2:end, :),'Color','g','LineWidth',2);
    %legend('Expt','Bkg','Fit','Location','Northeast');
    %     saveas(gcf,sprintf('P4_corr_fit_%d.png',i-1))

    plot(logspace(-3,0)*10,I_fit(1, :)+bkg(1, :),'Color','r','LineWidth',2);
    xlabel("q [nm^-^1]",'FontWeight','bold');
    ylabel('I(q) [cm^-^1]','FontWeight','bold');
    xlim(I_out.qcut*10);
    ylim([0.01 10]);
    %plot(logspace(-3,0)*10,I_fit_lowq+bkg,'--','Color','b','LineWidth',2);
    %legend('Expt','Bkg','Fit','LowQ Fit','Location','Northeast');
    %     saveas(gcf,'P4_lowq.png')

    %plot(logspace(-3,0)*10,I_fit_midq+bkg,'.','Color','b','MarkerSize',10);
    %legend('Expt','Bkg','Fit','LowQ Fit','MidQ Fit','Location','Northeast');
    %    saveas(gcf,'P4_midq.png')

    %plot(logspace(-3,0)*10,I_fit_highq(1,:)+bkg(1,:),'-.','Color','b','LineWidth',2);
    %legend('Expt','Bkg','Fit','LowQ Fit','MidQ Fit','HighQFit','Location','Northeast');
end
%% checking normality 
tiledlayout(1,2);
% histogram of experimental data
nexttile
histogram(residual(1,:),'LineWidth',2)
legend('Experimental Data')
xlabel('Weighted Residual','FontWeight','bold')
ylabel('Frequency','FontWeight','bold')
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
nexttile
histogram(residual(randi(100,1),:),'LineWidth',2, 'FaceColor', "#D95319")
legend('Bootstrapping Data')
xlabel('Weighted Residual','FontWeight','bold')
ylabel('Frequency','FontWeight','bold')

set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);

% calculate the skewness
data_skew_calc = skewness(residual(1,:));
bootstrap_skew_calc = skewness(residual(randi(100,1),:));
%% Fitting functions
function [A,Rg,NSLDpfit,bkgfit,std_Rg,std_NSLD,std_bkg] = fit_GaussCoilfunc(q,I_expt,SD,phi)
%% Calculate physical parameters
NSLDp = pureNSLD('PNIPAM'); %1/A^2
NSLDd = pureNSLD('D2O'); %1/A^2

%% Parameters for Fit
% Fitting parameters
% x = [Rg, A (prefactor), bkg]
% A = phi*d_SLD^2*N*v
[v,N] = mono_volume('PNIPAM',28050);
A0 = phi*(NSLDp-NSLDd)^2*N*v*1e8; %convert from 1/A to 1/cm because bkg is in 1/cm
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