close all
clc
clear

%% Inputs
% change folder to get bootstrapping result from reduced number of counts
filename = '16mMSol_0_merged.txt';
%a_sim = [0.001 0.0025 0.005 0.01 0.025 0.05 0.1 0.25 0.5 1];
% this is to get the bootstrapping error from reduced number of counts
a_sim = 1;
N_rep = 100;
saveout = 'y';
phi = 0.0685;

%% Generate simulated files
% Simulated files are named in this fashion: strcat(name,'_Nsim_',num2str(N_sim(i)),ext)
simulate_sigma(filename,a_sim);

%% Bootstrap from the simulated files
for i = 1:length(a_sim)
    [~,name,ext] = fileparts(filename);
    sim_file = strcat(name,'_asim_',num2str(a_sim(i)),ext);
    [I_rep,M_data] = bootstrap_SANS(sim_file,N_rep);
    if i == 1
        I_rep_mat = zeros(size(I_rep,1),size(I_rep,2),length(a_sim));
    end
    I_rep_mat(:,:,i) = I_rep;
end

%% Perform the fit
I_debye = zeros(size(I_rep_mat));
A = zeros(N_rep,length(a_sim)); Rg = zeros(N_rep,length(a_sim)); 
NSLDpfit = zeros(N_rep,length(a_sim)); bkg = zeros(N_rep,length(a_sim));
std_Rg = zeros(N_rep,length(a_sim)); std_NSLD = zeros(N_rep,length(a_sim)); 
std_bkg = zeros(N_rep,length(a_sim));
for i = 1:N_rep
    for j = 1:length(a_sim)
        sim_file = strcat(name,'_asim_',num2str(a_sim(j)),ext);
        M_data = csvread(sim_file);
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
% record Rg_avg and Rg_bs)std manually to Rgvsfrac.csv
%% Get the average error from nlparci
std_Rg_avg = mean(std_Rg);
std_NSLD_avg = mean(std_NSLD);
std_bkg_avg = mean(std_bkg);

% Plot the fit
figure()
for i = 1:length(a_sim)
    sim_file = strcat(name,'_asim_',num2str(a_sim(j)),ext);
    M_data = csvread(sim_file);
    subplot(2,5,i)
    box on;
    hold on;
    errorbar(M_data(:,1)*10,I_rep_mat(:,2,i),M_data(:,3),'','LineWidth',2,'Color','k');
    %plot(repmat(M_data(:,1)*10,1,N_rep),I_rep_mat(:,:,i),'s','LineWidth',2,'Color','k');
    plot(M_data(:,1)*10,I_debye(:,:,i),'LineWidth',2,'Color','r');
    xlabel('q [nm^-^1]','FontWeight','bold');
    ylabel('I(q) [cm^-^1]','FontWeight','bold');
    ylim([0.001 10]);
    set(gca,'FontSize',14,'TickLength',[0.03 0.03],'LineWidth',2,'yscale','log','xscale','log');
    set(gcf,'Color','w','units','normalized','outerposition',[0 0 1 1]);
    saveas(gcf,'16mMSol Iq_bootstrapped.png');
end

%% Plot the parameters as a function of count
% Rg
% Calculate what 5% error means for Rg
errorbound = 0.05*(Rg_avg(1,10)/10); % [nm]
figure()
box on;
hold on;
errorbar(a_sim,Rg_avg/10,Rg_bs_std/10,'s','MarkerSize',12,'LineWidth',2,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k');
plot(a_sim,ones(length(a_sim))*(Rg_avg(1,10)/10 - errorbound),'--','LineWidth',2,'Color','b')
plot(a_sim,ones(length(a_sim))*(Rg_avg(1,10)/10 + errorbound),'--','LineWidth',2,'Color','b')
xlim([0.0009 1.1])
xlabel('Fraction of Counts','FontWeight','bold');
ylabel('R_g [nm]','FontWeight','bold');
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2,'xscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
saveas(gcf,'Rg with count.png');
% %% Charlotte new section for scaling analysis
% hold off
% loglog(a_sim,2*Rg_bs_std/10,'o','LineWidth',2,'MarkerFaceColor','k','MarkerEdgeColor','k')
% hold on
% for i = 1:length(a_sim)
%     y_fit(i) = 0.1241*a_sim(i)^(-0.6109);
% end
% loglog(a_sim,2*y_fit/10,'LineWidth',2)
% xlabel('Fraction of Counts','FontWeight','bold');
% ylabel('Parameter Uncertainty','FontWeight','bold');
% set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2,'xscale','log');
% %set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
% hold off
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