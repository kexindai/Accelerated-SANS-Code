%% Data reduction for 16mMSol Rg ~ 26A
% Prior to running this MATLAB code, the user should run the
% stitch_16mMSol_new.m file to merge the large q and small q ranges
clc
clear
close all
% volume fraction of the PEG polymer in deuterated DMF
vp = 0.0685;
% the number is the number of samples in each folder of
% CD_2configTimeSlicing_16mMSol_500pts_X_merged, where X can vary from
% 0.01, 0.025, 0.05, 0.10, 0.25, 0.5, 1
number = 1; % I only used the first dataset but it is better for the code to fit all of the dataset in the folder
for i = 1:number
%clearvars I_debye
% read the merged file
q = importdata(sprintf('16mMSol_%d_merged.txt',i-1));
% perform the curve fitting and get confidence intervals and parameters
[A(i),Rg(i),bkg(i),ci_A(i),ci_Rg(i),ci_bkg(i)] = ...
    fit_GaussCoilfunc(q(:,1),q(:,2),q(:,3),vp);
I_debye = I_model(q(:,1),[Rg(i),A(i),bkg(i)]);
figure() 
% plot the data
errorbar(q(:,1)*10,q(:,2),q(:,3),'o','LineWidth',3,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
hold on
% Plot the fitting
loglog(q(:,1)*10,I_debye,'Color','r','LineWidth', 4)
xlim([10^(-1) 4.3])
ylim([10^(-4) 1.5])
xlabel('q [nm^-^1]','FontWeight','bold');
ylabel('I(q) [cm^-^1]','FontWeight','bold');
% save the plots
%saveas(gcf,sprintf('Debyefit_%d.png',i-1));
end
%% save the parameters
T = table(A',Rg',bkg',ci_A',ci_Rg',ci_bkg', 'VariableNames',{'A','Rg','bkg','ci_A','ci_Rg','ci_bkg'});
%writetable(T, 'Debye_paramsexp_0p5.csv')

%% this section tests for the correlation between errors in 16mM PEG solution
% Work flow
% 1. Take the model and the fitted parameters f(q, theta) where theta is
% the parmaters
% 2. Fix theta and add noise to the model. We will first test Gaussian
% noise (the assumption for MC bootstrapping)
% 3. We will create 100 replicates and fit the model with broad peak model
% 4. If nonlinear least squares is a biased estimator (errors are correlated), 
% then your theta’ will have expectation something other than theta

% the noise is gaussian centered at 0 but porportional to the dI
Num_rep = 100;
noise = q(:,3).*randn(length(q), Num_rep);
I_rep_data = repmat(I_debye, [1, Num_rep]) + noise;

% this plots the original data
figure()
loglog(q(:,1), q(:,2),'o')
hold on
loglog(q(:,1), I_rep_data + bkg,'s')
hold on
loglog(q(:,1), I_debye + bkg,'-','color','g','LineWidth',2)
hold on
%plot(q(:,1), ones(length(q))* bkg,'Color', 'b', 'LineWidth',2);
legend('original data','simulated data with noise', 'original data fit')
xlim([0 0.43]);
ylim([0.001 100]);
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
xlabel('q [nm^-^1]','FontWeight','bold');
ylabel('I(q) [cm^-^1]','FontWeight','bold');

% now we fit the simulated data again
A_test = zeros(Num_rep,1);
Rg_test = zeros(Num_rep,1);
bkg_test = zeros(Num_rep,1);
ci_A_test = zeros(Num_rep,1);
ci_Rg_test = zeros(Num_rep,1);
ci_bkg_test = zeros(Num_rep,1);
I_debye_test = zeros(Num_rep,length(q));
for i = 1:Num_rep
%clearvars I_debye
% read the merged file

% perform the curve fitting and get confidence intervals and parameters
[A_test(i),Rg_test(i),bkg_test(i),ci_A_test(i),ci_Rg_test(i),ci_bkg_test(i)] = ...
    fit_GaussCoilfunc(q(:,1),I_rep_data(:,i),q(:,3),vp);
I_debye_test(i,:) = I_model(q(:,1),[Rg_test(i),A_test(i),bkg_test(i)]); 

end
%% Calculate the expectation value from the simulation
E_A = mean(A_test);
E_Rg = mean(Rg_test);
E_bkg = mean(bkg_test);

% plot the histogram
figure()
histogram(A_test,'LineWidth',2);
xlabel('A values')
ylabel('Frequency')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);

figure()
histogram(Rg_test,'LineWidth',2);
xlabel('Rg values')
ylabel('Frequency')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);


%% plot the simulated data and the fit
figure()
loglog(q(:,1)*10,I_rep_data,'s')
set(gca,'FontSize',20,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
hold on
% Plot the fitting
loglog(q(:,1)*10,I_debye_test,'Color','r','LineWidth',2)
hold on
%loglog(q(:,1)*10, ones(length(q))* bkg,'b','LineWidth',3)
xlim([10^(-1) 4.3])
ylim([10^(-4) 1.5])
xlabel('q [nm^-^1]','FontWeight','bold');
ylabel('I(q) [cm^-^1]','FontWeight','bold');
% save the plots
%saveas(gcf,sprintf('Debyefit_%d.png',i-1));
%% Fitting functions
function [A,Rg,bkgfit,ci_A,ci_Rg,ci_bkg] = fit_GaussCoilfunc(q,I_expt,SD,phi)
%% Calculate physical parameters
%% Parameters for Fit
% Fitting parameters
% x = [Rg, A (prefactor), bkg]
% A = phi*d_SLD^2*N*v
% this is just initial guess so it is fine to leave it
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
x0 = [30; A0; 1]; %guesses/cm, better scaling for fit
[x_out,ci] = Debyefit(q,I_expt,SD,lb,ub,x0);

%% Extract physical quantities
% Radius of gyration
Rg = x_out(1); %Angstrom
% Prefactor
A = x_out(2);
% Scattering length density of polymer
%NSLDpfit = -sqrt((x_out(2)*1e-8)/(phi*N*v)) + NSLDd; %1/A^2
% Incoherent background scattering
bkgfit = x_out(3); %1/cm
% Confidence intervals for all fit parameters
ci_Rg = x_out(1)-ci(1,1);
ci_A = x_out(2)-ci(2,1);
ci_bkg = x_out(3)-ci(3,1);
% % Standard deviation of each fit parameter
% std_Rg = ci_Rg/1.96; %convert to standard deviation
% std_A = ci_A/1.96; %convert to standard deviation
% std_bkg = (ci_bkg/1.96); %convert to standard deviation
% std_NSLD = abs((NSLDpfit*0.5*std_A)/x_out(2)); %the unit conversion canceled out
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
