function [Z, xi , MSR] = glgel(M_data, I_rep_mat)
%% Input
qcut = [0.005 0.6]; % [A^-1] min and max q to keep (from 0.005 to 0.6 look good)
lowq = [0.005 0.011]; % [A^-1] very low q for the "clustering" region
highq = [0.24 0.36]; % [A^-1] high q for the region where scaling changes 
vhighq = [0.36 0.6]; % [A^-1] very high q for the background (flattened)
% Which qs to fit? 
% For Broad-peak Correlation Length model, use 
fitq = [0.005 0.6]; % [A^-1] q for fitting the "solvation" region (optional to use 0.6 as highest q)
vp = 0.05; % volume fraction of P4 (corresponds to 6.5%--see labnotebook)
model = 'glgel'; % correlation,glgel,fsgel

% %% Create data saving directories
% % Create a folder to hold the data so that never again do we have to
% % overwrite data by accident on different days
% hirespath = fullfile(pwd,'HiRes');
% if ~exist(hirespath, 'dir')
%     mkdir(hirespath);
% end

%% Read in the data
% Skip 2 rows but read all columns
 
% p4 = csvread('6p5_P4_stitched_sub.csv',2,0);
 buffer = csvread('buffer_stitched.txt',2,0);

% Check that the data is correct by plotting it
% figure()
% box on;
% hold on;
% errorbar(p4(:,1),p4(:,2),p4(:,3),'s','MarkerSize',6,'MarkerEdgeColor','k','LineWidth',1,'Color','k');
% errorbar(buffer(:,1),buffer(:,2),buffer(:,3),'.','MarkerSize',20,'Color','b');
% xlabel("q [nm^-^1]",'FontWeight','bold');
% ylabel('I(q) [cm^-^1]','FontWeight','bold');
% xlim(qcut);
% ylim([0.04 10]);
% title('6.5% P4: 1D SANS with buffer')
% legend('P4','100 mM Phosphate Buffer','Location','best');
% set(gcf,'Color','w','units','pixels','outerposition',[200 100 600 600]);
% set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',1,'xscale','log',...
%     'yscale','log');

% % Kratky-Porod Plot
% figure()
% box on;
% hold on;
% plot(p4(:,1)*10,(p4(:,1)*10).^2.*p4(:,2)/1e7,'s','MarkerSize',6,'MarkerEdgeColor','k','LineWidth',1,'Color','k');
% xlabel(strcat("q [nm^-^1]"),'FontWeight','bold');
% ylabel('q^2I(q) [nm^-^3]','FontWeight','bold');
% xlim(qcut*10);
% set(gcf,'Color','w','units','pixels','outerposition',[200 100 600 600]);
% set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',1,'xscale','log',...
%     'yscale','log');
% 
% % Orstein-Zernike Plot (1/I vs. Q^2)
% figure()
% box on;
% hold on;
% plot((p4(:,1)*10).^2,1./p4(:,2),'s','MarkerSize',6,'MarkerEdgeColor','k','LineWidth',1,'Color','k');
% xlabel('q^2 [nm^-^2]','FontWeight','bold');
% ylabel('I(q)^-^1 [cm]','FontWeight','bold');
% % xlim((qcut*10).^2);
% xlim([0.05 2.4].^2);
% set(gcf,'Color','w','units','pixels','outerposition',[200 100 600 600]);
% set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',1);

%% Truncate the data by getting rid of bad points at the edges of the plot
% (using the qcut input)
% Sample
p4_cut(:,1) = M_data(M_data(:,1)>qcut(1) & M_data(:,1)<qcut(2),1); % q [A^-1]
q = p4_cut(:,1);
p4_cut(:,2) = I_rep_mat(M_data(:,1)>qcut(1) & M_data(:,1)<qcut(2),1); % I [cm^-1]
p4_cut(:,3) = M_data(M_data(:,1)>qcut(1) & M_data(:,1)<qcut(2),3); % dI [cm^-1]
p4_cut(:,4) = M_data(M_data(:,1)>qcut(1) & M_data(:,1)<qcut(2),4); % dq [A^-1]
% Buffer
buffer_cut(:,1) = buffer(M_data(:,1)>qcut(1) & M_data(:,1)<qcut(2),1); % q [A^-1]
buffer_cut(:,2) = buffer(M_data(:,1)>qcut(1) & M_data(:,1)<qcut(2),2); % I [cm^-1]
buffer_cut(:,3) = buffer(M_data(:,1)>qcut(1) & M_data(:,1)<qcut(2),3); % dI [cm^-1]
buffer_cut(:,4) = buffer(M_data(:,1)>qcut(1) & M_data(:,1)<qcut(2),4); % dq [A^-1]

%% Background subtraction: solvent and incoherence
% Subtract solvent background, scaled by vol frac
I_sub_solv = p4_cut(:,2) - (1-vp)*buffer_cut(:,2); 
% Propagate the error from background subtraction
dI_sub_solv = sqrt(p4_cut(:,3).^2 + buffer_cut(:,3).^2);
% Subtract out the incoherent background (taken such that the high q goes
% to 0), then propagate error again.
% First truncate data to the very high q region
q_vhigh = q(q > vhighq(1) & q < vhighq(2));
I_sub_vhigh = I_sub_solv(q > vhighq(1) & q < vhighq(2));
% Take a flat average and get the standard deviation
bkg_avg = mean(I_sub_vhigh); 
bkg_std = std(I_sub_vhigh);
% Subtract the incoherent background
I_sub = I_sub_solv - bkg_avg*ones(size(I_sub_solv)); 
% Propagate the error
dI_sub = sqrt(dI_sub_solv.^2 + bkg_std^2);

%% Check subtracted curve
% figure()
% box on;
% hold on;
% errorbar(q,I_sub,dI_sub,'s','MarkerSize',6,'MarkerEdgeColor','k','LineWidth',1,'Color','k');
% xlabel("q [nm^-^1]",'FontWeight','bold');
% ylabel('I(q) [cm^-^1]','FontWeight','bold');
% title('6.5% P4: 1D SANS, subtracted')
% xlim(qcut);
% ylim([0.001 10]);
% set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',1,'xscale','log',...
%     'yscale','log');
% set(gcf,'Color','w','units','pixels','outerposition',[200 100 600 600]);
% 
% csvwrite('6p5_P4_stitched_sub.csv',[p4_cut(:,1),I_sub,dI_sub,p4_cut(:,4)]);

%% Further truncate for low q and high q
q_low = q(q > lowq(1) & q < lowq(2));
I_sub_low = I_sub(q > lowq(1) & q < lowq(2));
dI_sub_low = dI_sub(q > lowq(1) & q < lowq(2));

q_fit = q(q > fitq(1) & q < fitq(2));
I_sub_fit = I_sub(q > fitq(1) & q < fitq(2));
dI_sub_fit = dI_sub(q > fitq(1) & q < fitq(2));

q_high = q(q > highq(1) & q < highq(2));
I_sub_high = I_sub(q > highq(1) & q < highq(2));
dI_sub_high = dI_sub(q > highq(1) & q < highq(2));

%% Perform fit
if strcmp(model,'correlation')
    % Fit to correlation peak model: 
    % I(q) = A/q^n + C/(1+(|q-q0|xi)^m) + B [we do not fit B--we set it]
    % Guess in this order: [A, n, C, m, q0, xi]
    % B is going to be set such that the SANS data go to 0 at high q. [after
    % solvent subtraction!!]
    % For piece-wise fitting, we have low q as [A,n] and high q as [C,m,q0,xi]
    %%%%%%%%%%%%%%%%%% For FULL Correlation Peak Model %%%%%%%%%%%%%%%%%%%%%%%%
    % x0 = [1e-5; 2; 10; 2; 0.05; 22]; % [Angstrom]
    % lb = [0; 0; 0; 0; 0.01; 0]; 
    % ub = [0.1; 5; inf; 5; 0.07; 100];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % [A, n]
    x0_l = [1e-5; 2]; % [Inverse Angstrom q and cm I]
    lb_l = [0; 0]; 
    ub_l = [0.1; 5];
    % [A, n]
    x0_h = [1e-5; 4]; % [Inverse Angstrom q and cm I]
    lb_h = [0; 0]; 
    ub_h = [0.1; 10];
    options=optimoptions('lsqnonlin','MaxFunEvals',inf,'MaxIter',inf,...
        'TolFun',1e-12,'TolX',1e-12);

    % Fit the low q region
    [x_low,resnorml,residuall,~,~,~,jacobianl] = ...
        lsqnonlin(@(x) LSfunc(@(q_low,x) CorrLengthLowQ(q_low,x),q_low,I_sub_low,dI_sub_low,x),x0_l,lb_l,ub_l,options);
    %ci_l = nlparci(x_low,residuall,'jacobian',jacobianl);

    % Fit the high q region where the scaling changes from the rest of the fit;
    % fit this with the same Porod-type scaling as the low-q region
    [x_high,resnormh,residualh,~,~,~,jacobianh] = ...
        lsqnonlin(@(x) LSfunc(@(q_high,x) CorrLengthLowQ(q_high,x),q_high,I_sub_high,dI_sub_high,x),x0_h,lb_h,ub_h,options);
    %ci_h = nlparci(x_high,residualh,'jacobian',jacobianh);

    % Fit the high-q along with the low q, but with the values of A and n set
    % using the normal Correlation Length Model
    x0 = [x_low(1); x_low(2); 0.5; 2; 0.05; 22]; % [Inverse Angstrom q and cm I]
    lb = [x_low(1); x_low(2); 0; 0; 0.01; 0]; 
    ub = [x_low(1); x_low(2); inf; 5; 0.07; 100];
    %This takes a long time for some of my data
    [x_m,resnorm,residual,~,~,~,jacobian] = ...
        lsqnonlin(@(x) LSfunc(@(q,x) CorrLength(q,x),q_fit,I_sub_fit,dI_sub_fit,x),x0,lb,ub,options);
    %ci = nlparci(x_m,residual,'jacobian',jacobian);
    MSR = resnorml + resnormh + resnorm;
    
%     %% Output the parameters
%     fprintf('A = %.5f +/- %.5f\n',x_low(1),x_low(1)-ci_l(1,1));
%     fprintf('n = %.2f +/- %.1f\n',x_low(2),x_low(2)-ci_l(2,1));
%     fprintf('C = %.2f +/- %.2f\n',x(3),x(3)-ci(3,1));
%     fprintf('m = %.2f +/- %.1f\n',x(4),x(4)-ci(4,1));
%     fprintf('q0 = %.3f +/- %.3f A^-1\n',x(5),x(5)-ci(5,1));
%     fprintf('xi = %.0f +/- %.0f A\n',x(6),x(6)-ci(6,1));
%     fprintf('D = %.6f +/- %.5f\n',x_high(1),x_high(1)-ci_h(1,1));
%     fprintf('n = %.1f +/- %.1f\n',x_high(2),x_high(2)-ci_h(2,1));
% 
    %% Check fit
%     I_fit = CorrLength(logspace(-3,0),x);
%     I_fit_lowq = CorrLengthLowQ(logspace(-3,0),x_low);
%     I_fit_midq = CorrLengthMidQ(logspace(-3,0),x(3:end));
%     I_fit_highq = CorrLengthLowQ(logspace(-3,0),x_high);
%     bkg = bkg_avg*ones(1,50);

%     % Residuals
%     figure()
%     hold on;
%     box on;
%     plot(q_low*10,residuall,'.','MarkerSize',20,'Color','r');
%     plot(q_fit*10,residual,'.','MarkerSize',20,'Color','b');
%     xlabel("q [nm^-^1]",'FontWeight','bold');
%     ylabel('Residuals','FontWeight','bold');
%     xlim(fitq*10);
%     % ylim([-10 10]);
%     set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',1,'xscale','linear',...
%         'yscale','linear');
%     set(gcf,'Color','w','units','pixels','outerposition',[200 100 600 600]);
% 
%     figure()
%     box on;
%     hold on;
%     xlabel("q [nm^-^1]",'FontWeight','bold');
%     ylabel('I(q) [cm^-^1]','FontWeight','bold');
%     xlim(qcut*10);
%     ylim([0.01 10]);
%     set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',1,'xscale','log',...
%         'yscale','log');
%     set(gcf,'Color','w','units','pixels','outerposition',[200 100 600 600]);
%     errorbar(q*10,I_sub_solv,dI_sub_solv,'s','MarkerSize',10,'MarkerEdgeColor','k','LineWidth',1,'Color','k');
%     h = gobjects(3,1);
%     h(1) = plot(q*10,I_sub_solv,'s','MarkerSize',10,'MarkerEdgeColor','k','LineWidth',1,'Color','k');
%     saveas(gcf,'P4_expt.png')
%     h(2) = plot(logspace(-3,0)*10,bkg,'Color','b','LineWidth',2);
%     h(3) = plot(logspace(-3,0)*10,I_fit+bkg,'Color','r','LineWidth',2);
%     legend(h,'Expt','Bkg','Fit','Location','Northeast');
%     saveas(gcf,'P4_corr_fit.png')
%     h(4) = plot(logspace(-3,0)*10,I_fit_lowq+bkg,'--','Color','b','LineWidth',2);
%     legend(h, 'Expt','Bkg','Fit','LowQ Fit','Location','Northeast');
%     saveas(gcf,'P4_lowq.png')
%     h(5) = plot(logspace(-3,0)*10,I_fit_midq+bkg,'.','Color','b','MarkerSize',10);
%     legend(h,'Expt','Bkg','Fit','LowQ Fit','MidQ Fit','Location','Northeast');
%     saveas(gcf,'P4_midq.png')
%     h(6) = plot(logspace(-3,0)*10,I_fit_highq+bkg,'-.','Color','b','LineWidth',2);
%     legend(h,'Expt','Bkg','Fit','LowQ Fit','MidQ Fit','HighQFit','Location','Northeast');
%     saveas(gcf,'P4_highq.png')
%     saveas(gcf, 'P4_6p5_SANS_allfit.png')
%     export_fig(fullfile(hirespath,'P4_6p5_SANS_allfit'),'-r600', '-png');
    
elseif strcmp(model,'glgel')
    % Fitting to the Gauss-Lorentz gel:
    % I(q) = Igexp(-q^2Z^2/2)+Il/(1+q^2x^2)
    % Guess [Ig, Z [A], Il, x [A] ]
    x0 = [13; 228; 1.1; 23]; % [Inverse Angstrom q and cm I]
    lb = [-inf; 0; -inf; 0]; 
    ub = [inf; 1000; inf; 1000];
    options=optimoptions('lsqnonlin','MaxFunEvals',inf,'MaxIter',inf,...
        'TolFun',1e-12,'TolX',1e-12);
    [x,resnorm,residual,~,~,~,jacobian] = ...
        lsqnonlin(@(x) LSfunc(@(q,x) GaussLorentzGel(q,x),q_fit,I_sub_fit,dI_sub_fit,x),x0,lb,ub,options);
    MSR = resnorm;
    Z = x (3);
    xi = x (4);
    x_high = 0;
    %ci = nlparci(x,residual,'jacobian',jacobian);
    %% Output the parameters
%     fprintf('Ig = %.2f +/- %.2f cm^-1\n',x(1),x(1)-ci(1,1));
%     fprintf('Il = %.2f +/- %.1f cm^-1\n',x(3),x(3)-ci(3,1));
%     fprintf('Z = %.2f +/- %.2f Angstroms\n',x(2),x(2)-ci(2,1));
%     fprintf('x = %.2f +/- %.1f Angstroms\n',x(4),x(4)-ci(4,1));
    
%     %% Plot the fit
%     I_fit = GaussLorentzGel(logspace(-3,0),x);
%     bkg = bkg_avg*ones(1,50);
%     chi2 = sum(((GaussLorentzGel(q_fit,x)-I_sub_fit)./dI_sub_fit).^2);
%     % Residuals
%     figure()
%     hold on;
%     box on;
%     plot(q*10,residual,'.','MarkerSize',20,'Color','r');
%     xlabel(strcat("q [nm^-^1]"),'FontWeight','bold');
%     ylabel('Residuals','FontWeight','bold');
%     xlim(fitq*10);
%     % ylim([-10 10]);
%     set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',1,'xscale','linear',...
%         'yscale','linear');
%     set(gcf,'Color','w','units','pixels','outerposition',[200 100 600 600]);
% 
%     figure()
%     box on;
%     hold on;
%     xlabel(strcat("q [nm^-^1]"),'FontWeight','bold');
%     ylabel('I(q) [cm^-^1]','FontWeight','bold');
%     xlim(qcut*10);
%     ylim([0.01 10]);
%     set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',1,'xscale','log',...
%         'yscale','log');
%     set(gcf,'Color','w','units','pixels','outerposition',[200 100 600 600]);
%     errorbar(q*10,I_sub_solv,dI_sub_solv,'s','MarkerSize',10,'MarkerEdgeColor','k','LineWidth',1,'Color','k');
%     h = gobjects(3,1);
%     h(1) = plot(q*10,I_sub_solv,'s','MarkerSize',10,'MarkerEdgeColor','k','LineWidth',1,'Color','k');
%     h(2) = plot(logspace(-3,0)*10,bkg,'Color','b','LineWidth',2);
%     h(3) = plot(logspace(-3,0)*10,I_fit+bkg,'Color','r','LineWidth',2);
%     legend(h,'Expt','Bkg','Fit','Location','Northeast');
%     saveas(gcf,'P4_GLGel.png');
%     export_fig(fullfile(hirespath,'P4_GLGel'),'-r600', '-png');
elseif strcmp(model,'fsgel')
    % Fitting to the Gel (fine-scale polymer distribution) model:
    % I(q) = IL(1/(1+(((D+1/3)Q^2a1^2))^(D/2)) + IGexp(-Q^2a2^2)
    % Guess: [IG, IL, a1, a2, D]
    x0 = [11.5; 0.926; 13; 149; 3]; % [Inverse Angstrom q and cm I]
    lb = [0; 0; 0; 0; 0]; 
    ub = [20; 2; 500; 500; 5];
    options=optimoptions('lsqnonlin','MaxFunEvals',inf,'MaxIter',inf,...
        'TolFun',1e-12,'TolX',1e-12);
    [x,resnorm,residual,~,~,~,jacobian] = ...
        lsqnonlin(@(x) LSfunc(@(q,x) FSGel(q,x),q_fit,I_sub_fit,dI_sub_fit,x),x0,lb,ub,options);
    MSR = resnorm;
    %ci = nlparci(x,residual,'jacobian',jacobian);
%     %% Output the parameters
%     fprintf('IG = %.2f +/- %.2f cm^-1\n',x(1),x(1)-ci(1,1));
%     fprintf('IL = %.2f +/- %.2f cm^-1\n',x(2),x(2)-ci(2,1));
%     fprintf('a1 = %.2f +/- %.2f Angstroms\n',x(3),x(3)-ci(3,1));
%     fprintf('a2 = %.2f +/- %.1f Angstroms\n',x(4),x(4)-ci(4,1));
%     fprintf('Rg = %.2f +/- %.1f Angstroms\n',x(4)*sqrt(3),x(4)*sqrt(3)-ci(4,1)*sqrt(3));
%     fprintf('D = %.2f +/- %.1f\n',x(5),x(5)-ci(5,1));
    
%     %% Plot the fit
%     I_fit = FSGel(logspace(-3,0),x);
%     bkg = bkg_avg*ones(1,50);
%     chi2 = sum(((FSGel(q_fit,x)-I_sub_fit)./dI_sub_fit).^2);
%     % Residuals
%     figure()
%     hold on;
%     box on;
%     plot(q*10,residual,'.','MarkerSize',20,'Color','r');
%     xlabel(strcat("q [nm^-^1]"),'FontWeight','bold');
%     ylabel('Residuals','FontWeight','bold');
%     xlim(fitq*10);
%     % ylim([-10 10]);
%     set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',1,'xscale','linear',...
%         'yscale','linear');
%     set(gcf,'Color','w','units','pixels','outerposition',[200 100 600 600]);
% 
%     figure()
%     box on;
%     hold on;
%     xlabel(strcat("q [nm^-^1]"),'FontWeight','bold');
%     ylabel('I(q) [cm^-^1]','FontWeight','bold');
%     xlim(qcut*10);
%     ylim([0.01 10]);
%     set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',1,'xscale','log',...
%         'yscale','log');
%     set(gcf,'Color','w','units','pixels','outerposition',[200 100 600 600]);
%     h = gobjects(3,1);
%     errorbar(q*10,I_sub_solv,dI_sub_solv,'s','MarkerSize',10,'MarkerEdgeColor','k','LineWidth',1,'Color','k');
%     h(1) = plot(q*10,I_sub_solv,'s','MarkerSize',10,'MarkerEdgeColor','k','LineWidth',1,'Color','k');    
%     h(2) = plot(logspace(-3,0)*10,bkg,'Color','b','LineWidth',2);
%     h(3) = plot(logspace(-3,0)*10,I_fit+bkg,'Color','r','LineWidth',2);
%     legend('Expt','Bkg','Fit','Location','Northeast');
%     saveas(gcf,'P4_FSGel.png');
%     export_fig(fullfile(hirespath,'P4_FSGel'),'-r600', '-png');
end

%% LS function to minimize

function LS = LSfunc(func,q,I_expt,dI_expt,params)
%Obtain fit I data from correlation peak model (CorrLength)
I_fit = func(q,params); % any function in this code!!! :D :D :D 
%Calculate the residual
diff = I_fit - I_expt;
%Create the LS objective function
LS = diff./dI_expt;
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

%% MODEL: Gaussian Lorentzian Gel 
% This is a model that calculates scattering from a physical network,
% modeled as the sum of a low-q expoenntail decay plus a Lorentzian at
% higher q-values. Generally applied to gels.
% I(q) = Ig(0)exp(-q^2Z^2/2) + Il(0)/(1+q^2x^2)
% Guess in this order: [Ig, Z, Il, x] where the model claims that Z is the
% "static correlation" attributed to the frozen-in crosslinks of some gels,
% and x is the "dynamic correlation" attributed to the fluctuating polymer
% chain between crosslinks
function I_fit = GaussLorentzGel(q,params)
Ig = params(1); Il = params(3);
Z = params(2); x = params(4); % [in angstrom]
I_fit = Ig*exp((-q.^2*Z^2)./2) + Il./(1+q.^2*x^2);
end

%% MODEL: Gel (fine-scale polymer distribution)
% Unlike a concentrated polymer solution, the fine-scale polymer distribution 
% in a gel involves at least two characteristic length scales, a shorter 
% correlation length ( a1 ) to describe the rapid fluctuations in the 
% position of the polymer chains that ensure thermodynamic equilibrium, 
% and a longer distance (denoted here as a2 ) needed to account for the 
% static accumulations of polymer pinned down by junction points or clusters 
% of such points. The latter is derived from a simple Guinier function. 
% I(q) = IL(1/(1+(((D+1/3)Q^2a1^2))^(D/2)) + IGexp(-Q^2a2^2)
% Guess in this order: [IG, IL, a1, a2, D]
function I_fit = FSGel(q,params)
IG = params(1); IL = params(2); a1 = params(3); a2 = params(4); D = params(5);
I_fit = IL.*(1./(1+(((D+(1/3))*q.^2*a1^2)).^(D/2))) + IG*exp(-q.^2*a2^2);
end
end