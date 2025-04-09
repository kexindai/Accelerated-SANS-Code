close all
clc
clear
number = 1;
param_correlation = zeros(number,19);
param_GL = zeros(number,11);
param_FS = zeros(number,13);
param_DB = zeros(number,9);
param_OZLorentz = zeros(number,9);
%% Fitting all of the replicates in a for loop
for i = 1:number
    qcut = [0.005 0.6]; % [A^-1] min and max q to keep (from 0.005 to 0.6 look good)
    lowq = [0.005 0.011]; % [A^-1] very low q for the "clustering" region
    highq = [0.24 0.36]; % [A^-1] high q for the region where scaling changes
    vhighq = [0.36 0.6]; % [A^-1] very high q for the background (flattened)
    % Which qs to fit?
    % For Broad-peak Correlation Length model, use
    fitq = [0.005 0.6]; % [A^-1] q for fitting the "solvation" region (optional to use 0.6 as highest q)
    vp = 0.05; % volume fraction of P4 (corresponds to 6.5%--see labnoteb
    p4 = importdata(sprintf('P4_%d_merged.txt',i-1)); %csvread('6p5_P4_stitched.csv',2,0);
    % p4 = importdata('P4_merged.txt');
    % p4 = p4.data;
    buffer= importdata('Buffer_merged.txt'); %csvread('buffer_stitched.txt',2,0);
    buffer = buffer.data;
    % Truncate the data by getting rid of bad points at the edges of the plot
    % (using the qcut input)
    % Sample
    p4_cut(:,1) = p4(p4(:,1)>qcut(1) & p4(:,1)<qcut(2),1); % q [A^-1]
    q = p4_cut(:,1);
    p4_cut(:,2) = p4(p4(:,1)>qcut(1) & p4(:,1)<qcut(2),2); % I [cm^-1]
    p4_cut(:,3) = p4(p4(:,1)>qcut(1) & p4(:,1)<qcut(2),3); % dI [cm^-1]
    p4_cut(:,4) = p4(p4(:,1)>qcut(1) & p4(:,1)<qcut(2),3); % dq [A^-1]
    % Buffer
    buffer_cut(:,1) = buffer(p4(:,1)>qcut(1) & p4(:,1)<qcut(2),1); % q [A^-1]
    buffer_cut(:,2) = buffer(p4(:,1)>qcut(1) & p4(:,1)<qcut(2),2); % I [cm^-1]
    buffer_cut(:,3) = buffer(p4(:,1)>qcut(1) & p4(:,1)<qcut(2),3); % dI [cm^-1]
    buffer_cut(:,4) = buffer(p4(:,1)>qcut(1) & p4(:,1)<qcut(2),3); % dq [A^-1]
    
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
    figure('visible','off')
    box on;
    hold on;
    errorbar(q,I_sub,dI_sub,'s','MarkerSize',6,'MarkerEdgeColor','k','LineWidth',1,'Color','k');
    xlabel("q [nm^-^1]",'FontWeight','bold');
    ylabel('I(q) [cm^-^1]','FontWeight','bold');
    title('6.5% P4: 1D SANS, subtracted')
    xlim(qcut);
    ylim([0.001 10]);
    set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
    set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
    xlabel('q [nm^-^1]','FontWeight','bold');
    ylabel('I(q) [cm^-^1]','FontWeight','bold');
    
    
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
    csvwrite(sprintf('P4_stitched_sub_%d.csv',i-1),[q_fit,I_sub_fit, dI_sub_fit,p4_cut(:,4)]);
    %% Perform fit
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
    ci_l = nlparci(x_low,residuall,'jacobian',jacobianl);
    
    % Fit the high q region where the scaling changes from the rest of the fit;
    % fit this with the same Porod-type scaling as the low-q region
    [x_high,resnormh,residualh,~,~,~,jacobianh] = ...
        lsqnonlin(@(x) LSfunc(@(q_high,x) CorrLengthLowQ(q_high,x),q_high,I_sub_high,dI_sub_high,x),x0_h,lb_h,ub_h,options);
    ci_h = nlparci(x_high,residualh,'jacobian',jacobianh);
    
    % Fit the high-q along with the low q, but with the values of A and n set
    % using the normal Correlation Length Model
    % set the initial guess to 27 for xi
    x0 = [x_low(1); x_low(2); 0.5; 2; 0.02; 27]; % [Inverse Angstrom q and cm I]
    lb = [x_low(1); x_low(2); 0; 0; 0.01; 0];
    ub = [x_low(1); x_low(2); inf; 5; 0.07; 100];
    %This takes a long time for some of my data
    [x,resnorm,residual,~,~,~,jacobian] = ...
        lsqnonlin(@(x) LSfunc(@(q,x) CorrLength(q,x),q_fit,I_sub_fit,dI_sub_fit,x),x0,lb,ub,options);
    ci = nlparci(x,residual,'jacobian',jacobian);
    % LL = RSS/N
    MSR_correlation = sum((residual.*dI_sub_fit).^2)/length(residual);
    %Gaussian - maybe calculate the log likelihood instead?
    MSR_correlation_g = sum((residual.*dI_sub_fit).^2./(dI_sub_fit).^2);
    
    %% Output the parameters
    fprintf('A = %.5f +/- %.5f\n',x_low(1),x_low(1)-ci_l(1,1));
    fprintf('n = %.2f +/- %.1f\n',x_low(2),x_low(2)-ci_l(2,1));
    fprintf('C = %.2f +/- %.2f\n',x(3),x(3)-ci(3,1));
    fprintf('m = %.2f +/- %.1f\n',x(4),x(4)-ci(4,1));
    fprintf('q0 = %.3f +/- %.3f A^-1\n',x(5),x(5)-ci(5,1));
    fprintf('xi = %.0f +/- %.0f A\n',x(6),x(6)-ci(6,1));
    fprintf('D = %.6f +/- %.5f\n',x_high(1),x_high(1)-ci_h(1,1));
    fprintf('n = %.1f +/- %.1f\n',x_high(2),x_high(2)-ci_h(2,1));
    AIC_correlation = length(residual)*log(MSR_correlation) + 2*6;
    BIC_correlation = length(residual)*log(MSR_correlation) + 6 * log(length(residual));
    param_correlation(i,1) = x_low(1);
    param_correlation(i,2) = x_low(1)-ci_l(1,1);
    param_correlation(i,3) = x_low(2);
    param_correlation(i,4) = x_low(2)-ci_l(2,1);
    param_correlation(i,5) = x(3);
    param_correlation(i,6) = x(3)-ci(3,1);
    param_correlation(i,7) = x(4);
    param_correlation(i,8) = x(4)-ci(4,1);
    param_correlation(i,9) = x(5);
    param_correlation(i,10) = x(5)-ci(5,1);
    param_correlation(i,11) = x(6); % this is xi
    param_correlation(i,12) = x(6)-ci(6,1);
    param_correlation(i,13) = x_high(1);
    param_correlation(i,14) = x_high(1)-ci_h(1,1);
    param_correlation(i,15) = x_high(2);
    param_correlation(i,16) = x_high(2)-ci_h(2,1);
    param_correlation(i,17) = MSR_correlation;
    param_correlation(i,18) = AIC_correlation;
    param_correlation(i,19) = BIC_correlation;
    %% Check fit
    I_fit = CorrLength(logspace(-3,0),x);
    I_fit_lowq = CorrLengthLowQ(logspace(-3,0),x_low);
    I_fit_midq = CorrLengthMidQ(logspace(-3,0),x(3:end));
    I_fit_highq = CorrLengthLowQ(logspace(-3,0),x_high);
    bkg = bkg_avg*ones(1,50);
    
    % Residuals
    figure('visible','off')
    hold on;
    box on;
    plot(q_low*10,residuall,'.','MarkerSize',20,'Color','r');
    plot(q_fit*10,residual,'.','MarkerSize',20,'Color','b');
    xlabel("q [nm^-^1]",'FontWeight','bold');
    ylabel('Residuals','FontWeight','bold');
    xlim(fitq*10);
    ylim([-10 10]);
    set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
    set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
    
    figure('visible','on')
    box on;
    hold on;
    xlabel("q [nm^-^1]",'FontWeight','bold');
    ylabel('I(q) [cm^-^1]','FontWeight','bold');
    xlim(qcut*10);
    ylim([0.01 10]);
    set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
    %set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log');
    set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
    xlabel('q [nm^-^1]','FontWeight','bold');
    ylabel('I(q) [cm^-^1]','FontWeight','bold');
    errorbar(q*10,I_sub_solv,dI_sub_solv,'s','MarkerSize',10,'MarkerEdgeColor','k','LineWidth',1,'Color','k');
    %plot(q*10,I_sub_solv,'s','MarkerSize',10,'MarkerEdgeColor','k','LineWidth',1,'Color','m');
    h = gobjects(3,1);
    h(1) = plot(q*10,I_sub_solv,'s','MarkerSize',10,'MarkerEdgeColor','k','LineWidth',1,'Color','k');
    %     saveas(gcf,sprintf('P4_expt_%d.png',i-1))
    h(2) = plot(logspace(-3,0)*10,bkg,'Color','b','LineWidth',2);
    h(3) = plot(logspace(-3,0)*10,I_fit+bkg,'Color','r','LineWidth',2);
    legend(h,'Expt','Bkg','Fit','Location','Northeast');
    %     saveas(gcf,sprintf('P4_corr_fit_%d.png',i-1))
    h(4) = plot(logspace(-3,0)*10,I_fit_lowq+bkg,'--','Color','b','LineWidth',2);
    legend(h, 'Expt','Bkg','Fit','LowQ Fit','Location','Northeast');
    %     saveas(gcf,'P4_lowq.png')
    h(5) = plot(logspace(-3,0)*10,I_fit_midq+bkg,'.','Color','b','MarkerSize',10);
    legend(h,'Expt','Bkg','Fit','LowQ Fit','MidQ Fit','Location','Northeast');
    %    saveas(gcf,'P4_midq.png')
    h(6) = plot(logspace(-3,0)*10,I_fit_highq+bkg,'-.','Color','b','LineWidth',2);
    legend(h,'Expt','Bkg','Fit','LowQ Fit','MidQ Fit','HighQFit','Location','Northeast');
    %    saveas(gcf,'P4_highq.png')
    %saveas(gcf, sprintf('P4_6p5_SANS_allfit_%d.png',i-1))
    %    export_fig(fullfile(hirespath,'P4_6p5_SANS_allfit'),'-r600', '-png');
%% this section tests for the correlation between errors in P4
% Work flow
% 1. Take the model and the fitted parameters f(q, theta) where theta is
% the parmaters
% 2. Fix theta and add noise to the model. We will first test Gaussian
% noise (the assumption for MC bootstrapping)
% 3. We will create 100 replicates and fit the model with broad peak model
% 4. If nonlinear least squares is a biased estimator (errors are correlated), 
% then your theta’ will have expectation something other than theta

params = zeros(6,1);
params(1) = param_correlation(1); %A
params(2) = param_correlation(3); %n
params(3) = param_correlation(5); %C
params(4) = param_correlation(7); %m
params(5) = param_correlation(9); %q0
params(6) = param_correlation(11); %xi
% the noise is gaussian centered at 0 but porportional to the dI
Num_rep = 100;
noise = dI_sub_solv.*randn(length(q), Num_rep);
I_fit_data = CorrLength(q,params);
I_rep_data = repmat(I_fit_data, [1, Num_rep]) + noise;

% this plots the original data
figure()
loglog(q, I_sub_solv,'o')
hold on
loglog(q, I_rep_data + bkg_avg,'s')
hold on
loglog(q, I_fit_data + bkg_avg,'-','color','r','LineWidth',2)
hold on
plot(q, ones(length(q))* bkg_avg,'Color', 'b', 'LineWidth',2);
legend('original data','simulated data with noise', 'original data fit')
xlim(qcut);
ylim([0.001 100]);
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
xlabel('q [nm^-^1]','FontWeight','bold');
ylabel('I(q) [cm^-^1]','FontWeight','bold');

% now we fit the simulated data again

%% truncate the data
params_sim = zeros(6,100);
% I_fit_test = zeros(100,6);
for k = 1:Num_rep
    I_rep = I_rep_data(:,k);
    q_low = q(q > lowq(1) & q < lowq(2));
    I_sub_low = I_rep(q > lowq(1) & q < lowq(2));
    dI_sub_low = dI_sub(q > lowq(1) & q < lowq(2));
    
    q_fit = q(q > fitq(1) & q < fitq(2));
    I_sub_fit =  I_rep(q > fitq(1) & q < fitq(2));
    dI_sub_fit = dI_sub(q > fitq(1) & q < fitq(2));
    
    q_high = q(q > highq(1) & q < highq(2));
    I_sub_high =  I_rep(q > highq(1) & q < highq(2));
    dI_sub_high = dI_sub(q > highq(1) & q < highq(2));
    
    %% perform the fitting
    x0_l = [1e-5; 2]; % [Inverse Angstrom q and cm I]
    lb_l = [0; 0];
    ub_l = [0.1; 5];
    % [A, n]
    x0_h = [1e-5; 4]; % [Inverse Angstrom q and cm I]
    lb_h = [0; 0];
    ub_h = [0.1; 10];
    
    %% Here we fix A
    % [A, n]
    x0_l = [params(1); 2]; % [Inverse Angstrom q and cm I]
    lb_l = [params(1); 0];
    ub_l = [params(1); 5];
    % [A, n]
    x0_h = [params(1); 4]; % [Inverse Angstrom q and cm I]
    lb_h = [params(1); 0];
    ub_h = [params(1); 10];
    options=optimoptions('lsqnonlin','MaxFunEvals',inf,'MaxIter',inf,...
        'TolFun',1e-12,'TolX',1e-12);
    
    % Fit the low q region
    [x_low,resnorml,residuall,~,~,~,jacobianl] = ...
        lsqnonlin(@(x) LSfunc(@(q_low,x) CorrLengthLowQ(q_low,x),q_low,I_sub_low,dI_sub_low,x),x0_l,lb_l,ub_l,options);
    ci_l = nlparci(x_low,residuall,'jacobian',jacobianl);
    
    % Fit the high q region where the scaling changes from the rest of the fit;
    % fit this with the same Porod-type scaling as the low-q region
    [x_high,resnormh,residualh,~,~,~,jacobianh] = ...
        lsqnonlin(@(x) LSfunc(@(q_high,x) CorrLengthLowQ(q_high,x),q_high,I_sub_high,dI_sub_high,x),x0_h,lb_h,ub_h,options);
    ci_h = nlparci(x_high,residualh,'jacobian',jacobianh);
    
    % Fit the high-q along with the low q, but with the values of A and n set
    % using the normal Correlation Length Model
    % set the initial guess to 27 for xi
    x0 = [x_low(1); x_low(2); 0.5; 2; 0.02; 27]; % [Inverse Angstrom q and cm I]
    lb = [x_low(1); x_low(2); 0; 0; 0.01; 0];
    ub = [x_low(1); x_low(2); inf; 5; 0.07; 100];
    %This takes a long time for some of my data
    [x,resnorm,residual,~,~,~,jacobian] = ...
        lsqnonlin(@(x) LSfunc(@(q,x) CorrLength(q,x),q_fit,I_sub_fit,dI_sub_fit,x),x0,lb,ub,options);
    ci = nlparci(x,residual,'jacobian',jacobian);
    params_sim(:,k) = x;
    
    I_fit_test(k,:) = CorrLength(logspace(-3,0),x);
end
    % calculate percent difference
    bias_exp = mean(params_sim');
    param_percent_diff = (bias_exp' - params)./params * 100;
%% calculate the covariance matrix
cov_A_n = cov(params_sim(1,:), params_sim(2,:));
cov_A_C = cov(params_sim(1,:), params_sim(3,:));
cov_A_m0 = cov(params_sim(1,:), params_sim(4,:));
cov_A_q0 = cov(params_sim(1,:), params_sim(5,:));
cov_A_xi = cov(params_sim(1,:), params_sim(6,:));

cov_n_C = cov(params_sim(2,:), params_sim(3,:));
cov_n_m0 = cov(params_sim(2,:), params_sim(4,:));
cov_n_q0 = cov(params_sim(2,:), params_sim(5,:));
cov_n_xi = cov(params_sim(2,:), params_sim(6,:));

cov_C_m0 = cov(params_sim(3,:), params_sim(4,:));
cov_C_q0 = cov(params_sim(3,:), params_sim(5,:));
cov_C_xi = cov(params_sim(3,:), params_sim(6,:));

cov_m0_q0 = cov(params_sim(4,:), params_sim(5,:));
cov_m0_xi = cov(params_sim(4,:), params_sim(6,:));
cov_q0_xi = cov(params_sim(5,:), params_sim(6,:));
%% plot the covariance of each parameters in scatter plots

sz = 50;
% A and n
subplot(5,3,1)
scatter(params_sim(1,:), params_sim(2,:), sz, 'filled');

set(gca,'FontSize',20,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
xlabel('A values')
ylabel('n values')
% A and C
subplot(5,3,2)
scatter(params_sim(1,:), params_sim(3,:), sz, 'filled');
xlabel('A values')
ylabel('C values')
set(gca,'FontSize',20,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
% A and m0
subplot(5,3,3)
scatter(params_sim(1,:), params_sim(4,:), sz, 'filled');
xlabel('A values')
ylabel('m_{0} values')
set(gca,'FontSize',20,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
% A and q0
subplot(5,3,4)
scatter(params_sim(1,:), params_sim(5,:), sz, 'filled');
xlabel('A values')
ylabel('q_{0} values')
set(gca,'FontSize',20,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
% A and xi
subplot(5,3,5)
scatter(params_sim(1,:), params_sim(6,:), sz, 'filled');
xlabel('A values')
ylabel('\xi values')
set(gca,'FontSize',20,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
% n and C
subplot(5,3,6)
scatter(params_sim(2,:), params_sim(3,:), sz, 'filled');
xlabel('n values')
ylabel('C values')
set(gca,'FontSize',20,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
% n and m0
subplot(5,3,7)
scatter(params_sim(2,:), params_sim(4,:), sz, 'filled');
xlabel('n values')
ylabel('m_{0} values')
set(gca,'FontSize',20,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
% n and q0
subplot(5,3,8)
scatter(params_sim(2,:), params_sim(5,:), sz, 'filled');
xlabel('n values')
ylabel('q_{0} values')
set(gca,'FontSize',20,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
% n and xi
subplot(5,3,9)
scatter(params_sim(2,:), params_sim(6,:), sz, 'filled');
xlabel('n values')
ylabel('\xi values')
set(gca,'FontSize',20,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
% C and m0
subplot(5,3,10)

scatter(params_sim(3,:), params_sim(4,:), sz, 'filled');
xlabel('C values')
ylabel('m_{0} values')
set(gca,'FontSize',20,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
% C and q0
subplot(5,3,11)

scatter(params_sim(3,:), params_sim(5,:), sz, 'filled');
xlabel('C values')
ylabel('q_{0} values')
set(gca,'FontSize',20,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
% C and xi
subplot(5,3,12)

scatter(params_sim(3,:), params_sim(6,:), sz, 'filled');
xlabel('C values')
ylabel('\xi values')
set(gca,'FontSize',20,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
% m0 and q0
subplot(5,3,13)

scatter(params_sim(4,:), params_sim(5,:), sz, 'filled');
xlabel('m_{0} values')
ylabel('q_{0} values')
set(gca,'FontSize',20,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
% m0 and xi
subplot(5,3,14)

scatter(params_sim(4,:), params_sim(6,:), sz, 'filled');
xlabel('m_{0} values')
ylabel('\xi values')
set(gca,'FontSize',20,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
% q0 and xi
subplot(5,3,15)

scatter(params_sim(5,:), params_sim(6,:), sz, 'filled');
xlabel('q_{0} values')
ylabel('\xi values')
set(gca,'FontSize',20,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
%% plotting the data and the fit
figure()
loglog(q*10, I_rep_data + bkg_avg,'s')
hold on
loglog(logspace(-3,0)*10, I_fit_test + bkg_avg, 'LineWidth', 2)
hold on
plot(logspace(-3,0)*10,bkg_avg,'Color','b','LineWidth',2);
xlim(qcut*10);
ylim([0.001 100]);
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
xlabel('q [nm^-^1]','FontWeight','bold');
ylabel('I(q) [cm^-^1]','FontWeight','bold');

figure()
histogram(params_sim(1,:),'LineWidth',2);
xlabel('A values')
ylabel('Frequency')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);

figure()
histogram(params_sim(2,:),'LineWidth',2);
xlabel('n values')
ylabel('Frequency')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
figure()
histogram(params_sim(3,:),'LineWidth',2);
xlabel('C values')
ylabel('Frequency')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
figure()
histogram(params_sim(4,:),'LineWidth',2);
xlabel('m_{0} values')
ylabel('Frequency')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
figure()
histogram(params_sim(5,:),'LineWidth',2);
xlabel('q_{0} values')
ylabel('Frequency')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
figure()
histogram(params_sim(6,:),'LineWidth',2);
xlabel('\xi values')
ylabel('Frequency')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);

    %% Fitting to the Gauss-Lorentz gel:
    % I(q) = Igexp(-q^2Z^2/2)+Il/(1+q^2x^2)
    % Guess [Ig, Z [A], Il, x [A] ]
    x0 = [13; 228; 1.1; 22]; % [Inverse Angstrom q and cm I]
    lb = [-inf; 0; -inf; 0];
    ub = [inf; 1000; inf; 1000];
    options=optimoptions('lsqnonlin','MaxFunEvals',inf,'MaxIter',inf,...
        'TolFun',1e-12,'TolX',1e-12);
    [x_GL,resnorm,residual_GL,~,~,~,jacobian] = ...
        lsqnonlin(@(x) LSfunc(@(q,x) GaussLorentzGel(q,x),q_fit,I_sub_fit,dI_sub_fit,x),x0,lb,ub,options);
    MSR_GL = sum((residual_GL.*dI_sub_fit).^2)/length(residual);
    ci_GL = nlparci(x_GL,residual_GL,'jacobian',jacobian);
    %% Output the parameters
    %     fprintf('Ig = %.2f +/- %.2f cm^-1\n',x(1),x(1)-ci(1,1));
    %     fprintf('Il = %.2f +/- %.1f cm^-1\n',x(3),x(3)-ci(3,1));
    %     fprintf('Z = %.2f +/- %.2f Angstroms\n',x(2),x(2)-ci(2,1));
    %     fprintf('x = %.2f +/- %.1f Angstroms\n',x(4),x(4)-ci(4,1));
    AIC_GL = length(residual_GL)*log(MSR_GL) + 2*4;
    BIC_GL = length(residual_GL)*log(MSR_GL) + 4*log(length(residual_GL));
    param_GL(i,1) = x_GL(1);
    param_GL(i,2) = x_GL(1)-ci_GL(1,1);
    param_GL(i,3) = x_GL(2);
    param_GL(i,4) = x_GL(2)-ci_GL(2,1);
    param_GL(i,5) = x_GL(3);
    param_GL(i,6) = x_GL(3)-ci_GL(3,1);
    param_GL(i,7) = x_GL(4);
    param_GL(i,8) = x_GL(4)-ci_GL(4,1);
    param_GL(i,9) = MSR_GL;
    param_GL(i,10) = AIC_GL;
    param_GL(i,11) = BIC_GL;

    %% Plot the fit
    I_fit_GL = GaussLorentzGel(logspace(-3,0),x_GL);
    bkg = bkg_avg*ones(1,50);
    chi2 = sum(((GaussLorentzGel(q_fit,x)-I_sub_fit)./dI_sub_fit).^2);
    % Residuals
    figure('visible','off')
    hold on;
    box on;
    plot(q*10,residual_GL,'.','MarkerSize',20,'Color','r');
    xlabel(strcat("q [nm^-^1]"),'FontWeight','bold');
    ylabel('Residuals','FontWeight','bold');
    xlim(fitq*10);
    % ylim([-10 10]);
    set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
    set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
    %
    figure('visible','off')
    box on;
    hold on;
    xlabel(strcat("q [nm^-^1]"),'FontWeight','bold');
    ylabel('I(q) [cm^-^1]','FontWeight','bold');
    xlim(qcut*10);
    ylim([0.01 10]);
    set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
    set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
    errorbar(q*10,I_sub_solv,dI_sub_solv,'s','MarkerSize',10,'MarkerEdgeColor','k','LineWidth',1,'Color','k');
    h = gobjects(3,1);
    h(1) = plot(q*10,I_sub_solv,'s','MarkerSize',10,'MarkerEdgeColor','k','LineWidth',1,'Color','k');
    h(2) = plot(logspace(-3,0)*10,bkg,'Color','b','LineWidth',2);
    h(3) = plot(logspace(-3,0)*10,I_fit_GL+bkg,'Color','r','LineWidth',2);
    legend(h,'Expt','Bkg','Fit','Location','Northeast');
    %saveas(gcf,sprintf('P4_GLGel_%d.png',i-1));
 %% 
    % Fitting to the Gel (fine-scale polymer distribution) model:
    % I(q) = IL(1/(1+(((D+1/3)Q^2a1^2))^(D/2)) + IGexp(-Q^2a2^2)
    % Guess: [IG, IL, a1, a2, D]
    x0 = [11.5; 0.926; 13; 149; 3]; % [Inverse Angstrom q and cm I]
    lb = [0; 0; 0; 0; 0];
    ub = [20; 2; 500; 500; 5];
    options=optimoptions('lsqnonlin','MaxFunEvals',inf,'MaxIter',inf,...
        'TolFun',1e-12,'TolX',1e-12);
    [x_FS,resnorm,residual_FS,~,~,~,jacobian] = ...
        lsqnonlin(@(x) LSfunc(@(q,x) FSGel(q,x),q_fit,I_sub_fit,dI_sub_fit,x),x0,lb,ub,options);
    ci_FS = nlparci(x_FS,residual_FS,'jacobian',jacobian);
    MSR_FS = sum((residual_FS.*dI_sub_fit).^2)/length(residual_FS);
    AIC_FS = length(residual_FS)*log(MSR_FS) + 2*5;
    BIC_FS = length(residual_FS)*log(MSR_FS) + 5*log(length(residual_FS));
    param_FS(i,1) = x_FS(1);
    param_FS(i,2) = x_FS(1)-ci_FS(1,1);
    param_FS(i,3) = x_FS(2);
    param_FS(i,4) = x_FS(2)-ci_FS(2,1);
    param_FS(i,5) = x_FS(3);
    param_FS(i,6) = x_FS(3)-ci_FS(3,1);
    param_FS(i,7) = x_FS(4);
    param_FS(i,8) = x_FS(4)-ci_FS(4,1);
    param_FS(i,9) = x_FS(4);
    param_FS(i,10) = x_FS(4)-ci_FS(4,1);
    param_FS(i,11) = MSR_FS;
    param_FS(i,12) = AIC_FS;
    param_FS(i,13) = BIC_FS;
    figure('visible','off')
    I_fit_FS = FSGel(logspace(-3,0),x_FS);
    bkg = bkg_avg*ones(1,50);
    chi2 = sum(((FSGel(q_fit,x)-I_sub_fit)./dI_sub_fit).^2);
    box on;
    hold on;
    xlabel(strcat("q [nm^-^1]"),'FontWeight','bold');
    ylabel('I(q) [cm^-^1]','FontWeight','bold');
    xlim(qcut*10);
    ylim([0.01 10]);
    set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
    set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
    errorbar(q*10,I_sub_solv,dI_sub_solv,'s','MarkerSize',10,'MarkerEdgeColor','k','LineWidth',1,'Color','k');
    h = gobjects(3,1);
    h(1) = plot(q*10,I_sub_solv,'s','MarkerSize',10,'MarkerEdgeColor','k','LineWidth',2,'Color','k');
    h(2) = plot(logspace(-3,0)*10,bkg,'Color','b','LineWidth',2);
    h(3) = plot(logspace(-3,0)*10,I_fit_FS+bkg,'Color','r','LineWidth',2);
    legend(h,'Expt','Bkg','Fit','Location','Northeast');
    saveas(gcf,sprintf('P4_FSGel_%d.png',i-1));
    
       %% Fitting to Debye Bueche
% I(q) = scale*L^3/(1+(q*L)^2)^2 + bgk scale = 8*pi*phi*(1-phi)*delrho^2
% scale = params(1); L = params(2); bkg = params(3);
    x0 = [1; 20; 0.001];
    lb = [0; 0; 0];
    ub = [inf; inf; inf];
    options=optimoptions('lsqnonlin','MaxFunEvals',inf,'MaxIter',inf,...
        'TolFun',1e-12,'TolX',1e-12);
    [x_DB,resnorm,residual_DB,~,~,~,jacobian] = ...
        lsqnonlin(@(x) LSfunc(@(q,x) DB(q,x),q_fit,I_sub_fit,dI_sub_fit,x),x0,lb,ub,options);
    ci_DB = nlparci(x_DB,residual_DB,'jacobian',jacobian);
    MSR_DB = sum((residual_DB.*dI_sub_fit).^2)/length(residual_DB);
    AIC_DB = length(residual_DB)*log(MSR_DB) + 2*3;
    BIC_DB = length(residual_DB)*log(MSR_DB) + 3*log(length(residual_DB));
    param_DB(i,1) = x_DB(1);
    param_DB(i,2) = x_DB(1)-ci_DB(1,1);
    param_DB(i,3) = x_DB(2);
    param_DB(i,4) = x_DB(2)-ci_DB(2,1);
    param_DB(i,5) = x_DB(3);
    param_DB(i,6) = x_DB(3)-ci_DB(3,1);
    param_DB(i,7) = MSR_DB;
    param_DB(i,8) = AIC_DB;
    param_DB(i,9) = BIC_DB;
    figure('visible','off')
    I_fit_DB = DB(logspace(-3,0),x_DB);
    bkg = bkg_avg*ones(1,50);
    chi2 = sum(((OZLorentz(q_fit,x)-I_sub_fit)./dI_sub_fit).^2);
    box on;
    hold on;
    xlabel(strcat("q [nm^-^1]"),'FontWeight','bold');
    ylabel('I(q) [cm^-^1]','FontWeight','bold');
    xlim(qcut*10);
    ylim([0.01 10]);
    set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
    set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
    errorbar(q*10,I_sub_solv,dI_sub_solv,'s','MarkerSize',10,'MarkerEdgeColor','k','LineWidth',1,'Color','k');
    h = gobjects(3,1);
    h = gobjects(3,1);
    h(1) = plot(q*10,I_sub_solv,'s','MarkerSize',10,'MarkerEdgeColor','k','LineWidth',2,'Color','k');
    h(2) = plot(logspace(-3,0)*10,bkg,'Color','b','LineWidth',3);
    h(3) = plot(logspace(-3,0)*10,I_fit_DB+bkg,'Color','r','LineWidth',3);
    legend(h,'Expt','Bkg','Fit','Location','Northeast');
    %saveas(gcf,sprintf('P4_DB_%d.png',i-1));
    
%% Fitting to the Ornstein Zernike and squared Lorentz
    % I(q) = scale/(1+(qL)^2) + bkg
    % L is the screening length
    x0 = [1; 20; 0.001];
    lb = [0; 0; 0];
    ub = [inf; inf; inf];
    options=optimoptions('lsqnonlin','MaxFunEvals',inf,'MaxIter',inf,...
        'TolFun',1e-12,'TolX',1e-12);
    [x_OZLorentz,resnorm,residual_OZLorentz,~,~,~,jacobian] = ...
        lsqnonlin(@(x) LSfunc(@(q,x) OZLorentz(q,x),q_fit,I_sub_fit,dI_sub_fit,x),x0,lb,ub,options);
    ci_OZLorentz = nlparci(x_OZLorentz,residual_OZLorentz,'jacobian',jacobian);
    MSR_OZLorentz = sum((residual_OZLorentz.*dI_sub_fit).^2)/length(residual_OZLorentz);
    AIC_OZLorentz = length(residual_OZLorentz)*log(MSR_OZLorentz) + 2*3;
    BIC_OZLorentz = length(residual_OZLorentz)*log(MSR_OZLorentz) + 3*log(length(residual_OZLorentz));
    param_OZLorentz(i,1) = x_OZLorentz(1);
    param_OZLorentz(i,2) = x_OZLorentz(1)-ci_OZLorentz(1,1);
    param_OZLorentz(i,3) = x_OZLorentz(2);
    param_OZLorentz(i,4) = x_OZLorentz(2)-ci_OZLorentz(2,1);
    param_OZLorentz(i,5) = x_OZLorentz(3);
    param_OZLorentz(i,6) = x_OZLorentz(3)-ci_OZLorentz(3,1);
    param_OZLorentz(i,7) = MSR_OZLorentz;
    param_OZLorentz(i,8) = AIC_OZLorentz;
    param_OZLorentz(i,9) = BIC_OZLorentz;
    figure('visible','off')
    I_fit_OZLorentz = OZLorentz(logspace(-3,0),x_OZLorentz);
    bkg = bkg_avg*ones(1,50);
    chi2 = sum(((OZLorentz(q_fit,x)-I_sub_fit)./dI_sub_fit).^2);
    box on;
    hold on;
    xlabel(strcat("q [nm^-^1]"),'FontWeight','bold');
    ylabel('I(q) [cm^-^1]','FontWeight','bold');
    xlim(qcut*10);
    ylim([0.01 10]);
    set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
    set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
    errorbar(q*10,I_sub_solv,dI_sub_solv,'s','MarkerSize',10,'MarkerEdgeColor','k','LineWidth',1,'Color','k');
    h = gobjects(3,1);
    h(1) = plot(q*10,I_sub_solv,'s','MarkerSize',10,'MarkerEdgeColor','k','LineWidth',2,'Color','k');
    h(2) = plot(logspace(-3,0)*10,bkg,'Color','b','LineWidth',3);
    h(3) = plot(logspace(-3,0)*10,I_fit_OZLorentz+bkg,'Color','r','LineWidth',3);
    legend(h,'Expt','Bkg','Fit','Location','Northeast');
    %saveas(gcf,sprintf('P4_OZLorentz_%d.png',i-1));
end
%% Saving the parameters
% correlation model
T_correlation = table(param_correlation);
T_new_correlation = splitvars(T_correlation, 'param_correlation', 'NewVariableNames',{'A','A_error','n','n_error',...
    'C','C_error','m','m_error','q0','q0_error', 'xi', 'xi_error','D','D_error','n2','n2 error','MSR','AIC','BIC'});
%writetable(T_new_correlation, 'P4_correlation_parameters_0p01.csv');
%csvwrite('P4_GL_parameters.csv',param_GL);
%csvwrite('P4_FS_parameters.csv',param_FS);
%csvwrite('P4_DB_parameters.csv',param_DB);
%csvwrite('P4_OZLorentz_parameters.csv',param_OZLorentz);

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

%% MODEL: Ornstein Zernike and squared Lorentz
% I(q) = scale/(1+(qL)^2) + bkg
% L is the screening length
function I_fit = OZLorentz(q,params)
scale = params(1); L = params(2); bkg = params(3);
I_fit = scale./( 1 + (q*L).^2) + bkg;
end

%% MODEL: Debye Bueche
% I(q) = scale*L^3/(1+(q*L)^2)^2 + bgk scale = 8*pi*phi*(1-phi)*delrho^2
function I_fit = DB(q,params)
scale = params(1); L = params(2); bkg = params(3);
I_fit = scale*L^3./(1+(q*L).^2).^2 + bkg;
end

