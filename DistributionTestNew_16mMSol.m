% this file reads the low Q and high Q raw data and tries to compare with gamma
% distribution
% low Q value q = 1.506607E-02
% High q value q = 2.500345E-01
clc
clear all
close all
%%
Number_Rep = 100; % for 0.01
Ilow = zeros(Number_Rep,3);
Ihigh = zeros(Number_Rep,3);
for i = 1:Number_Rep
sample_qlow = importdata(sprintf('16mMSolq_1_%d_1D.txt',i-1));
% intensity
Ilow(i,1) = sample_qlow.data(sample_qlow.data(:,1) == 1.506607E-02,2);
% dI
Ilow(i,2) = sample_qlow.data(sample_qlow.data(:,1) == 1.506607E-02,3);
% dQ
Ilow(i,3) = sample_qlow.data(sample_qlow.data(:,1) == 1.506607E-02,4);
sample_qhigh = importdata(sprintf('16mMSolq_2_%d_1D.txt', i-1));
Ihigh(i,1) = sample_qhigh.data(sample_qhigh.data(:,1) == 2.500345E-01,2);
Ihigh(i,2) = sample_qhigh.data(sample_qhigh.data(:,1) == 2.500345E-01,3);
Ihigh(i,3) = sample_qhigh.data(sample_qhigh.data(:,1) == 2.500345E-01,4);
end


%% this part calculates from the data the mean and std. This is because at low Q the error was larger than we expect due to additional error added from our data reduction. 
rng("default")
m_calc = mean(Ilow(:,1));
v_raw = mean(Ilow(:,2).^2); % the variance from data reduction code
v_calc = (std(Ilow(:,1)))^2;
percentdiff_calc = (v_raw - v_calc)/v_calc; % we should report this
%%
a = m_calc^2/v_calc;
b = v_calc/m_calc;
I_sim_gam = gamrnd(a,b,[Number_Rep,1]);
I_sim_norm = normrnd(m_calc,sqrt(v_calc),[Number_Rep 1]);
mu = log((m_calc^2)/sqrt(v_calc+m_calc^2));
sigma = sqrt(log(v_calc/(m_calc^2)+1));
I_sim_lognorm = lognrnd(mu, sigma, [Number_Rep 1]);
%param = fitdist(Ilow(:,1),'Gamma');
%I_sim_fit = gamrnd(param.a,param.b,[Number_Rep,1]);
T_save_low = table(percentdiff_calc, a, b);
writetable(T_save_low, 'Low Q output values.csv')
%% this section compares the different distribution for low Q

figure(1) % compares gamma distribution to SANS data
h_data = histogram(Ilow(:,1),'Normalization','probability');
hold on
h_gamma = histogram(I_sim_gam , h_data.BinEdges, 'Normalization','probability');
%h_gamma.NumBins = 50;
legend('Low Q SANS Data','Simulation from Gamma Dis')
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
saveas(gcf, 'lowQ_gamma_calc_new.png')
figure(2)
histogram(Ilow(:,1),'Normalization','probability')
hold on
h_normal = histogram(I_sim_norm,h_data.BinEdges,'Normalization','probability');
legend('Low Q SANS Data','Simulation from Normal Dis')
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
saveas(gcf, 'lowQ_normal_calc_new.png')
% figure(3)
% histfit(Ilow(:,1),10,'gamma')
% hold on
% legend('Low Q SANS Data','Fitted Gamma Dist')
% set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2);
% set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
% saveas(gcf, 'lowQ_FittedGamma_new.png')
figure(4)
histogram(Ilow(:,1),'Normalization','probability')
hold on
histogram(I_sim_lognorm, h_data.BinEdges,'Normalization','probability')
legend('Low Q SANS Data','Simulation from Log Normal Distribution')
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
saveas(gcf, 'lowQ_LogNorm_calc_new.png')

%% Plot qq plot of dQ and dI
figure
qqplot(Ilow(:,3), Ilow(:,2))
figure
qqplot(Ilow(:,2))
%% this section simulates the high Q region
Number_Rep = 100;
m_high = mean(Ihigh(:,1));
v_high = (std(Ihigh(:,1)))^2;
v_high_raw = mean(Ihigh(:,2).^2);
percentdiff_high = (v_high - v_high_raw)/v_high; % we should report this
a_high = m_high^2/v_high;
b_high = v_high/m_high;
I_sim_gam_high = gamrnd(a_high,b_high,[Number_Rep,1]);
I_sim_norm_high = normrnd(m_high,sqrt(v_high),[Number_Rep 1]);
mu_high = log((m_high^2)/sqrt(v_high+m_high^2));
sigma_high = sqrt(log(v_high/(m_high^2)+1));
I_sim_lognorm_high = lognrnd(mu_high, sigma_high, [Number_Rep 1]);
param_high = gamfit(Ihigh(:,1));
a_fit_high = param_high(1);
b_fit_high = param_high(2);
T_save_high = table(percentdiff_high,a_fit_high, a_high, b_fit_high, b_high);
writetable(T_save_high, 'High Q output values.csv')

figure() % compares gamma distribution to SANS data
h_data_high = histogram(Ihigh(:,1),'Normalization','probability')
hold on
histogram(I_sim_gam_high ,h_data_high.BinEdges, 'Normalization','probability')
legend('High Q SANS Data','Simulation from Gamma Dis')
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
saveas(gcf, 'highQ_Gamma_calc_new.png')
figure()
histogram(Ihigh(:,1),'Normalization','probability')
hold on
histogram(I_sim_norm_high,h_data_high.BinEdges, 'Normalization','probability')
legend('High Q SANS Data','Simulation from Normal Dis')
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
saveas(gcf, 'highQ_Normal_calc_new.png')
%%
figure()
histfit(Ihigh(:,1),10,'gamma')
legend('High Q SANS Data','Fitted Gamma Distribution')
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
saveas(gcf, 'highQ_Gamma_Fitted_new.png')
%%
figure()
histogram(Ihigh(:,1),'Normalization','probability')
hold on
histogram(I_sim_lognorm_high,h_data_high.BinEdges, 'Normalization','probability')
legend('High Q SANS Data','Simulation from Log Normal Distribution')
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
saveas(gcf, 'highQ_LogNormal_calc_new.png')