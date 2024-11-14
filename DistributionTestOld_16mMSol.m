% this file reads the low Q and high Q raw data and tries to compare with gamma
% distribution
% low Q value q = 1.506607E-02
% High q value q = 2.500345E-01
clc
clear all
close all
%%
Number_Rep = 100; % for 0.01
Ilow = zeros(Number_Rep,2);
Ihigh = zeros(Number_Rep,2);
for i = 1:Number_Rep
sample_qlow = importdata(sprintf('16mMSolq_1_%d_1D.txt',i-1));
Ilow(i,1) = sample_qlow.data(sample_qlow.data(:,1) == 1.506607E-02,2);
Ilow(i,2) = sample_qlow.data(sample_qlow.data(:,1) == 1.506607E-02,3);
sample_qhigh = importdata(sprintf('16mMSolq_2_%d_1D.txt', i-1));
Ihigh(i,1) = sample_qhigh.data(sample_qhigh.data(:,1) == 2.500345E-01,2);
Ihigh(i,2) = sample_qhigh.data(sample_qhigh.data(:,1) == 2.500345E-01,3);
end


Data_pts = 100;
m = Ilow(1,1);
v = (Ilow(1,2))^2;
a = m^2/v;
b = v/m;
I_sim_gam = gamrnd(a,b,[Data_pts,1]);
I_sim_norm = normrnd(m,sqrt(v),[Data_pts 1]);
mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));
I_sim_lognorm = lognrnd(mu, sigma, [Data_pts 1]);
%param = fitdist(Ilow(:,1),'Gamma');
%I_sim_fit = gamrnd(param.a,param.b,[Data_pts,1]);

%% this section compares the different distribution for low Q

figure(1) % compares gamma distribution to SANS data
histogram(Ilow(:,1),'Normalization','probability')
hold on
histogram(I_sim_gam ,'Normalization','probability')
legend('Low Q SANS Data','Simulation from Gamma Dis')
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);

figure(2)
histogram(Ilow(:,1),'Normalization','probability')
hold on
histogram(I_sim_norm,'Normalization','probability')
legend('Low Q SANS Data','Simulation from Normal Dis')
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
%%
% figure(3)
% histfit(Ilow(:,1),10,'gamma')
% hold on
% legend('Low Q SANS Data','Fitted Gamma Dist')
% set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2);
% set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
%%
figure(4)
histogram(Ilow(:,1),'Normalization','probability')
hold on
histogram(I_sim_lognorm,'Normalization','probability')
legend('Low Q SANS Data','Simulation from Log Normal Distribution')
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);


%% this section simulates the high Q region
Data_pts = 100;
m_high = Ihigh(1,1);
v_high = (Ihigh(1,2))^2;
a_high = m_high^2/v_high;
b_high = v_high/m_high;
I_sim_gam_high = gamrnd(a_high,b_high,[Data_pts,1]);
I_sim_norm_high = normrnd(m_high,sqrt(v_high),[Data_pts 1]);
mu_high = log((m_high^2)/sqrt(v_high+m_high^2));
sigma_high = sqrt(log(v_high/(m_high^2)+1));
I_sim_lognorm_high = lognrnd(mu_high, sigma_high, [Data_pts 1]);
param_high = gamfit(Ihigh(:,1));

figure() % compares gamma distribution to SANS data
histogram(Ihigh(:,1),'Normalization','probability')
hold on
histogram(I_sim_gam_high ,'Normalization','probability')
legend('High Q SANS Data','Simulation from Gamma Dis')
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);

figure()
histogram(Ihigh(:,1),'Normalization','probability')
hold on
histogram(I_sim_norm_high,'Normalization','probability')
legend('High Q SANS Data','Simulation from Normal Dis')
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);

%%
figure()
histfit(Ihigh(:,1),10,'gamma')
legend('High Q SANS Data','Fitted Gamma Distribution')
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
%%
figure()
histogram(Ihigh(:,1),'Normalization','probability')
hold on
histogram(I_sim_lognorm_high,'Normalization','probability')
legend('High Q SANS Data','Simulation from Log Normal Distribution')
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2);
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);