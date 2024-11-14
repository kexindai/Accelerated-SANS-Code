clc
clear
param_F127_mean = csvread('Bootstrap_F127_predict_mean.csv'); 
param_F127_std = csvread('Bootstrap_F127_predict_std.csv'); 
param_exp = csvread('F-127paramest.csv',1,0);
LCL_R = param_exp(end,3) - param_exp(end,3) * 0.05;
UCL_R = param_exp(end,3) + param_exp(end,3) * 0.05;
LCL_Rg = param_exp(end,4) - param_exp(end,4) * 0.05;
UCL_Rg = param_exp(end,4) + param_exp(end,4) * 0.05;

param_mean = reshape(param_F127_mean, 8, [])';
param_std = reshape(param_F127_std, 8, [])';
frac = [0.01 0.025 0.05 0.1 0.25 0.5 1];
% this figure plots prediction from data 0.01 with experimental data this
% plots R or Rc
figure(1)
errorbar(frac, param_exp(:,3), param_exp(:,9), 's','MarkerSize',10,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2)
hold on
errorbar(frac, param_mean(:,2),param_std(:,2),'o','MarkerSize',10,'Color','[0 0.4470 0.7410]','MarkerEdgeColor','[0 0.4470 0.7410]','LineWidth',2)
hold on
plot(frac,UCL_R*ones(length(frac)),'--','LineWidth',2)
hold on
plot(frac,LCL_R*ones(length(frac)),'--','LineWidth',2)
hold on

legend('Experiment','Bootstrap')
xlim([0.008 1])
ylim([15 25])
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2,'xscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
xlabel('fraction of counts','FontWeight','bold');
ylabel('R [Å]','FontWeight','bold');
%% this plots Rg
figure(2)
errorbar(frac, param_exp(:,4), param_exp(:,10), 's','MarkerSize',10,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2)
hold on
errorbar(frac, param_mean(:,3),param_std(:,3),'o','MarkerSize',10,'Color','[0 0.4470 0.7410]','MarkerEdgeColor','[0 0.4470 0.7410]','LineWidth',2)
hold on
plot(frac,UCL_Rg*ones(length(frac)),'--','LineWidth',2)
hold on
plot(frac,LCL_Rg*ones(length(frac)),'--','LineWidth',2)
hold on
legend('Experiment','Bootstrap')
xlim([0.008 1])
ylim([0 20])
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2,'xscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
xlabel('fraction of counts','FontWeight','bold');
ylabel('Rg [Å]','FontWeight','bold');

%% this compares the predicted scattering data from 0.01 frac to experimental full counts

data_sim = importdata('F127_0_merged_asim_100.txt');
data_exp = importdata('F127_exp.txt');
figure()
h = gobjects(2, 1);
h(1) = errorbar(data_sim(:,1),data_sim(:,2), data_sim(:,3),'o','MarkerSize',10,'Color','[0 0.4470 0.7410]','MarkerEdgeColor','[0 0.4470 0.7410]','LineWidth',2);
hold on
h(2) = errorbar(data_exp(:,1),data_exp(:,2), data_exp(:,3),'s','MarkerSize',10,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2);

xlabel('q [nm^-^1]','FontWeight','bold');
ylabel('I(q) [cm^-^1]','FontWeight','bold');
legend(h([2 1]),{'Experimental Data','Simulated Data'})
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);

figure()
loglog(data_exp(:,1)*10, data_exp(:,3),'s','MarkerSize',10,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2)
hold on
loglog(data_sim(:,1)*10, data_sim(:,3),'o','MarkerSize',10,'Color','[0 0.4470 0.7410]','MarkerEdgeColor','[0 0.4470 0.7410]','LineWidth',2)

xlabel('q [nm^-^1]','FontWeight','bold');
ylabel('I(q) [cm^-^1]','FontWeight','bold');
legend('Experimental Error','Simulated Error')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);

%% this compares the boostrap error, CI of experimental data, and experimental mean and STD from true replicates
frac = [0.01 0.025 0.05 0.1 0.25 0.5 1];
Exp_param = csvread('F-127paramest.csv',1,0);
LCL_R = Exp_param(end,3) - Exp_param(end,3) * 0.05;
UCL_R = Exp_param(end,3) + Exp_param(end,3) * 0.05;
LCL_Rg = Exp_param(end,4) - Exp_param(end,4) * 0.05;
UCL_Rg = Exp_param(end,4) + Exp_param(end,4) * 0.05;
BSparam_0p01 = csvread('Bootstrap_F127_0p01.csv');
BSparam_0p025 = csvread('Bootstrap_F127_0p025.csv');
BSparam_0p05 = csvread('Bootstrap_F127_0p05.csv');
BSparam_0p1 = csvread('Bootstrap_F127_0p1.csv');
BSparam_0p25 = csvread('Bootstrap_F127_0p25.csv');
BSparam_0p5 = csvread('Bootstrap_F127_0p5.csv');
BSparam_1 = csvread('Bootstrap_F127_1.csv');
BSparam_mean = [mean(BSparam_0p01); mean(BSparam_0p025); mean(BSparam_0p05);
           mean(BSparam_0p1); mean(BSparam_0p25); mean(BSparam_0p5); mean(BSparam_1)];
BSparam_std = [std(BSparam_0p01); std(BSparam_0p025); std(BSparam_0p05);
           std(BSparam_0p1); std(BSparam_0p25); std(BSparam_0p5); std(BSparam_1)];
figure() % plot R
errorbar(frac, Exp_param(:,3), Exp_param(:,9), 's','MarkerSize',10,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2)
hold on
errorbar(frac, BSparam_mean(:,2), BSparam_std(:,2), 'o','MarkerSize',10,'Color','[0 0.4470 0.7410]','MarkerEdgeColor','[0 0.4470 0.7410]','LineWidth',2)
hold on
plot(frac,UCL_R*ones(length(frac)),'--','LineWidth',2)
hold on
plot(frac,LCL_R*ones(length(frac)),'--','LineWidth',2)

legend('Experiment','Bootstrap')
xlim([0.008 1])
ylim([15 25])
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2,'xscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
xlabel('fraction of counts','FontWeight','bold');
ylabel('R [Å]','FontWeight','bold');
figure() % plot Rg
errorbar(frac, Exp_param(:,4), Exp_param(:,10), 's','MarkerSize',10,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2)
hold on
errorbar(frac, BSparam_mean(:,3), BSparam_std(:,3), 'o','MarkerSize',10,'Color','[0 0.4470 0.7410]','MarkerEdgeColor','[0 0.4470 0.7410]','LineWidth',2)
hold on
plot(frac,UCL_Rg*ones(length(frac)),'--','LineWidth',2)
hold on
plot(frac,LCL_Rg*ones(length(frac)),'--','LineWidth',2)
legend('Experiment','Bootstrap')
xlim([0.008 1])
ylim([0 20])
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2,'xscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
xlabel('fraction of counts','FontWeight','bold');
ylabel('Rg [Å]','FontWeight','bold');

