clc
clear 
close all
%% this section reads all the boostrap files under Bootstrap folder
frac = [0.01 0.025 0.05 0.1 0.25 0.5 1];
exp = readtable('F127 exp param_rev.xlsx');
%%
data_{1} = csvread('Bootstrap_F127_0p01.csv');
data_{2} = csvread('Bootstrap_F127_0p025.csv');
data_{3} = csvread('Bootstrap_F127_0p05.csv');
data_{4} = csvread('Bootstrap_F127_0p1.csv');
data_{5} = csvread('Bootstrap_F127_0p25.csv');
data_{6} = csvread('Bootstrap_F127_0p5.csv');
data_{7} = csvread('Bootstrap_F127_1.csv');

% we want to plot Rc and Rg with mean and STDR
Rc = zeros(length(frac),1);
Rc_error_bootstrap = zeros(length(frac),1);
Rc_error_ci = zeros(length(frac),1);
Rg = zeros(length(frac),1);
Rg_error_bootstrap = zeros(length(frac),1);
Rg_error_ci = zeros(length(frac),1);

for i = 1:length(frac)
    Rc (i,1) = mean(data_{i}(:,2));
    Rc_exp(i,1) = data_{i}(1,2);
    Rc_error_bootstrap (i,1) = std(data_{i}(:,2));
    Rc_error_ci (i,1) = data_{i}(1,10);
    Rg (i,1) = mean(data_{i}(:,3));
    Rg_exp(i,1) = data_{i}(1,3);
    Rg_error_bootstrap (i,1) = std(data_{i}(:,3));
    Rg_error_ci (i,1) = data_{i}(1,11);
end
%% 10 and 11 are the CI for Rc and Rg
UCL_Rg = data_{7}(1,3) + data_{7}(1,3) * 0.1;
LCL_Rg = data_{7}(end,3) - data_{7}(end,3) * 0.1;
UCL_Rc = data_{7}(end,2) + data_{7}(end,2) * 0.05;
LCL_Rc = data_{7}(end,2) - data_{7}(end,2) * 0.05;
 % plots the RC FIRST
% this parts plots the exp only 
figure()
errorbar(frac,Rc_exp, Rc_error_ci,'s','MarkerSize',15,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',4)
hold on
errorbar(frac,Rc, Rc_error_bootstrap,'^','MarkerSize',15,'Color','[0 0.4470 0.7410]','MarkerEdgeColor','[0 0.4470 0.7410]','LineWidth',4)
hold on
errorbar(frac, exp{:,2}, exp{:,3},'o', 'MarkerSize',15,'Color','[0.8500 0.3250 0.0980]','MarkerEdgeColor','[0.8500 0.3250 0.0980]','LineWidth',4)
hold on
plot(frac,UCL_Rc*ones(length(frac)),'--','LineWidth',3, 'Color', [0.4940 0.1840 0.5560])
hold on
plot(frac,LCL_Rc*ones(length(frac)),'--','LineWidth',3, 'Color', [0.4940 0.1840 0.5560])
hold on
legend('Experiment with 95% CI','Bootstrap','Experimental Average and STD')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
xlabel('Fraction Of Counts','FontWeight','bold');
ylabel('R [Å]','FontWeight','bold');
xlim([0.009 1.0])
ylim([10 70])

%% plots the RG
figure()
errorbar(frac, Rg_exp, Rg_error_ci,'s','MarkerSize',15,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',4)
hold on
errorbar(frac, Rg, Rg_error_bootstrap,'^','MarkerSize',15,'Color','[0 0.4470 0.7410]','MarkerEdgeColor','[0 0.4470 0.7410]','LineWidth',4)
hold on
errorbar(frac, exp{:,4}, exp{:,5},'o', 'MarkerSize',15,'Color','[0.8500 0.3250 0.0980]','MarkerEdgeColor','[0.8500 0.3250 0.0980]','LineWidth',4)
hold on
plot(frac,UCL_Rg*ones(length(frac)),'--','LineWidth',3, 'Color', [0.4940 0.1840 0.5560])
hold on
plot(frac,LCL_Rg*ones(length(frac)),'--','LineWidth',3, 'Color', [0.4940 0.1840 0.5560])
hold on
legend('Experiment with 95% CI','Bootstrap','Experimental Average and STD')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
xlabel('Fraction Of Counts','FontWeight','bold');
ylabel('Rg [Å]','FontWeight','bold');
xlim([0.009 1.0])
ylim([-0.05 55])

%% this additional plots the scaling of error from MC bootstrapping and compares with experimental error from fitting
% Rc
figure()
plot(frac, Rc_error_ci,'s','MarkerSize',15,'MarkerEdgeColor','k','LineWidth',4)
hold on
plot(frac, exp{:,3},'o','MarkerSize',15,'MarkerEdgeColor','[0.8500 0.3250 0.0980]','LineWidth',4)
hold on
plot(frac, Rc_error_bootstrap,'^','MarkerSize',15,'MarkerEdgeColor','[0 0.4470 0.7410]','LineWidth',4)
hold on
scale_Rc = lsqnonlin(@(scale)fitscalingfactor(frac, Rc_error_bootstrap , scale), 1); 
scale_Rc = 1;
frac_fit = linspace(0.01,1,1001);
plot(frac_fit,1./sqrt(frac_fit)*scale_Rc,'-','LineWidth',4, 'color',[0 0.4470 0.7410])
xlabel('Fraction Of Counts','FontWeight','bold');
ylabel('R Standard Deviation','FontWeight','bold');
legend('Bootstrap','Scaling')
legend('95% Confidence Interval','Experimental Replicates','Bootstrap','Scaling')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','lin');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
ylim([0 20])
% Rg
figure()
plot(frac, Rg_error_ci,'s','MarkerSize',15,'MarkerEdgeColor','k','LineWidth',4)
hold on
plot(frac, exp{:,5},'o','MarkerSize',15,'MarkerEdgeColor','[0.8500 0.3250 0.0980]','LineWidth',4)
hold on
plot(frac, Rg_error_bootstrap,'^','MarkerSize',15,'MarkerEdgeColor','[0 0.4470 0.7410]','LineWidth',4)
hold on

%scale = sqrt(0.01)*Rg_error_bootstrap(1);
scale_Rg = lsqnonlin(@(scale)fitscalingfactor(frac, Rg_error_bootstrap , scale), 0.35); 

plot(frac_fit,1./sqrt(frac_fit)*scale_Rg,'-','LineWidth',4, 'color',[0 0.4470 0.7410])
xlabel('Fraction Of Counts','FontWeight','bold');
ylabel('Rg Standard Deviation','FontWeight','bold');
legend('Bootstrap','Scaling')
legend('95% Confidence Interval','Experimental Replicates','Bootstrap','Scaling')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','lin');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
%ylim([0 6])
%% AIC and BIC

modeldiff = csvread('F-127modeldifferentiation.csv',1,0);
figure()
plot(modeldiff(:,1), modeldiff(:,2),'s','MarkerSize',20,'LineWidth',3, 'Color', 'k')
hold on
plot(modeldiff(:,1), modeldiff(:,3),'^','MarkerSize',20,'LineWidth',3, 'Color','[0.8500 0.3250 0.0980]')
hold on
plot(modeldiff(:,1), modeldiff(:,4),'o','MarkerSize',20,'LineWidth',3, 'Color','[0 0.4470 0.7410]')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
xlabel('Fraction Of Counts','FontWeight','bold');
ylabel('AIC','FontWeight','bold');
lgd = legend('sphere','fuzzy sphere','spherical micelle');
lgd.FontSize = 20;
xlim([0.009 1.0])
ylim([-300 100])
xticks([10^(-2) 10^(-1) 10^0])
xticklabels({'10^{-2}','10^{-1}','1'})
figure()
plot(modeldiff(:,1), modeldiff(:,5),'s','MarkerSize',20,'LineWidth',3, 'Color', 'k')
hold on
plot(modeldiff(:,1), modeldiff(:,6),'^','MarkerSize',20,'LineWidth',3, 'Color', '[0.8500 0.3250 0.0980]')
hold on
plot(modeldiff(:,1), modeldiff(:,7),'o','MarkerSize',20,'LineWidth',3, 'Color', '[0 0.4470 0.7410]')
hold on
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
xlabel('Fraction Of Counts','FontWeight','bold');
ylabel('BIC','FontWeight','bold');
lgd = legend('sphere','fuzzy sphere','spherical micelle');
lgd.FontSize = 20;
xlim([0.009 1.0])
ylim([-300 100])
xticks([10^(-2) 10^(-1) 10^0])
xticklabels({'10^{-2}','10^{-1}','1'})

% function to minimize
function dy = fitscalingfactor(frac, dI, scale)
    dy = (log(1./sqrt(frac)*scale) - log(dI));
end