%% this compares the error of CI and bootstrap xi and q0
clc
clear
close all
param_sim = csvread('BootstrapError.csv',1,0);
correlation = readmatrix('correlation.csv');
UCL_q0 = correlation(end,10) + correlation(end,10) * 0.05;
LCL_q0 = correlation(end,10) - correlation(end,10) * 0.05;
UCL_xi = correlation(end,12) + correlation(end,12) * 0.05;
LCL_xi = correlation(end,12) - correlation(end,12) * 0.05;
%param_sim = table2array(data);
frac = param_sim(:,1);
figure()
% this plots q0
errorbar(frac,correlation(:,10)*10,correlation(:,11)*10,'s','MarkerSize',15,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',4);
hold on
errorbar(frac,correlation(:,23)*10,correlation(:,24)*10,'o','MarkerSize',15,'Color',[0.8500 0.3250 0.0980],'MarkerEdgeColor',[0.8500 0.3250 0.0980],'LineWidth',4);
hold on
errorbar(frac,param_sim(:,5+1)*10,param_sim(:,13+1)*10,'^','MarkerSize',15,'Color',[0 0.4470 0.7410],'LineWidth',4);
hold on
plot(frac,UCL_q0*ones(length(frac))*10,'--','LineWidth',3, 'Color', [0.4940 0.1840 0.5560])
hold on
plot(frac,LCL_q0*ones(length(frac))*10,'--','LineWidth',3,'Color', [0.4940 0.1840 0.5560])
hold on
lgd_1 = legend('Experiment with 95% CI','Experimental Average with STD','Bootstrap');
lgd_1.FontSize = 20;
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
xlabel('Fraction Of Counts','FontWeight','bold');
ylabel('q_{0} [nm^{-1}]','FontWeight','bold');
xlim([0.008 1])
ylim([0 0.05*10])

% this plots xi. need to figure out why the bootstrap average and std are biased!
figure()
errorbar(frac,correlation(:,12)/10,correlation(:,13)/10,'s','MarkerSize',15,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',4);
hold on
errorbar(frac,correlation(:,21)/10,correlation(:,22)/10,'o','MarkerSize',15,'Color',[0.8500 0.3250 0.0980],'MarkerEdgeColor',[0.8500 0.3250 0.0980],'LineWidth',4);
hold on
errorbar(frac,param_sim(:,7)/10,param_sim(:,15)/10,'^','MarkerSize',15,'Color',[0 0.4470 0.7410],'LineWidth',4);
hold on
plot(frac,UCL_xi*ones(length(frac))/10,'--','LineWidth',3, 'Color', [0.4940 0.1840 0.5560])
hold on
plot(frac,LCL_xi*ones(length(frac))/10,'--','LineWidth',3, 'Color', [0.4940 0.1840 0.5560])
hold on
lgd_2 = legend('Experiment with 95% CI','Experimental Average with STD','Bootstrap')
lgd_2.FontSize = 20;
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
xlabel('Fraction Of Counts','FontWeight','bold');
ylabel('\xi [nm]','FontWeight','bold');
xlim([0.008 1])
ylim([2 4.5])

%% this section plots the parameter error scaling
% xi
figure()
plot(frac, correlation(:,13)/10, 's','MarkerSize',15,'MarkerEdgeColor','k', 'LineWidth',4)
hold on
plot(frac, correlation(:,22)/10,'o','MarkerSize',15,'MarkerEdgeColor','[0.8500 0.3250 0.0980]', 'LineWidth',4)
hold on
plot(frac, param_sim(:,15)/10,'^','MarkerSize',15,'MarkerEdgeColor','[0 0.4470 0.7410]', 'LineWidth',4)
hold on
scale = 1./sqrt(0.01)/param_sim(1,15)/10;
scale_xi = lsqnonlin(@(scale)fitscalingfactor(frac, param_sim(:,15)/10 , scale), 1/scale); 
frac_fit = linspace(0.01, 1,1001);
plot(frac_fit,1./sqrt(frac_fit)*scale_xi,'-','LineWidth',4, 'color',[0 0.4470 0.7410])
xlabel('Fraction Of Counts','FontWeight','bold');
ylabel('\xi Standard Deviation','FontWeight','bold');
legend('95% Confidence Interval','Experimental Replicates','Bootstrap','Scaling')
%legend('MC Bootstrapping','Scaling')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','lin');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);

%% q0
figure()
plot(frac, correlation(:,11)*10,'s','MarkerSize',15,'MarkerEdgeColor','k','LineWidth',4)
hold on
plot(frac, correlation(:,24)*10,'o','MarkerSize',15,'MarkerEdgeColor','[0.8500 0.3250 0.0980]','LineWidth',4)
hold on
plot(frac, param_sim(:,13+1)*10,'^','MarkerSize',15,'MarkerEdgeColor','[0 0.4470 0.7410]','LineWidth',4)
hold on
scale = 1./sqrt(0.01)/param_sim(1,13+1)*10;
scale_q0 = lsqnonlin(@(scale)fitscalingfactor(frac, param_sim(:,13+1)*10 , scale), 1/scale); 
frac_fit = linspace(0.01, 1,1001);
plot(frac_fit,1./sqrt(frac_fit)*scale_q0,'-','LineWidth',4, 'color',[0 0.4470 0.7410])
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','lin');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
xlabel('Fraction Of Counts','FontWeight','bold');
ylabel('q_{0} Standard Deviation','FontWeight','bold');
legend('95% Confidence Interval','Experimental Replicates','Bootstrap','Scaling')
%legend('MC Bootstrapping','Scaling')


function dy = fitscalingfactor(frac, dI, scale)
    dy = log(1./sqrt(frac)*scale) - log(dI);
end