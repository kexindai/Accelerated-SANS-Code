clc
clear
data = readtable('Rgvsfrac.csv');
data = table2array(data);
UCL = data(end, 2) + data(end, 2) * 0.05;
LCL = data(end, 2) - data(end, 2) * 0.05;
figure()
errorbar(data(:,1),data(:,2),data(:,3),'s','MarkerSize',15,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',4)
hold on
errorbar(data(:,1),data(:,4),data(:,5),'^','MarkerSize',15,'Color','[0 0.4470 0.7410]','MarkerEdgeColor','[0 0.4470 0.7410]','LineWidth',4)
hold on
errorbar(data(:,1),data(:,8),data(:,9),'o','MarkerSize',15,'Color','[0.8500 0.3250 0.0980]','MarkerEdgeColor','[0.8500 0.3250 0.0980]','LineWidth',4)
hold on
plot(data(:,1),UCL*ones(length(data(:,1)),1),'--','LineWidth',3, 'color',[0.4940 0.1840 0.5560])
hold on
plot(data(:,1),LCL*ones(length(data(:,1)),1),'--','LineWidth',3,'color',[0.4940 0.1840 0.5560])
legend('Experiment with 95% CI','Bootstrap','Experimental Average with STD')
xlim([0.008 1])
ylim([25 35])
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
xlabel('Fraction Of Counts','FontWeight','bold');
ylabel('Rg [Å]','FontWeight','bold');



%% the scaling of Rg STD from bootstrapping
figure()
plot(data(:,1),data(:,3),'s','MarkerSize',15,'MarkerEdgeColor','k','LineWidth',4)
hold on
plot(data(:,1), data(:,5),'^','MarkerSize',15,'Color','[0 0.4470 0.7410]','MarkerEdgeColor','[0 0.4470 0.7410]','LineWidth',4)
hold on
plot(data(:,1), data(:,9),'o','MarkerSize',15,'MarkerEdgeColor','[0.8500 0.3250 0.0980]','LineWidth',4)
hold on
scale = 1./sqrt(0.01)/data(1,5);
scale_Rg = lsqnonlin(@(scale)fitscalingfactor(data(:,1), data(:,5) , scale), 1/scale); 
frac_fit = linspace(0.01, 1,1001);
plot(frac_fit,1./sqrt(frac_fit)*scale_Rg,'-','LineWidth',4, 'color','[0 0.4470 0.7410]')
xlabel('Fraction Of Counts','FontWeight','bold');
ylabel('Rg Standard Deviation','FontWeight','bold');
legend('95% Confidence Interval', 'Bootstrap','Experimental Replicates','Scaling')
%legend('MC Bootstrapping','Scaling')
ylim([-0.1 2.5])
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','lin');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
function dy = fitscalingfactor(frac, dI, scale)
    dy = log(1./sqrt(frac)*scale) - log(dI);
end
