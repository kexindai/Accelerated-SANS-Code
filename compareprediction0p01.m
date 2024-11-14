close all
clc
clear
data_exp = importdata('16mMSol_0_merged.txt');
data_sim = importdata('16mMSol_0_merged_asim_100.txt');
figure()
h = gobjects(2, 1);
h(1) = errorbar(data_sim(:,1)*10, data_sim(:,2),data_exp(:,3),'o','MarkerSize',10,'Color','[0 0.4470 0.7410]','MarkerEdgeColor','[0 0.4470 0.7410]','LineWidth',2);
hold on
h(2) = errorbar(data_exp(:,1)*10, data_exp(:,2),data_exp(:,3),'o','MarkerSize',10,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2);
xlim([10^(-1) 4.3])
ylim([10^(-4) 1.5])
xlabel('q [nm^-^1]','FontWeight','bold');
ylabel('I(q) [cm^-^1]','FontWeight','bold');
legend(h([2 1]),{'Experimental Data','Simulated Data'})
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
figure()
loglog(data_exp(:,1)*10, data_exp(:,3),'o','MarkerSize',10,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2)
hold on
loglog(data_sim(:,1)*10, data_sim(:,3),'o','MarkerSize',10,'Color','[0 0.4470 0.7410]','MarkerEdgeColor','[0 0.4470 0.7410]','LineWidth',2)
xlim([10^(-1) 4.3])
ylim([10^(-4) 1.5])
xlabel('q [nm^-^1]','FontWeight','bold');
ylabel('I(q) [cm^-^1]','FontWeight','bold');
legend('Experimental Error','Simulated Error')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);