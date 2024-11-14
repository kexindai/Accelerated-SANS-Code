clc
clear
%% 16mMSol
Sol = importdata('16mMSol_merged.txt')
Sol_0p5 = importdata('16mMSol_0_merged_0p5.txt')
figure(1)
%errorbar(Sol.data(:,1),Sol.data(:,2),Sol.data(:,3),'o')
loglog(Sol.data(:,1),Sol.data(:,2),'o')
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',1,'xscale','log',...
    'yscale','log');
hold on
%errorbar(Sol_0p5(:,1),Sol_0p5(:,2),Sol_0p5(:,3),'o')
loglog(Sol_0p5(:,1),Sol_0p5(:,2),'.')
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',1,'xscale','log',...
    'yscale','log');
xlabel("q [A^-^1]",'FontWeight','bold');
ylabel('I(q) [cm^-^1]','FontWeight','bold');
figure()
loglog(Sol.data(:,1),Sol.data(:,3),'o')
hold on
loglog(Sol_0p5(:,1),Sol_0p5(:,3),'o')
%% P4
figure(3)
P4 = importdata('P4_merged.txt');
loglog(P4.data(:,1),P4.data(:,2),'o')
figure(4)
loglog(P4.data(:,1),P4.data(:,3),'o')
%% F127
clc
clear
figure(5)
F127 = importdata('F127_merged.txt');
loglog(F127.data(:,1),F127.data(:,2),'o')
figure(6)
loglog(F127.data(:,1),F127.data(:,3),'o')