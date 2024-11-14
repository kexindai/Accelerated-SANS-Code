close all
clc
clear
data1 = importdata('F127q_2_0_0p05.txt');
data2 = importdata('F127q_2_0_0p1.txt');
data_merged_1 = importdata('F127_0p05_merged.txt')
data_merged_2 = importdata('F127_0p1_merged.txt')
solvent2 = importdata('D2Oq_2_1D.txt');
dI_1 = sqrt(data1.data(:,3).^2 + solvent2.data(:,3).^2);
dI_2 = sqrt(data2.data(:,3).^2 + solvent2.data(:,3).^2);
figure(1)
loglog(data1.data(:,1), dI_1,'o')
hold on
loglog(data2.data(:,1), dI_2,'o')
hold on
% % loglog(data_merged_1(:,1), data_merged_1(:,3),'o')
% % % hold on
% % % loglog(data_merged_2(:,1), data_merged_2(:,3),'o')
legend('error of 0.05 calculated','error of 0.1 calculated','error of 0.05 merged','error of 0.1 merged')
edit = importdata('F127_0_merged.txt');
edit(326:902,3) = dI_2;
i = 1;
%writematrix(edit,sprintf('F127_%d_merged.txt',i-1),'Delimiter','space') 
%%
figure(2)
errorbar(data1.data(:,1),data1.data(:,2),data1.data(:,3),'o')
hold on
errorbar(data2.data(:,1),data2.data(:,2),data2.data(:,3),'o')
legend('0.05','0.1')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);


figure(3)
errorbar(data_merged_1(:,1),data_merged_1(:,2),data_merged_1(:,3),'o')
hold on
errorbar(data_merged_2(:,1),data_merged_2(:,2),data_merged_2(:,3),'o')
legend('0.05','0.1')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);

%%
figure(4)
errorbar(data_merged_1(:,1),data_merged_1(:,2),data_merged_1(:,3),'o')
hold on
errorbar(data1.data(:,1),data1.data(:,2),data1.data(:,3),'o')
legend('merged','raw')
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log','yscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);

%% compare error of 0.05 and 0.1 in intensity
figure(5)
plot(data_merged_1(:,1), data_merged_1(:,3),'o')
hold on
plot(data_merged_2(:,1), data_merged_2(:,3),'o')
legend('0.05','0.1')