clc
clear
F127_new = importdata('F127_0_merged.txt');
F127_old = importdata('Pluronics F-127_Solvent_Reduced_500.txt');
%errorbar(F127_new.data(:,1),F127_new.data(:,2),F127_new.data(:,3),'o')
errorbar(F127_new(:,1),F127_new(:,2),F127_new(:,3),'o')
hold on
errorbar(F127_old.data(:,1),F127_old.data(:,2),F127_old.data(:,3),'o')
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',1,'xscale','log',...
    'yscale','log');
