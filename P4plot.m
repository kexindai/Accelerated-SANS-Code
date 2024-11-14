clc
clear
close all
correlation = csvread('correlation.csv',1,0);
UCL_q0 = correlation(end,10) + correlation(end,10) * 0.05;
LCL_q0 = correlation(end,10) - correlation(end,10) * 0.05;
UCL_xi = correlation(end,12) + correlation(end,12) * 0.05;
LCL_xi = correlation(end,12) - correlation(end,12) * 0.05;
FSGEL = csvread('FSGEL.csv',1,0);
GLGEL = csvread('GLGEL.csv',1,0);
OZLorentz = csvread('OZLorentz.csv',1,0);
DB = csvread('DB.csv',1,0);
frac = correlation(:,1);
AIC_correlation = correlation(:,19);
BIC_correlation = correlation(:,20);
AIC_FSGEL = FSGEL(:,end-1);
BIC_FSGEL = FSGEL(:,end);
AIC_GLGEL = GLGEL(:,end-1);
BIC_GLGEL = GLGEL(:,end);
AIC_DB = DB(:,end-1);
BIC_DB = DB(:,end);
AIC_OZLorentz = OZLorentz(:,end-1);
BIC_OZLorentz = OZLorentz(:,end-1);
%% this plots the correlation length and q0
figure()
errorbar(frac,correlation(:,10)*10,correlation(:,11)*10,'s','MarkerSize',10,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2);
hold on
plot(frac,UCL_q0*ones(length(frac))*10,'--','LineWidth',2)
hold on
plot(frac,LCL_q0*ones(length(frac))*10,'--','LineWidth',2)
hold on
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
xlabel('Fraction Of Counts','FontWeight','bold');
ylabel('q_{0} [nm^{-1}]','FontWeight','bold');
xlim([0.009 1])
ylim([0 0.05*10])
figure()
errorbar(frac,correlation(:,12)/10,correlation(:,13)/10,'s','MarkerSize',10,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k','LineWidth',2);
hold on
plot(frac,UCL_xi*ones(length(frac))/10,'--','LineWidth',2)
hold on
plot(frac,LCL_xi*ones(length(frac))/10,'--','LineWidth',2)
hold
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
xlabel('Fraction Of Counts','FontWeight','bold');
ylabel('\xi [nm]','FontWeight','bold');
xlim([0.009 1])
ylim([2 3.5])


%% Plot AIC and BIC
figure()
plot(frac,AIC_correlation,'p','MarkerSize',20, 'LineWidth',3)
hold on
plot(frac,AIC_FSGEL,'^','MarkerSize',20,'LineWidth',3)
hold on
plot(frac,AIC_GLGEL,'s','MarkerSize',20,'LineWidth',3)
hold on
plot(frac,AIC_DB,'x','MarkerSize',20,'LineWidth',3)
hold on
plot(frac,AIC_OZLorentz,'o','MarkerSize',20,'LineWidth',3)
hold on
lgd = legend('Broad Peak Model', 'Fine Scale Polymer Gel Model', 'Gauss Lorentz Gel model', 'Debye Bueche','Ornstein Zernike and squared Lorentz');
lgd.FontSize = 20;
xlabel('Fraction Of Counts','FontWeight','bold');
ylabel('AIC','FontWeight','bold');
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);

figure()
plot(frac,BIC_correlation,'p','MarkerSize',20, 'LineWidth',3)
hold on
plot(frac,BIC_FSGEL,'^','MarkerSize',20, 'LineWidth',3)
hold on
plot(frac,BIC_GLGEL,'s','MarkerSize',20, 'LineWidth',3)
hold on
plot(frac,BIC_DB,'x','MarkerSize',20, 'LineWidth',3)
hold on
plot(frac,BIC_OZLorentz,'o','MarkerSize',20, 'LineWidth',3)
hold on
lgd_2 = legend('Broad Peak Model', 'Fine Scale Polymer Gel Model', 'Gauss Lorentz Gel model', 'Debye Bueche','Ornstein Zernike and squared Lorentz');
lgd_2.FontSize = 20;
xlabel('Fraction Of Counts','FontWeight','bold');
ylabel('BIC','FontWeight','bold');
set(gca,'FontSize',30,'TickLength',[0.03 0.03],'LineWidth',3,'xscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);