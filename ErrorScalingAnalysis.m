%% this matlab code will plot the error from 1D SANS data as each fraction of counts
clc
clear
close all
lowq = readmatrix('P4_Error_lowQ.csv');
midq = readmatrix('P4_Error_midQ.csv');
highq = readmatrix('P4_Error_highQ.csv');
%%% this reads the low Q error q < 0.033 choose q = 1.506607E-02 to avoid
%%% negative I values at low Q looked at raw file

figure(1)
plot(lowq(:,1),lowq(:,3),'o','MarkerSize',10,'Color','k','MarkerFaceColor','k','LineWidth',2)
hold on
rescale_low = 1./sqrt(lowq(1,1))/lowq(1,3);
x_low = lsqnonlin(@(scale)fitscalingfactor(lowq(:,1),lowq(:,3) , scale), 1/rescale_low); 
frac = logspace(-2,0,50);
plot(frac,1./sqrt(frac)*x_low,'-','LineWidth',2, 'color',[0 0.4470 0.7410])
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2,'xscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
legend('low Q experimental error','low Q scaling')
xlabel('fraction of counts','FontWeight','bold');
ylabel('dI','FontWeight','bold');
%%% this reads the mid Q error q = [0.033, 0.049]; q = 0.04 looked at
%%% merged file


figure(2)
plot(midq(:,1),midq(:,3),'o','MarkerSize',10,'Color','k','MarkerFaceColor','k','LineWidth',2)
hold on
rescale_mid = 1./sqrt(midq(1,1))/midq(1,3);
% now instead we fit for the scaling factor
x_mid = lsqnonlin(@(scale)fitscalingfactor(midq(:,1),midq(:,3) , scale), 1/rescale_mid); 
plot(frac,1./sqrt(frac)*x_mid,'-','LineWidth',2, 'color',[0.8500 0.3250 0.0980])
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2,'xscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
legend('mid Q experimental error','mid Q scaling')
xlabel('fraction of counts','FontWeight','bold');
ylabel('dI','FontWeight','bold');


%%% this reads the high Q error q > 0.049 q = 2.500345E-01	looked at raw
%%% file

figure(3)

plot(highq(:,1),highq(:,3),'o','MarkerSize',10,'Color','k','MarkerFaceColor','k','LineWidth',2)
hold on
rescale_high = 1./sqrt(highq(1,1))/highq(1,3);
x_high = lsqnonlin(@(scale)fitscalingfactor(highq(:,1),highq(:,3) , scale), 1/rescale_high); 
plot(frac,1./sqrt(frac)*x_high,'-','LineWidth',2, 'color',[0.9290 0.6940 0.1250])

set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2,'xscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
legend('high Q experimental error','high Q scaling')
xlabel('fraction of counts','FontWeight','bold');
ylabel('dI','FontWeight','bold');
%% plot everything on the same plot
figure()
plot(lowq(:,1),lowq(:,3),'o','MarkerSize',10,'Color',[0 0.4470 0.7410],'MarkerFaceColor',[0 0.4470 0.7410],'LineWidth',2)
hold on
% this is to rescale everything relative to 0.01
% now we try to fit for rescale_low instead

plot(frac,1./sqrt(frac)*x_low,'-','LineWidth',2, 'color',[0 0.4470 0.7410])
hold on
plot(midq(:,1),midq(:,3),'o','MarkerSize',10,'Color',[0.8500 0.3250 0.0980],'MarkerFaceColor',[0.8500 0.3250 0.0980],'LineWidth',2)
hold on
plot(frac,1./sqrt(frac)*x_mid,'-','LineWidth',2, 'color',[0.8500 0.3250 0.0980])
hold on
plot(highq(:,1),highq(:,3),'o','MarkerSize',10,'Color',[0.9290 0.6940 0.1250],'MarkerFaceColor',[0.9290 0.6940 0.1250],'LineWidth',2)
hold on
plot(frac,1./sqrt(frac)*x_high,'-','LineWidth',2, 'color',[0.9290 0.6940 0.1250])
hold on
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2,'xscale','log','yscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
legend('low Q experimental error','low Q scaling','mid Q experimental error','mid Q scaling', 'high Q experimental error',...
    'high Q scaling')
xlabel('fraction of counts','FontWeight','bold');
ylabel('dI','FontWeight','bold');
ylim([0 10])

function dy = fitscalingfactor(frac, dI, scale)
    dy = log(1./sqrt(frac)*scale) - log(dI);
end