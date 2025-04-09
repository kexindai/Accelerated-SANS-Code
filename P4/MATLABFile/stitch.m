clc
clear
close all
%%
%addpath(genpath('CD_2configTimeSlicing_P4_500pts_1_mergedQ'))
%% import data
for i = 1:100
P4_q2 = importdata(sprintf('P4q_2_%d_1D.txt', i-1));
P4_q1 = importdata(sprintf('P4q_1_%d_1D.txt',i-1));
%% Merge the data
overlap = [0.033, 0.049];
[~, ind_1_LB] = min(abs(P4_q1.data(:,1) - overlap(1)));
[~, ind_1_UB] = min(abs(P4_q1.data(:,1) - overlap(2)));
[~, ind_2_LB] = min(abs(P4_q2.data(:,1) - overlap(1)));
[~, ind_2_UB] = min(abs(P4_q2.data(:,1) - overlap(2)));
q1 = P4_q1.data(ind_1_LB:ind_1_UB,:);
q2 = P4_q2.data(ind_2_LB:ind_2_UB,:);
a_scale = lsqnonlin(@ (a) scale(q1,q2,a) , 0.5, 0, 1);
%%
% figure()
% loglog(P4_q1.data(:,1),P4_q1.data(:,2))
% hold on
% loglog(P4_q2.data(:,1),a_scale.*P4_q2.data(:,2))
%%
qcut_1 = [P4_q1.data(1,1),P4_q1.data(ind_1_UB,1)];
qcut_2 = [P4_q2.data(ind_2_LB,1), P4_q2.data(end,1)];

p4_cut_1(:,1) = P4_q1.data(P4_q1.data(:,1)>qcut_1(1) & P4_q1.data(:,1)<qcut_1(2),1); % q [A^-1]
p4_cut_1(:,2) = P4_q1.data(P4_q1.data(:,1)>qcut_1(1) & P4_q1.data(:,1)<qcut_1(2),2); % I [cm^-1]
p4_cut_1(:,3) = P4_q1.data(P4_q1.data(:,1)>qcut_1(1) & P4_q1.data(:,1)<qcut_1(2),3); % dI [cm^-1]
p4_cut_1(:,4) = P4_q1.data(P4_q1.data(:,1)>qcut_1(1) & P4_q1.data(:,1)<qcut_1(2),4); % dq [A^-1]
p4_cut_2(:,1)= P4_q2.data(P4_q2.data(:,1)>qcut_2(1) & P4_q2.data(:,1)<qcut_2(2),1); % q [A^-1]
p4_cut_2(:,2) = P4_q2.data(P4_q2.data(:,1)>qcut_2(1) & P4_q2.data(:,1)<qcut_2(2),2); % I [cm^-1]
p4_cut_2(:,3) = P4_q2.data(P4_q2.data(:,1)>qcut_2(1) & P4_q2.data(:,1)<qcut_2(2),3); % dI [cm^-1]
p4_cut_2(:,4) = P4_q2.data(P4_q2.data(:,1)>qcut_2(1) & P4_q2.data(:,1)<qcut_2(2),4); % dq [A^-1]
p4_cut_2(:,2) = p4_cut_2(:,2)*a_scale;
P4 = [p4_cut_1; p4_cut_2];
%%
writematrix(P4,sprintf('P4_%d_merged.txt',i-1),'Delimiter','space')  

% %% Plot the merged file and the error
errorbar (P4(:,1),P4(:,2),P4(:,3),'s','MarkerSize',6,'MarkerEdgeColor','k','LineWidth',1,'Color','k')
xlabel("q [A^-^1]",'FontWeight','bold');
ylabel('I(q) [cm^-^1]','FontWeight','bold');
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',1,'xscale','log',...
    'yscale','log');
saveas(gcf,sprintf('P4_%d_merged_errorbar.png',i-1'))
plot (P4(:,1),P4(:,3),'s','MarkerSize',6,'MarkerEdgeColor','k','LineWidth',1,'Color','k')
xlabel("q [A^-^1]",'FontWeight','bold');
ylabel('dI(q) [cm^-^1]','FontWeight','bold');
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',1);
saveas(gcf,sprintf('P4_%d_merged_error.png',i-1'))
end
function ls = scale(q1,q2,a)
    ls = sum(q1(:,2) - a * q2(:,2)).^2; 
end


