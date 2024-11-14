clc
clear
close all
%% import data
rho_PEO = 1.125;
rho_PPO = 1.06;
N_PEO = 65;
N_PPO = 198;
Mw_EO = 44.05;
Mw_PO = 58;
phi_PEO = (N_PEO*Mw_EO/rho_PEO)/(N_PEO*Mw_EO/rho_PEO + N_PPO* Mw_PO/rho_PPO);
rho = phi_PEO*rho_PEO + (1 - phi_PEO) *rho_PPO;
vp = rho*1.5/100; %1.5 w/v %
for i = 1:1
clearvars -except i vp
sample_q1 = importdata(sprintf('F127q_1_%d_1D.txt',i-1));
sample_q2 = importdata(sprintf('F127q_2_%d_1D.txt', i-1));
solvent_q1 = importdata(sprintf('D2Oq_1_1D.txt'));
solvent_q2 = importdata(sprintf('D2Oq_2_1D.txt'));
q1(:,1) = sample_q1.data(length(sample_q1.data) - length(solvent_q1.data) +1:end,1);
q1(:,2) =  sample_q1.data(length(sample_q1.data) - length(solvent_q1.data) +1:end,2)  - (1-vp)*solvent_q1.data(:,2) ;
q1(:,3) = sqrt(sample_q1.data(length(sample_q1.data) - length(solvent_q1.data) +1:end,3).^2  + solvent_q1.data(:,3).^2);
q2(:,1) = sample_q2.data(1:length(solvent_q2.data),1);
q2(:,2) =  sample_q2.data(1:length(solvent_q2.data),2)  - (1-vp)*solvent_q2.data(:,2) ;
% I think something is wrong with this. 0.05 and 0.1 are weird
q2(:,3) = sqrt(sample_q2.data(1:length(solvent_q2.data),3).^2  + solvent_q2.data(:,3).^2);
%% Merge the data
overlap = [0.033, 0.049];
[~, ind_1_LB] = min(abs(q1(:,1) - overlap(1)));
[~, ind_1_UB] = min(abs(q1(:,1) - overlap(2)));
[~, ind_2_LB] = min(abs(q2(:,1) - overlap(1)));
[~, ind_2_UB] = min(abs(q2(:,1) - overlap(2)));
q1_fit = q1(ind_1_LB:ind_1_LB + length(ind_2_LB:ind_2_UB) - 1,:);
q2_fit = q2(ind_2_LB:ind_2_UB,:);
a_scale = lsqnonlin(@ (a) scale(q1_fit,q2_fit,a) , 0.5, 0, 1);
%%
% figure()
% loglog(q1.data(:,1),q1.data(:,2))
% hold on
% loglog(q2.data(:,1),a_scale.*q2.data(:,2))
%%
qcut_1 = [q1(212,1),q1(ind_1_UB,1)];
qcut_2 = [q2(ind_2_LB,1), q2(end,1)];

q_cut_1(:,1) = q1 (q1 (:,1)>qcut_1(1) & q1 (:,1)<qcut_1(2),1); % q [A^-1]
q_cut_1(:,2) = q1 (q1 (:,1)>qcut_1(1) & q1 (:,1)<qcut_1(2),2); % I [cm^-1]
q_cut_1(:,3) = q1 (q1 (:,1)>qcut_1(1) & q1 (:,1)<qcut_1(2),3); % dI [cm^-1]
%q_cut_1(:,4) = q1 (q1 (:,1)>qcut_1(1) & q1 (:,1)<qcut_1(2),4); % dq [A^-1]
q_cut_2(:,1) = q2 (q2 (:,1)>qcut_2(1) & q2 (:,1)<qcut_2(2),1); % q [A^-1]
q_cut_2(:,2) = q2 (q2 (:,1)>qcut_2(1) & q2 (:,1)<qcut_2(2),2); % I [cm^-1]
q_cut_2(:,3) = q2 (q2 (:,1)>qcut_2(1) & q2 (:,1)<qcut_2(2),3); % dI [cm^-1]
%q_cut_2(:,4) = q2 (q2 (:,1)>qcut_2(1) & q2 (:,1)<qcut_2(2),4); % dq [A^-1]
q_cut_2(:,2) = q_cut_2(:,2)*a_scale;

q = [q_cut_1; q_cut_2];
%%
writematrix(q,sprintf('F127_%d_merged.txt',i-1),'Delimiter','space')  

%% Plot to check
figure('visible','off')
loglog (q(:,1),q(:,2),'o')
saveas(gcf,sprintf('F127_%d_merged_noerrorbar.png',i-1))
figure('visible','off')
errorbar (q(:,1),q(:,2),q(:,3),'s','MarkerSize',6,'MarkerEdgeColor','k','LineWidth',1,'Color','k')
xlabel("q [A^-^1]",'FontWeight','bold');
ylabel('I(q) [cm^-^1]','FontWeight','bold');
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',1,'xscale','log',...
    'yscale','log');
saveas(gcf,sprintf('F127_%d_merged.png',i-1))
figure('visible','off')
plot (q(:,1),q(:,3),'s','MarkerSize',6,'MarkerEdgeColor','k','LineWidth',1,'Color','k')
xlabel("q [A^-^1]",'FontWeight','bold');
ylabel('dI(q) [cm^-^1]','FontWeight','bold');
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',1);
saveas(gcf,sprintf('F127_%d_merged_error.png',i-1'))
end
function ls = scale(q1,q2,a)
    ls = sum(q1(:,2) - a * q2(:,2)).^2; 
end


