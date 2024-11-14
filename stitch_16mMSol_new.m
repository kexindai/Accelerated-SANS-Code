clc
clear
close all
vp = 0.0685;
for i = 1:40
clearvars -except i vp
sample_q2 = importdata(sprintf('16mMSolq_2_%d_1D.txt', i-1));
sample_q1 = importdata(sprintf('16mMSolq_1_%d_1D.txt',i-1));
solvent_q1 = importdata(sprintf('dDMFq_1_1D.txt', i-1));
solvent_q2 = importdata(sprintf('dDMFq_2_1D.txt',i-1));
q1(:,1) = sample_q1.data(:,1);
q1(:,2) =  sample_q1.data(:,2)  - (1-vp)*solvent_q1.data(:,2) ;
q1(:,3) = sqrt(sample_q1.data(:,3).^2  + solvent_q1.data(:,3).^2);
q2(:,1) = sample_q2.data(:,1);
q2(:,2) =  sample_q2.data(:,2)  - (1-vp)*solvent_q2.data(:,2) ;
q2(:,3) = sqrt(sample_q2.data(:,3).^2  + solvent_q2.data(:,3).^2);
%% Merge the data
overlap = [0.033, 0.049];
[~, ind_1_LB] = min(abs(q1(:,1) - overlap(1)));
[~, ind_1_UB] = min(abs(q1(:,1) - overlap(2)));
[~, ind_2_LB] = min(abs(q2(:,1) - overlap(1)));
[~, ind_2_UB] = min(abs(q2(:,1) - overlap(2)));
q1_fit = q1(ind_1_LB:ind_1_UB,:);
q2_fit = q2(ind_2_LB:ind_2_UB,:);
a_scale = lsqnonlin(@ (a) scale(q1_fit,q2_fit,a) , 0.5, 0, 1);
%%
% figure()
% loglog(q1 (:,1),q1 (:,2))
% hold on
% loglog(q2 (:,1),a_scale.*q2 (:,2))
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
writematrix(q,sprintf('16mMSol_%d_merged.txt',i-1),'Delimiter','space')  

%% Plot to check
q = importdata('16mMSol_0_merged.txt');
%%
figure()
loglog (q(:,1)*10,q(:,2),'o')
axis([0.1 1.5 10^(-4) 4.3])
saveas(gcf,sprintf('16mMSol_%d_merged_noerrorbar.png',i-1))
%%
figure()
errorbar (q(:,1),q(:,2),q(:,3),'s','MarkerSize',6,'MarkerEdgeColor','k','LineWidth',1,'Color','k')
xlabel("q [A^-^1]",'FontWeight','bold');
ylabel('I(q) [cm^-^1]','FontWeight','bold');
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',1,'xscale','log',...
    'yscale','log');
saveas(gcf,sprintf('16mMSol_%d_merged.png',i-1))
%%
figure()
plot (q(:,1),q(:,3),'s','MarkerSize',6,'MarkerEdgeColor','k','LineWidth',1,'Color','k')
xlabel("q [A^-^1]",'FontWeight','bold');
ylabel('dI(q) [cm^-^1]','FontWeight','bold');
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',1);
saveas(gcf,sprintf('16mMSol_%d_merged_error.png',i-1'))
end
function ls = scale(q1,q2,a)
    ls = sum(q1(:,2) - a * q2(:,2)).^2; 
end


