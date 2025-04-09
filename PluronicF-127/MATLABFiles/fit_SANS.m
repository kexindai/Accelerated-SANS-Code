close all
clc
clear

%% Inputs
% here we choose one replicates for each fraction of counts for bootstrapping
% the parameters are [ndensity; Rc; Rg; d;  A; b; sld_core; sld_shell];
Nsample = 1;
paramsave = zeros(Nsample,16);
k = 11;
filename = sprintf('F127_%d_merged_0p01.txt', k-1);
a_sim = 1;
N_rep = 100;
saveout = 'y';

%% Generate simulated files
% Simulated files are named in this fashion: strcat(name,'_Nsim_',num2str(N_sim(i)),ext)
%simulate_sigma(filename,a_sim);

%% Bootstrap from the simulated files
for i = 1:length(a_sim)
    [~,name,ext] = fileparts(filename);
    sim_file = strcat(name,ext);
    [I_rep,M_data] = bootstrap_SANS(sim_file, N_rep);
    if i == 1
        I_rep_mat = zeros(size(I_rep,1),size(I_rep,2),length(a_sim));
    end
    I_rep_mat(:,:,i) = I_rep;
end

%% Plot the bootstrap and original file
figure()
semilogx(M_data(:,1), I_rep(:,2:end),'o')
hold on
semilogx(M_data(:,1),M_data(:,2),'.','Color','k', 'MarkerSize',20)
xlabel('q [nm^-^1]','FontWeight','bold');
ylabel('I(q) [cm^-^1]','FontWeight','bold');
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2,'xscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);
%% Perform the fit F127
saveparam_spheremicelle = zeros(N_rep,18);
for i = 1:N_rep
    for j = 1:length(a_sim)
        sim_file = strcat(name, ext);
        M_data = importdata(sim_file);
        [saveparam_spheremicelle(i,:), I_fit(:,i), residual(:,i), ci] = micelle_bootstrap(M_data, I_rep_mat(:,i));
    end
end

%% plot the fit for F127
figure()
close all
loglog(M_data(22:end-2,1), I_rep(22:end-2,2:end),'o')
hold on
loglog(M_data(22:end-2,1),M_data(22:end-2,2),'.','Color','k', 'MarkerSize',20)
hold on
loglog(M_data(22:end-2,1), I_fit(:,2:end),'-','Color','g')
hold on
loglog(M_data(22:end-2,1), I_fit(:,1),'-','Color','r','LineWidth',2)
xlabel('q [nm^-^1]','FontWeight','bold');
ylabel('I(q) [cm^-^1]','FontWeight','bold');
set(gca,'FontSize',16,'TickLength',[0.03 0.03],'LineWidth',2,'xscale','log');
set(gcf,'Color','w','units','pixels','outerposition',[50 50 600 600]);

%%
csvwrite('Bootstrap_F127_0p01.csv', saveparam_spheremicelle)
