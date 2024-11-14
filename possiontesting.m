clc
clear
Nsample = 1;
paramsave = zeros(Nsample,16);
k=1;
filename = sprintf('P4_%d_merged.txt', k-1);
a_sim = 1; % this means we are not rescaling the error. Use as it is in the raw data.
%a_sim = [0.001 0.0025 0.005 0.01 0.025 0.05 0.1 0.25 0.5 1];
N_rep = 100;
saveout = 'y';
phi = 0.05;
data = importdata(filename);
pd = makedist('Poisson','lambda',data(1,2));
pd_n = makedist('Normal','mu', data(1,2) , 'sigma', data(1,3))
x = -10:1:100;
p = pdf(pd,x);
p_n = pdf(pd_n, x)
figure (1)
plot(x,p)
hold on
plot(x,p_n)
legend('Poisson','Normal')

poisson_val = poissrnd(data(1,2),1,100);
normal_val = normrnd(data(1,2),data(1,3),1,100);
%%
Nsample = 1;
paramsave = zeros(Nsample,16);
k=1;
    % change folder to get bootstrapped error
filename = sprintf('P4_%d_merged.txt', k-1);
%N_input = 1000; % number of counts for the input file starting point
% a_sim = linspace(0.1,1,10); % number of counts to simulate
% a_sim(end+1) = 0.001;
% a_sim = logspace(-3,0,10);
a_sim = 1; % this means we are not rescaling the error. Use as it is in the raw data.
%a_sim = [0.001 0.0025 0.005 0.01 0.025 0.05 0.1 0.25 0.5 1];
N_rep = 100;
saveout = 'y';
phi = 0.05;
data = importdata(filename);
x = randi([0 100],1,100);
y = poisspdf(x,data(1,2));

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