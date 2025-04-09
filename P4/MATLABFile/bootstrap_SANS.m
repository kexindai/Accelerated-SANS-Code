% Inputs: (1) filename, string
%         (2) N_rep, number of replicates
% Output: (1) I_rep: a matrix of size N_q x N_rep, where N_q is the number
%             of q-values
%         (2) M_data: original data matrix
function [I_rep,M_data] = bootstrap_SANS(filename,N_rep)
%% Read in the data
% Columns are (1) q [inverse Angstrom] 
%             (2) I [inverse cm]
%             (3) dI [inverse cm] >> this is "sigma"
%             (4) dq [inverse Angstrom]
M_data = importdata(filename);
I = M_data(:,2);
dI = M_data(:,3);
I_rep = zeros(size(M_data,1),N_rep);
% this is the parameters for gamma distribution
a = I.^2./dI.^2;
b = dI.^2./I;
rng("default")
%% Create replicates
% Populate the M = 1 slice using the original data
I_rep(:,1) = M_data(:,2);
% Create replicas that are N(Icut,SDcut^2). Note, MATLAB reads the Gaussian
% distribution as N(Icut,SDcut), so it uses the SD, not the variance.
% Obviously, do not create replicas if M = 1.
if N_rep > 1
    for i = 2:N_rep
        % change to gamma distribution. the least square needs to change
        % too if we do gamma
        %I_rep(:,i) = gamrnd(a,b, size(I));
        I_rep(:,i) = normrnd(I,dI,size(I));
    end
end