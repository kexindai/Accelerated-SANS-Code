% Inputs: (1) filename, string
%         (2) a_sim: multiplicative factor of number of counts to simulate (vector)
%         (3) saveout: 'y','n' to save the output or not
function simulate_sigma(filename,a_sim)
[~,name,ext] = fileparts(filename);
%% Simulate the set of uncertainties
% We know that sigma ~ 1/sqrt(N), where N = number of counts, so if we
% assume that the input data had N_input counts, and the new data has N_new
% counts, where N_new = aN_input, we can calculate the new uncertainties by 
% multiplying by a factor of sqrt(a)

% Columns are (1) q [inverse Angstrom] 
%             (2) I [inverse cm]
%             (3) dI [inverse cm] >> this is "sigma"
%             (4) dq [inverse Angstrom]
M = importdata(filename);
for i = 1:length(a_sim)
    M_sim = zeros(size(M));
    M_sim(:,1:2) = M(:,1:2);
    %M_sim(:,4) = M(:,4);
    M_sim(:,3) = M(:,3)/sqrt(a_sim(i));
    csvwrite(strcat(name,'_asim_',num2str(a_sim(i)),ext),M_sim);
end