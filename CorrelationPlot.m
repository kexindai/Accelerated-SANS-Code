clc
clear 
close all
%% Code testing for random process
%%% note that the image is not different from the actual data!
rng('default')
X = randn(563,100);
Corr_rand = corr(X');
figure(1)
imagesc(Corr_rand);
colorbar
set(gca,'XScale','lin', 'YScale','lin','YMinorTick','on','Ydir','normal',...
    'FontSize', 20, 'LineWidth',2) 
% initialize
Number_Rep = 100;
sample_qlow_I_data = zeros(563, Number_Rep);
sample_qhigh_I_data = zeros(579, Number_Rep);
sample_qlow_dI_data = zeros(563, Number_Rep);
sample_qhigh_dI_data = zeros(579, Number_Rep);
% Here we read all the data and store them
for i = 1:Number_Rep
sample_qlow = importdata(sprintf('P4q_1_%d_1D.txt',i-1));
sample_qlow_I_data(:, i) = sample_qlow.data(:,2);
sample_qlow_dI_data(:, i) = sample_qlow.data(:,3);

sample_qhigh = importdata(sprintf('P4q_2_%d_1D.txt', i-1));
sample_qhigh_I_data(:, i) = sample_qhigh.data(:,2);
sample_qhigh_dI_data(:, i) = sample_qhigh.data(:,3);
end

%% Here we need to first find the average. Separate by low q and high q here
%%% for low q
q_low = sample_qlow.data(:, 1); % extra the q
average_lowq = mean(sample_qlow_I_data,2);

var_lowq = std(sample_qlow_I_data, 0, 2).^2; % calculate the variance
diff_I_lowq = sample_qlow_I_data - repmat(average_lowq, 1, 100);

Correlation_lowq = corr(diff_I_lowq');
Corre_func = zeros(length(diff_I_lowq),1);


v = zeros(length(diff_I_lowq), length(diff_I_lowq));

% this calculates the dot product of the residual across all the replicates
% for each q. the index j prevents double counting
for i = 1:563
    for j = 1:563-(i-1)
        v(j, i) = dot(diff_I_lowq(j,:), diff_I_lowq(j+i-1,:));
    end
end


% this normalize the correlation function by total variance across all q
% values and replicates
for i = 1:563
    Corre_func (i,1) = sum(v(:,i))/sum(v(:,1));
    
end

off_diagonal_correlation = zeros(length(diff_I_lowq), 1);
for i = 1:563
    for j = 1:563-i+1
        off_diagonal_correlation(j) = off_diagonal_correlation(j) + Correlation_lowq(i+j-1,i);
    end
end
for i = 1:563
    off_diagonal_correlation(i) = off_diagonal_correlation(i) / (563 - i + 1);
end
figure(1)
plot(q_low - q_low(1), Corre_func,'LineWidth',2)
set(gca,'XScale','lin', 'YScale','lin','YMinorTick','on','Ydir','normal',...
    'FontSize', 20, 'LineWidth',2) 
xlabel('\Delta q [Å]')
ylabel('Correlation Function')
figure(2)
plot(q_low, off_diagonal_correlation,'LineWidth',2)
set(gca,'XScale','lin', 'YScale','lin','YMinorTick','on','Ydir','normal',...
    'FontSize', 20, 'LineWidth',2) 
xlabel('q [Å]')
ylabel('Correlation Function')

figure(3)
imagesc(q_low, q_low, Correlation_lowq)
colorbar
xlabel('q [Å]')
ylabel('q [Å]')
% xlim([0.0037 0.005])
% ylim([0.0037 0.005])
set(gca,'XScale','log', 'YScale','log','YMinorTick','on','Ydir','normal',...
    'FontSize', 20, 'LineWidth',2) 
%%% for high q
q_high = sample_qhigh.data(:, 1);
average_highq = sum(sample_qhigh_I_data(:,1:Number_Rep), 2)/Number_Rep;
diff_I_highq = sample_qhigh_I_data - repmat(average_highq, 1, 100);
Correlation_highq = corr(diff_I_highq');

Corre_func = zeros(length(diff_I_highq),1);


v = zeros(length(diff_I_highq), length(diff_I_highq));

% this calculates the dot product of the residual across all the replicates
% for each q. the index j prevents double counting
for i = 1:length(diff_I_highq)
    for j = 1:length(diff_I_highq)-(i-1)
        v_high(j, i) = dot(diff_I_highq(j,:), diff_I_highq(j+i-1,:));
    end
end


% this normalize the correlation function by total variance across all q
% values and replicates
for i = 1:length(diff_I_highq)
    Corre_func_high (i,1) = sum(v_high(:,i))/sum(v_high(:,1));
    
end
figure(4)
plot(q_high - q_high(1), Corre_func_high,'LineWidth',2)
set(gca,'XScale','lin', 'YScale','lin','YMinorTick','on','Ydir','normal',...
    'FontSize', 20, 'LineWidth',2) 
xlabel('\Delta q [Å]')
ylabel('Correlation Function')
%clims = [0 10];
figure(5)
imagesc(q_high, q_high, Correlation_highq)
colorbar
xlabel('q [Å]')
ylabel('q [Å]')
set(gca,'XScale','log', 'YScale','log','YMinorTick','on','Ydir','normal',...
    'FontSize', 20, 'LineWidth',2) 
% xlim([0.03 0.05])
% ylim([0.03 0.05])