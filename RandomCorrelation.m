clc
clear 
close all
% Generate a random matrix with normally distributed entries
rng(0); % For reproducibility
n = 100; % Number of samples
m = 563;   % Number of variables
data = randn(n, m);

% Calculate the correlation matrix
corr_matrix = corrcoef(data);

% Display the correlation matrix
disp('Correlation matrix:');
disp(corr_matrix);

% Visualize the correlation matrix
figure;
imagesc(corr_matrix);
colorbar;
xlabel('Variables');
ylabel('Variables');
title('Correlation Matrix of Randomly Generated Data');
