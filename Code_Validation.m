% Testing
%%
clc
clear
covariance_matrix_of_gaussian_process = zeros(200, 200);
for i = 1:200
  for j = 1:200
      % kernel
    covariance_matrix_of_gaussian_process(i,j) = exp(-(i-j)^2 / 25);
  end
end
mean_of_gaussian_process = zeros(200,1);
replicates = mvnrnd(mean_of_gaussian_process, covariance_matrix_of_gaussian_process, 200);
% Input matrix, specified as an n-by-k matrix. The rows of X correspond to observations, and the columns correspond to variables.

Corre_matrix = corr(replicates);
figure(1)
imagesc(Corre_matrix)
colorbar
figure(2)
imagesc(covariance_matrix_of_gaussian_process)
colorbar
%%
v = zeros(200,200);
Corre_func = zeros(200,1);
replicates_transpose = replicates';
for i = 1:200
    for j = 1:200-(i-1)
        v(j, i) = dot(replicates_transpose(j,:), replicates_transpose(j+i-1, :));
    end
end

for i = 1:200
    Corre_func (i,1) = sum(v(:,i))/sum(v(:,1));
end
figure(3)
i = 1:200;
plot(i, Corre_func(i), 'LineWidth',2)
hold on
plot(i, exp(-(i-1).^2 / 25),  'LineWidth',2)
xlabel('Distance')
ylabel('Correlation Function')
legend('Correlation Function Calculated from 200 Replicates', 'Gaussian Process Kernel')
set(gca,'XScale','lin', 'YScale','lin','YMinorTick','on','Ydir','normal',...
    'FontSize', 20, 'LineWidth',2) 