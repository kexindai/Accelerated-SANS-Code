clc
clear
%%
% Sample data
x = [1, 2, 3, 4, 5; 4, 5, 6, 7, 8];
y = [-x + 10];  % Example: y = -x + 10

Corr_coff = corr(x', y');
imagesc(Corr_coff)
colorbar
% Calculate mean and variance of y
%%
mean_y = mean(y);
var_y = var(y)';

% Number of data points
N = length(x);

% Define distance bins
distance_bins = linspace(0, max(x) - min(x), 100);

% Preallocate array for correlation function
correlation_function = zeros(size(distance_bins));

% Calculate correlation for each distance bin
for i = 1:length(distance_bins)
    d = distance_bins(i);
    sum_cov = 0;
    count = 0;
    
    for j = 1:N
        for k = j+1:N
            % Calculate the distance between points
            dist = abs(x(j) - x(k));
            
            % Check if the distance is within the current bin
            if abs(dist - d) < (distance_bins(2) - distance_bins(1)) / 2
                % Calculate the covariance
                sum_cov = sum_cov + (y(j) - mean_y) * (y(k) - mean_y);
                count = count + 1;
            end
        end
    end
    
    % Calculate the average covariance for the current distance
    if count > 0
        correlation_function(i) = sum_cov / (count * var_y);
    else
        correlation_function(i) = NaN;  % No data for this bin
    end
end

% Remove NaNs
valid_idx = ~isnan(correlation_function);
distance_bins = distance_bins(valid_idx);
correlation_function = correlation_function(valid_idx);

% Display the results
disp('Distance bins:');
disp(distance_bins);
disp('Correlation Function:');
disp(correlation_function);

% Plot the correlation function
figure;
plot(distance_bins, correlation_function, 'o-', 'LineWidth', 2);
xlabel('Distance');
ylabel('Correlation');
title('Correlation Function');
grid on;
