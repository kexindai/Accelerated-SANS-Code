clc
clear
x = [0 1 2];
y = x;
z = [10 12 14];
z_mean = mean(z);
var = std(z).^2;
for i = 1:length(x)
    for j = 1:length(y)
        d(i,j) = sqrt((x(i) - x(j)).^2 + (y(i) - y(j)).^2);
        correlation(i,j) = ((z(i) - z_mean)*(z(j) - z_mean))/var;
    end
end


