function [] = phaseplot(thetas)
[l, num_rhythms] = size(thetas);
figure
for j = 1:num_rhythms
    subplot(num_rhythms, 1, j);
    %plot(1:1:1000, thetas(1:1000, j));
    plot(1:1:l, thetas(:, j));
end