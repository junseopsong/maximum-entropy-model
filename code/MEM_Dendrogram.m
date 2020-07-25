clear all; clc;

load 'LocalMinima.txt';
load 'EnergyBarrier.txt';

minimaEnergy = squeeze(LocalMinima(:,4))';
n = size(LocalMinima, 1);
N = 9;

transitionEnergy = zeros(n);
for i = 1:n
    for j = 1:n
        if i ~= j
            transitionEnergy(i,j) = EnergyBarrier(i,j) + max(minimaEnergy(i), minimaEnergy(j));
        end
    end
end

d = nonzeros(tril(transitionEnergy))';

%%
figure(1);
H = dendrogram(linkage(d, 'complete'), 'data', minimaEnergy);
set(gca, 'FontName', 'Arial', 'FontSize', 14);%, 'XColor', 'w');
set(H, 'LineWidth', 2, 'Color', 'k');
ylim([floor(min(minimaEnergy)) ceil(max(d))]);

%%
figure(2);
for j = 1:N
    b(j) = 2^(j-1);
end
for i = 1:n
    minimaState(:,i) = bitand(LocalMinima(i,2), b(:)) ~= 0;
end
dat = double(minimaState);
dat = [dat, dat(:,end)];
dat = [dat; dat(end,:)];
surf(dat);
colormap([0.3 0.3 0.3;1 1 1]); caxis([0 1]);
set(gca, 'Xtick', (1:n)+0.5, 'XtickLabel', 1:n, 'Ytick', (1:N) + 0.5, 'YtickLabel', 1:N);
set(gca, 'FontName', 'Arial', 'FontSize', 14);
axis([1 n+1 1 N+1]); view([0 -90]);

%%
figure(3);
for i = 1:n
    subplot(2, ceil(n/2), i);
    dat = double(reshape(minimaState(:,i),3,3)');
    dat = [dat, dat(:,end)];
    dat = [dat; dat(end,:)];
    surf(dat);
    colormap([0.3 0.3 0.3;1 1 1]); caxis([0 1]);
    axis([1 3+1 1 3+1]); view([0 -90]); axis square;
    axis off; set(gca, 'FontName', 'Arial', 'FontSize', 14);
    title(num2str(i));
end
