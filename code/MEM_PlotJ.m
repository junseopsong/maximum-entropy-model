clear all; clc;

load J.txt;
load h.txt;
load ColorMap_BR.mat;

c = zeros(9,9);
link = [1 2; 1 4; 2 1; 2 3; 2 5; 3 2; 3 6; 4 1; 4 5; 4 7; 5 2; 5 4; 5 6; 5 8; 6 3; 6 5; 6 9; 7 4; 7 8; 8 5; 8 7; 8 9; 9 6; 9 8];
for i = 1:size(link,1)
    c(link(i,1), link(i,2)) = 1;
end

f = figure(1);
imagesc(c); axis square;
set(f, 'ColorMap', [1 1 1; 0 0 0]);
set(gca,'FontName','Arial','FontSize',14);
set(gca,'XTick',1:9,'YTick',1:9);

f = figure(2);
imagesc(J, [-10 10]); cb = colorbar(); axis square;
set(f, 'ColorMap', cmap); %set(cb,'YTick',-8:8);
set(gca,'FontName','Arial','FontSize',14);
set(gca,'XTick',1:9,'YTick',1:9);

f = figure(3);
imagesc(h', [-12 12]); cb = colorbar();
set(f, 'ColorMap', cmap); %set(cb,'YTick',-3:3);
set(gca,'FontName','Arial','FontSize',14);
set(gca,'XTick',[],'YTick',[]);
set(gcf,'Units','normal')
set(gca,'Position',[0.1 0.1 0.7 0.7])
