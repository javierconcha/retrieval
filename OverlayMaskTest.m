imL8cropmask = imread(...
    '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/LC80160302013262LGN00/MOBELM/LC80160302013262LGN00_ONelm131126testWaterMask2.tif');


%%
figure('name',date)
fs = 16;
ms = 16;
set(gcf,'color','white')
imagesc(CHLmaplog10);
set(gca,'fontsize',fs)
axis equal
axis image
axis off
h = colorbar;
set(h,'fontsize',fs,'Location','southoutside')
set(h,'Position',[.2 .1 .6 .05])
title(h,'L8 retrieved C_a [mg m^{-3}]','FontSize',fs)
set(gca, 'Units', 'normalized', 'Position', [0 0.1 1 1])

mask = logical(imL8cropmask);

gray = cat(3, 0.5*ones(size(CDOMmaplog10)), 0.5*ones(size(CDOMmaplog10)), 0.5*ones(size(CDOMmaplog10)));

hold on
h0 = imshow(gray);
set(h0, 'AlphaData', ~mask)
hold off
%%
[M,N] = size(CDOMmaplog10);
block_size = 50;
P = ceil(M / block_size);
Q = ceil(N / block_size);
alpha_data = checkerboard(block_size, P, Q) > 0;
alpha_data = alpha_data(1:M, 1:N);

mask = logical(imL8cropmask);
set(h0, 'AlphaData', alpha_data);
%%




set(h0, 'AlphaData', ~mask);

%%

Data=magic(100);
c=[1 10 100 1000 10000 10000];
figure
contourf(log(Data(:,:)),log(c));
colormap(bone);  %Color palate named "bone"
caxis(log([c(1) c(length(c))]));
colorbar('FontSize',11,'YTick',log(c),'YTickLabel',c);


