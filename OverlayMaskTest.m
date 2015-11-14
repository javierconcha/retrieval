imL8cropmask = imread(...
    '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/LC80160302013262LGN00/MOBELM/LC80160302013262LGN00_ONelm131126testWaterMask2.tif');


%%
% define mask
mask = logical(imL8cropmask);
gray = cat(3,   0.5*ones(size(CDOMmaplog10)), ...
                0.5*ones(size(CDOMmaplog10)),...
                0.5*ones(size(CDOMmaplog10)));
            
figure('name',date)
fs = 16;
ms = 16;
set(gcf,'color','white')
imagesc(gray); % display color of the mask first
hold on
h0 = imagesc(CHLmaplog10); % display Map imaage second
set(h0, 'AlphaData', mask) % Apply transparency to the mask
set(gca,'fontsize',fs)
axis equal
axis image
axis off
h = colorbar('FontSize',11);
set(h,'fontsize',fs,'Location','southoutside')
set(h,'Position',[.2 .1 .6 .05])
title(h,'L8 retrieved C_a [mg m^{-3}]','FontSize',fs)
set(gca, 'Units', 'normalized', 'Position', [0 0.1 1 1])
hold off
y=get(h,'XTick');
set(h,'XTickLabel',10.^y)

