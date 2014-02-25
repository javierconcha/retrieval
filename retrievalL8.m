% Retrieve the concentration from a image using a LUT from HL
cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/retrieval/
%% L5 image
% % % % im0115big = imread('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L5images/LT50160302000115AAA02/LT50160302000115AAA02radtif.tif');
% % % % im0115bigmask = imread('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L5images/LT50160302000115AAA02/LT50160302000115AAA02watermask.tif');
% % % % im0115bigmask(im0115bigmask>0)=1;
% % % % 
% % % % L5bands = [0.485,0.560,0.660,0.830,1.650,2.220];
% % % % %% for displaying
% % % % im0115R = im0115big(:,:,3);
% % % % im0115G = im0115big(:,:,2);
% % % % im0115B = im0115big(:,:,1);
% % % % 
% % % % im0115R = im0115R/max(im0115R(:));% normalize for maximum
% % % % im0115G = im0115G/max(im0115G(:));
% % % % im0115B = im0115B/max(im0115B(:));
% % % % 
% % % % im0115RGB(:,:,1) = histeq(im0115R);% normalize for maximum
% % % % im0115RGB(:,:,2) = histeq(im0115G);
% % % % im0115RGB(:,:,3) = histeq(im0115B);
% % % % 
% % % % figure
% % % % fs = 15;
% % % % set(gcf,'color','white')
% % % % imagesc(double(im0115RGB))
% % % % axis equal
% % % % hold on
% % % % 
% % % % % plot(2778,2123,'g.')
% % % % % plot(2814,2123,'g.')
% % % % % plot(2814,2159,'g.')
% % % % % plot(2778,2159,'g.')
% % % % 
% % % % xA = 2500;
% % % % xB = 2600;
% % % % yA = 2000;
% % % % yB = 2100;
% % % % 
% % % % 
% % % % plot(xA,yA,'g.')
% % % % plot(xB,yA,'g.')
% % % % plot(xB,yB,'g.')
% % % % plot(xA,yB,'g.')
% % % % %% water mask for big image
% % % % figure
% % % % imagesc(im0115bigmask)
% % % % colormap(gray)
% % % % axis equal
% % % % 
% % % % %% ROI middle of the lake
% % % % % RR = im0115big(2778:2814,2123:2159,:);
% % % % RR = im0115big(xA:xB,yA:yB,:);
% % % % ROImiddleOfTheLake = reshape(RR,[size(RR,1)*size(RR,2) size(RR,3)]);
% % % % 
% % % % 
% % % % figure
% % % % fs = 15;
% % % % set(gcf,'color','white')
% % % % plot(L5bands,ROImiddleOfTheLake')
% % % % title('ROI middle of the lake','fontsize',fs)
% % % % xlabel('wavelength [ \mum]','fontsize',fs)
% % % % ylabel('Radiance [W/m^2/sr/um]','fontsize',fs)
% % % % set(gca,'fontsize',fs)
% % % % %% Deglint
% % % % figure
% % % % fs = 15;
% % % % set(gcf,'color','white')
% % % % plot(ROImiddleOfTheLake(:,4),ROImiddleOfTheLake(:,1),'.')
% % % % title('ROI middle of the lake','fontsize',fs)
% % % % xlabel('NIR','fontsize',fs)
% % % % ylabel('Band 1','fontsize',fs)
% % % % set(gca,'fontsize',fs)
% % % % axis equal
% % % % 
% % % % 
% % % % 
% % % % figure
% % % % fs = 15;
% % % % set(gcf,'color','white')
% % % % plot(ROImiddleOfTheLake(:,4),ROImiddleOfTheLake(:,2),'.')
% % % % title('ROI middle of the lake','fontsize',fs)
% % % % xlabel('NIR Band','fontsize',fs)
% % % % ylabel('Band 2','fontsize',fs)
% % % % set(gca,'fontsize',fs)
% % % % axis equal
% % % % 
% % % % figure
% % % % fs = 15;
% % % % set(gcf,'color','white')
% % % % plot(ROImiddleOfTheLake(:,4),ROImiddleOfTheLake(:,3),'.')
% % % % title('ROI middle of the lake','fontsize',fs)
% % % % xlabel('NIR','fontsize',fs)
% % % % ylabel('Band 3','fontsize',fs)
% % % % set(gca,'fontsize',fs)
% % % % axis equal
% % % % 
% % % % figure
% % % % fs = 15;
% % % % set(gcf,'color','white')
% % % % plot(ROImiddleOfTheLake(:,4),ROImiddleOfTheLake(:,5),'.')
% % % % title('ROI middle of the lake','fontsize',fs)
% % % % xlabel('NIR','fontsize',fs)
% % % % ylabel('Band 5','fontsize',fs)
% % % % set(gca,'fontsize',fs)
% % % % axis equal
% % % % 
% % % % figure
% % % % fs = 15;
% % % % set(gcf,'color','white')
% % % % plot(ROImiddleOfTheLake(:,4),ROImiddleOfTheLake(:,6),'.')
% % % % title('ROI middle of the lake','fontsize',fs)
% % % % xlabel('NIR','fontsize',fs)
% % % % ylabel('Band 7','fontsize',fs)
% % % % set(gca,'fontsize',fs)
% % % % axis equal
% % % % 
% % % % 
% % % % [bi1,minNIR] = glint_removal(ROImiddleOfTheLake(:,1),ROImiddleOfTheLake(:,4));
% % % % [bi2,~] = glint_removal(ROImiddleOfTheLake(:,2),ROImiddleOfTheLake(:,4));
% % % % [bi3,~] = glint_removal(ROImiddleOfTheLake(:,3),ROImiddleOfTheLake(:,4));
% % % % [bi5,~] = glint_removal(ROImiddleOfTheLake(:,5),ROImiddleOfTheLake(:,4));
% % % % [bi7,~] = glint_removal(ROImiddleOfTheLake(:,6),ROImiddleOfTheLake(:,4));
% % % % 
% % % % bi1
% % % % bi2
% % % % bi3
% % % % bi5
% % % % bi7
% % % % minNIR
% % % % 
% % % % %% L5 water pixels 
% % % % imnew = reshape(im0115big,[size(im0115big,1)*size(im0115big,2) size(im0115big,3)]);
% % % % masknew = reshape(im0115bigmask,[size(im0115bigmask,1)*size(im0115bigmask,2) size(im0115bigmask,3)]);
% % % % 
% % % % waterpixelsL5 = imnew(masknew==1,:);
% % % % waterpixelsL5 = double(waterpixelsL5);
% % % % 
% % % % % %% display All water pixels
% % % % % 
% % % % % figure
% % % % % fs = 15;
% % % % % set(gcf,'color','white')
% % % % % plot(L5bands,waterpixelsL5(1:1000000,:))
% % % % % title('Reflectance water L5 image','fontsize',fs)
% % % % % xlabel('wavelength [\mu m]','fontsize',fs)
% % % % % ylabel('radiance [W/m^2/sr/um]','fontsize',fs)
% % % % % set(gca,'fontsize',fs)
% % % % 
% % % % 
% % % % %%
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % %% read tif
% % % % % filename = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L5images/LT50160302000115AAA02/LT50160302000115AAA02_ON130408resampledtif.twf';
% % % % % info = imfinfo(filename);
% % % % % t = Tiff(filename, 'r'); 
% % % % % % t.getTag('GeoKeyDirectoryTag')
% % % % % %%
% % % % % % Get the dimensions of the image to calculate coordinates.
% % % % % numRows = t.getTag('ImageLength');
% % % % % numCols = t.getTag('ImageWidth');
% % % % % % Get the ID number of the tile containing the coordinates.
% % % % % tileNum = t.computeStrip([numRows numCols],1);
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% L8 image cropped

% imL8crop = imread(...
%     '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/LC80170302013237LGN00/LC80170302013237LGN00_ONelm130917resampledtif.tif');

% Bright=city pixel;Dark=Reflectance Measurement
% imL8crop = imread(...
%     '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/LC80160302013262LGN00/LC80160302013262LGN00_ONelm131209test.tif');
  

% 01/12/14 ELM with two pixels with known reflectances - traditional ELM
% imL8crop = imread(...
%     '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/LC80160302013262LGN00/LC80160302013262LGN00_ONelm140112testtif.tif');

% 01/21/14 ELM with two pixels with known reflectances
imL8crop = imread(...
    '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/LC80160302013262LGN00/LC80160302013262LGN00_ONelm140121testtif.tif');

%%%% Mask
% imL8cropmask = imread(...
%     '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/LC80170302013237LGN00/LC80170302013237LGN00_ONmaskmin0p5resampled.tif');
imL8cropmask = imread(...
    '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/LC80160302013262LGN00/LC80160302013262LGN00_ONelm131126testWaterMask2.tif');

imL8cropmask(imL8cropmask>0)=1;


L8bands = [0.4430,0.4826,0.5613,0.6546,0.8646,1.6090,2.2010];

%%%% water pixels. Convert each band in columns.

imnew = reshape(imL8crop,[size(imL8crop,1)*size(imL8crop,2) size(imL8crop,3)]);
masknew = reshape(imL8cropmask,[size(imL8cropmask,1)*size(imL8cropmask,2) size(imL8cropmask,3)]);



waterpixels = imnew(masknew==1,:);
waterpixels = double(waterpixels);

% added 01-11-14. plot radiance curves
% radiance image
imL8radcrop = imread(...
    '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L8images/LC80160302013262LGN00/LC80160302013262LGN00rad_ONelm.tif');

imradnew = reshape(imL8radcrop,[size(imL8radcrop,1)*size(imL8radcrop,2) size(imL8radcrop,3)]);

waterradpixels = imradnew(masknew==1,:);
waterradpixels = double(waterradpixels);
%% for displaying
imL8cropRGB(:,:,1)=imadjust(imL8crop(:,:,4));
imL8cropRGB(:,:,2)=imadjust(imL8crop(:,:,3));
imL8cropRGB(:,:,3)=imadjust(imL8crop(:,:,2));


impos = double(imL8cropRGB);
% impos(impos<0)=0;% only positive values

maskRGB(:,:,1)=double(imL8cropmask);
maskRGB(:,:,2)=double(imL8cropmask);
maskRGB(:,:,3)=double(imL8cropmask);
 

impos = impos.*maskRGB;


figure
set(gcf,'color','white')
imagesc(impos)
axis equal

%% mask display
figure
set(gcf,'color','white')
imshow(imadjust(imL8cropmask))



%% Stats water pixels
meanwp = mean(waterpixels,1);
stdwp = std(waterpixels,1);
maxwp = max(waterpixels,[],1);
minwp = min(waterpixels,[],1);


figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,meanwp,'k')
title('Reflectance water L8 image','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
hold on
plot(L8bands,meanwp+stdwp,'g')
plot(L8bands,meanwp-stdwp,'g')
plot(L8bands,maxwp,'r')
plot(L8bands,minwp,'r')
xlim([min(L8bands) max(L8bands)])

format short
disp('---------------------------------------------------')
disp('Basic Stats      Min       Max       Mean     Stdev  ') 
disp('---------------------------------------------------')

bands = [1 2 3 4 5 6 7];

for i = 1:size(minwp,2)
    str = sprintf('     band %i  %2.6f  %2.6f  %2.6f  %2.6f  %2.6f  %2.6f',bands(i),minwp(i), maxwp(i), meanwp(i), stdwp(i));
    disp(str)
end

nbins = 100;

figure
subplot(2,4,1)
fs = 15;
set(gcf,'color','white')
hist(waterpixels(:,1),nbins)
title('band 1','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,4,2)
fs = 15;
set(gcf,'color','white')
hist(waterpixels(:,2),nbins)
title('band 2','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,4,3)
fs = 15;
set(gcf,'color','white')
hist(waterpixels(:,3),nbins)
title('band 3','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,4,4)
fs = 15;
set(gcf,'color','white')
hist(waterpixels(:,4),nbins)
title('band 4','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,4,5)
fs = 15;
set(gcf,'color','white')
hist(waterpixels(:,5),nbins)
title('band 5','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,4,6)
fs = 15;
set(gcf,'color','white')
hist(waterpixels(:,6),nbins)
title('band 6','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,4,7)
fs = 15;
set(gcf,'color','white')
hist(waterpixels(:,7),nbins)
title('band 7','fontsize',fs)
set(gca,'fontsize',fs)

p=mtit('Histogram Water Pixels per Band',...
 	     'fontsize',fs+1,'xoff',0,'yoff',.025);
     
%% display All water pixels

figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,waterpixels')
title('Reflectance water L8 image','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
% ylim([0 0.18])
 hold on
 plot(L8bands,meanwp,'g','linewidth',2)
 
%% display All water pixels RADIANCE values

figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,waterradpixels')
title('Radiance water L8 image','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('Radiance [W/m^2/sr]','fontsize',fs)
set(gca,'fontsize',fs)
% ylim([0 0.18])
%  hold on
%  plot(L8bands,meanwp+stdwp,'g','linewidth',2) 

%% negative values
im = double(imL8crop);
imneg = zeros(size(im));
imneg(im<0)=im(im<0);% only negatives

imnegmask = zeros(size(imL8crop));% for displaying
imnegmask(im<0)=1; % negative values are white
imnegmask(im>=0)=0; % positive values are black
imnegmask = imnegmask+0.5*repmat(~double(imL8cropmask),[1 1 size(imL8crop,3)]); % for the land appear gray

% figure
% subplot(3,3,1)
% fs = 15;
% set(gcf,'color','white')
% imshow(imnegmask(:,:,1))
% title('band 1','fontsize',fs)
% set(gca,'fontsize',fs)
% 
% subplot(3,3,2)
% fs = 15;
% set(gcf,'color','white')
% imshow(imnegmask(:,:,2))
% title('band 2','fontsize',fs)
% set(gca,'fontsize',fs)
% 
% subplot(3,3,3)
% fs = 15;
% set(gcf,'color','white')
% imshow(imnegmask(:,:,3))
% title('band 3','fontsize',fs)
% set(gca,'fontsize',fs)
% 
% subplot(3,3,4)
% fs = 15;
% set(gcf,'color','white')
% imshow(imnegmask(:,:,4))
% title('band 4','fontsize',fs)
% set(gca,'fontsize',fs)

figure
subplot(2,2,1)
fs = 15;
set(gcf,'color','white')
imshow(imnegmask(:,:,5))
title('band 5','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,2,2)
fs = 15;
set(gcf,'color','white')
imshow(imnegmask(:,:,6))
title('band 6','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,2,3)
fs = 15;
set(gcf,'color','white')
imshow(imnegmask(:,:,7))
title('band 7','fontsize',fs)
set(gca,'fontsize',fs)


p=mtit('Negative Values',...
 	     'fontsize',fs+1,'xoff',0,'yoff',.025);
     
     
format short
disp('----------------------------------------------------------')
disp('Basic Stats      Min       Max       Mean     Stdev    N  ') 
disp('----------------------------------------------------------')

bands = [5 6 7];

for i = 1:size(bands,2)
    data = imneg(:,:,bands(i));
    str = sprintf('     band %i  %2.6f  %2.6f  %2.6f  %2.6f  %2.6f  %2.6f %6.0i'...
        ,bands(i),min(data(:)), max(data(:)), mean(data(:)), std(data(:)),sum(sum(data<0)));
    disp(str)
end
     
%% Display negative values
bn = 7;
waterpixels_neg = waterpixels(waterpixels(:,bn)<0,:);

figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,waterpixels_neg')
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
str = sprintf('Curves with negatives values in Band %i',bn);
title(str,'fontsize',fs)
% ylim([0 0.18])

%% Display high values in the different bands
im = double(imL8crop);
bn = 5;
imhighNIR = zeros(size(im,1),size(im,2));
cond2 = im(:,:,bn)> (meanwp(bn)+2*stdwp(bn)) & imL8cropmask ~= 0;
imhighNIR(cond2)=1;% only high NIR

imhighNIRmask = imhighNIR+0.5*repmat(~double(imL8cropmask),[1 1 1]); % for the land appear gray


figure
fs = 15;
set(gcf,'color','white')
imshow(imhighNIRmask)
str = sprintf('High values for band %i',bn);
title(str,'fontsize',fs)
set(gca,'fontsize',fs)

% display All water pixels with no high values in NIR

figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,waterpixels(waterpixels(:,bn)<meanwp(bn)+2*stdwp(bn),:)')
str = sprintf('Reflectance water L8 image with low values band %i',bn);
title(str,'fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
% ylim([0 0.18])
hold on
plot(L8bands,meanwp+stdwp,'g','linewidth',2)

figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,waterpixels(waterpixels(:,bn)>meanwp(bn)+2*stdwp(bn),:)')
str = sprintf('Reflectance water L8 image with high values band %i',bn);
title(str,'fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
% ylim([0 0.18])
hold on
plot(L8bands,meanwp+stdwp,'g','linewidth',2)

%% Pixels that B5 > B3
im = double(imL8crop);
imB5greaterthanB3 = zeros(size(im,1),size(im,2));
cond3 = (im(:,:,5)> im(:,:,3) )& (imL8cropmask ~= 0);
imB5greaterthanB3(cond3)=1;% only high NIR

imB5greaterthanB3mask = imB5greaterthanB3+0.5*repmat(~double(imL8cropmask),[1 1 1]); % for the land appear gray


figure
fs = 15;
set(gcf,'color','white')
imshow(imB5greaterthanB3mask)
str = sprintf('High values for band %i',bn);
title(str,'fontsize',fs)
set(gca,'fontsize',fs)

figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,waterpixels(waterpixels(:,5)>waterpixels(:,3),:)')
str = sprintf('Reflectance water L8 image with high values band %i',bn);
title(str,'fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
% ylim([0 0.18])
hold on
plot(L8bands,meanwp+stdwp,'g','linewidth',2)

%% Pixels para incluir en el IGARSS14 abstract
% water pixels reflectance
waterpixelsamples = waterpixels(1:70:end,:);
figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,waterpixelsamples(waterpixelsamples(:,5)<waterpixelsamples(:,3),:)')
str = sprintf('Reflectance water pixels L8');
title(str,'fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
% ylim([0 .25])

% water pixels Radiance
% figure
% fs = 15;
% set(gcf,'color','white')
% plot(L8bands,waterradpixels(waterpixels(:,5)<waterpixels(:,3),:)')
% str = sprintf('Radiance water pixels L8');
% title(str,'fontsize',fs)
% xlabel('wavelength [\mu m]','fontsize',fs)
% ylabel('Radiance [W/m^2/sr]','fontsize',fs)
% set(gca,'fontsize',fs)



%% LUTs from HydroLight
% LUT1temp = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/ReformedLUT.txt');
% LUT2temp = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/ReformedLUT1.txt');
% % LUTconc = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/ReformedLUTConc.txt');
% % LUTconc = LUTconc';%ex: 1000x3
% 
% 
% wl120 = LUT1temp(:,1);
% wl140 = LUT2temp(:,1);
% 
% clear LUT1temp LUT2temp

% LUT1 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/ReformedLUT1speclyb.txt');
% LUT2 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/ReformedLUT2speclyb.txt');
% LUT1 = LUT1';% ex: 1000x6
% LUT2 = LUT2';% ex: 1000x6

% % % before 04/21/13
% % rr = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/LUTjavierL8.txt');
% % LUT3 = rr(:,2:end)';
% % LUTconc3 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/concentration_list.txt');
% % 
% % % LUT1000 = LUT3(LUTconc(:,2)<=24,:);
% % % LUT1000conc = LUTconc(LUTconc(:,2)<=24,:);
% % %
% % 
% % 
% % 
% % % new 04/21/13
% % rr = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/LUTL8130421.txt');
% % LUT1200 = rr(:,2:end)';
% % LUT1200conc = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/concentration_list130421.txt');
% % 
% % % new 05/08/13
% % rr = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/LUTL8130508.txt');
% % LUT1568 = rr(:,2:end)';
% % LUT1568conc = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/concentration_list130508.txt');


% new 05/08/13
% rr = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/LUTL5130511.txt');
% LUT2240 = rr(:,2:end)';
% LUT2240conc = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/concentration_list130511.txt');

% % new 10/21/13
% rr = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/LUTL8131021.txt');
% LUT2240 = rr(:,2:end)';
% LUT2240conc = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/concentration_list130511.txt');

% new 12/18/13
% rr = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/LUTL8130919.txt'); % Created for 09/19/13 image!!!
% LUT = rr(:,2:end)';
% LUTconc = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/concentration_list130919_2.txt');

% new 02/06/14
rr = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/LUTL8130919_3.txt'); % Created for 09/19/13 image!!!
LUT = rr(:,2:end)';
LUTconc = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/concentration_list130919_2.txt');



% LUT1000 = LUT1200(LUT1200conc(:,1)<=100,:);
% LUT1000conc = LUT1200conc(LUT1200conc(:,1)<=100,:);

% cond = LUT1200conc(:,1)<=100;
% LUT1000 = LUT1200(cond,:);
% LUT1000conc = LUT1200conc(cond,:);


%% Display LUT
figure
fs = 15;
set(gcf,'color','white')
plot(L8bands,LUT)
title('Reflectance LUT from HydroLight','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
%  hold on
%  plot(L8bands,meanwp+stdwp,'g','linewidth',2)
% % 
% % figure
% % fs = 15;
% % set(gcf,'color','white')
% % plot(L8bands,LUT2)
% % title('Reflectance LUT2 - HydroLight','fontsize',fs)
% % xlabel('wavelength [\mu m]','fontsize',fs)
% % ylabel('reflectance','fontsize',fs)
% % set(gca,'fontsize',fs)
%% Darkest 
% % Darkest in red band
% [~,darkestindex] = min(waterpixels(:,3));
% waterdarkestred = waterpixels(darkestindex,:);
% 
% [~,LUTdarkestindex] = min(LUT1000(:,3));
% LUTdarkestred = LUT1000(LUTdarkestindex,:);
% LUT1000conc(LUTdarkestindex,:)
% 
% % Darkest in Magnitude
% wpnorm = sqrt(sum(waterpixels.^2,2)); % for the waterpixels
% [normmin,index] = min(wpnorm);
% waterpixelspermag = waterpixels(index,:);
% 
% LUT1000norm = sqrt(sum(LUT1000.^2,2)); % for the LUT
% [normmin,index] = min(LUT1000norm);
% LUT1000conc(index,:)
% 
% figure
% fs = 15;
% set(gcf,'color','white')
% plot(L8bands,waterdarkestred')
% hold on
% plot(L8bands,waterpixelspermag','c')
% plot(L8bands,LUTdarkestred','r')
% plot(L8bands,LUT1000(1,:)','k')
% 
% legend('L8 per b3','L8 per mag','HL:0.25;0.25;14','HL:0.25;0.25;0.25')
% title('Darkest pixel - L8 vs HL','fontsize',fs)
% xlabel('wavelength [\mu m]','fontsize',fs)
% ylabel('reflectance','fontsize',fs)
% set(gca,'fontsize',fs)
% % ylim([0 0.2])

%% Test the optimization algorhythm
LUTconc = LUTconc;
CDOMconc = unique(LUTconc(:,3))
SMconc   = unique(LUTconc(:,2))
CHLconc  = unique(LUTconc(:,1))

disp('--------------------------------------------------------------------------')
disp('Running Optimization Routine')
    [XResultstest,residual] = opt(LUT(:,1:5),LUT(:,1:5),LUTconc);
disp('Optimization Routine finished Successfully')

% E_RMS
disp('--------------------------------------------------')
E_Chl = sqrt(sum((XResultstest(:,1)-LUTconc(:,1)).^2)/size(XResultstest,1));
E_Chl = E_Chl*100/68;
str = sprintf('E_Chl  = %2.2f %%',E_Chl);
disp(str)

E_SM = sqrt(sum((XResultstest(:,2)-LUTconc(:,2)).^2)/size(XResultstest,1));
E_SM = E_SM*100/24;
str = sprintf('E_SM   = %2.2f %%',E_SM);
disp(str)

E_CDOM = sqrt(sum((XResultstest(:,3)-LUTconc(:,3)).^2)/size(XResultstest,1));
E_CDOM = E_CDOM*100/14;
str = sprintf('E_CDOM = %2.2f %%',E_CDOM);
disp(str)
%% Residual Histogram
figure
set(gcf,'color','white')
subplot(2,3,1)
hist(residual(:,1))
title('band 1')

subplot(2,3,2)
hist(residual(:,2))
title('band 2')

subplot(2,3,3)
hist(residual(:,3))
title('band 3')

subplot(2,3,4)
hist(residual(:,4))
title('band 4')

subplot(2,3,5)
hist(residual(:,5))
title('band 5')
%% Display data vs retrieved

figure
fs = 15;
set(gcf,'color','white')
plot(LUTconc(:,1),XResultstest(:,1),'.')
xLimits = get(gca,'XLim');  %# Get the range of the x axis
yLimits = get(gca,'YLim');  %# Get the range of the y axis
hold on
plot(xLimits,xLimits,'k')
ylim(xLimits)
xlim(xLimits)
title('CHL Real vs retrieved','fontsize',fs)
xlabel('real','fontsize',fs)
ylabel('retrieved','fontsize',fs)
set(gca,'fontsize',fs)
axis equal

figure
fs = 15;
set(gcf,'color','white')
plot(LUTconc(:,2),XResultstest(:,2),'.')
xLimits = get(gca,'XLim');  %# Get the range of the x axis
yLimits = get(gca,'YLim');  %# Get the range of the y axis
hold on
plot(xLimits,xLimits,'k')
ylim(xLimits)
xlim(xLimits)
title('SM Real vs retrieved','fontsize',fs)
xlabel('real','fontsize',fs)
ylabel('retrieved','fontsize',fs)
set(gca,'fontsize',fs)
axis equal

figure
fs = 15;
set(gcf,'color','white')
plot(LUTconc(:,3),XResultstest(:,3),'.')
xLimits = get(gca,'XLim');  %# Get the range of the x axis
yLimits = get(gca,'YLim');  %# Get the range of the y axis
hold on
plot(xLimits,xLimits,'k')
ylim(xLimits)
xlim(xLimits)
title('CDOM Real vs retrieved','fontsize',fs)
xlabel('real','fontsize',fs)
ylabel('retrieved','fontsize',fs)
set(gca,'fontsize',fs)
axis equal

%% Retrieval %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('--------------------------------------------------------------------------')
disp('Running Optimization Routine')
    XResults = opt(waterpixels(:,1:5),LUT(:,1:5),LUTconc);
disp('Optimization Routine finished Successfully')
%% Mapping Concentrations log scale

landmask = imL8cropmask==1;

ConcRet = zeros(size(masknew),3);
ConcRet(masknew==1,:) = XResults; % Concentration Retrieved

CHLmap  = reshape(ConcRet(:,1),[size(imL8cropmask,1) size(imL8cropmask,2) size(imL8cropmask,3)]);
SMmap   = reshape(ConcRet(:,2),[size(imL8cropmask,1) size(imL8cropmask,2) size(imL8cropmask,3)]);
CDOMmap = reshape(ConcRet(:,3),[size(imL8cropmask,1) size(imL8cropmask,2) size(imL8cropmask,3)]);

CHLmaplog10 = log10(CHLmap);
CHLmaplog10(CHLmaplog10==-Inf)=-12;
% CHLmaplog10masked = bsxfun(@times, CHLmaplog10, landmask);

SMmaplog10 = log10(SMmap);
SMmaplog10(SMmaplog10==-Inf)=-12;
% SMmaplog10masked = bsxfun(@times, SMmaplog10, landmask);

CDOMmaplog10 = log10(CDOMmap);
CDOMmaplog10(CDOMmaplog10==-Inf)=-12;
% CDOMmaplog10masked = bsxfun(@times, CDOMmaplog10, landmask);

figure
fs = 15;
set(gcf,'color','white')
imagesc(CHLmaplog10)
title('CHL map log scale','fontsize',fs)
set(gca,'fontsize',fs)
axis('equal')
colorbar

figure
fs = 15;
set(gcf,'color','white')
imagesc(SMmaplog10)
title('SM map log scale','fontsize',fs)
set(gca,'fontsize',fs)
axis('equal')
colorbar

figure
fs = 15;
set(gcf,'color','white')
imagesc(CDOMmaplog10)
title('CDOM map log scale','fontsize',fs)
set(gca,'fontsize',fs)
axis('equal')
colorbar
%% Mapping Concentrations linear scale
figure
fs = 15;
set(gcf,'color','white')
imagesc(CHLmap)
title('CHL map linear scale','fontsize',fs)
set(gca,'fontsize',fs)
axis('equal')
colorbar

figure
fs = 15;
set(gcf,'color','white')
imagesc(SMmap)
title('SM map linear scale','fontsize',fs)
set(gca,'fontsize',fs)
axis('equal')
colorbar

figure
fs = 15;
set(gcf,'color','white')
imagesc(CDOMmap)
title('CDOM map linear scale','fontsize',fs)
set(gca,'fontsize',fs)
axis('equal')
colorbar
%% Histogram of concentrations log scale
nbins = 50;
figure
subplot(1,3,1)
fs = 15;
set(gcf,'color','white')
hist(CHLmaplog10(CHLmaplog10~=-Inf),nbins)
title('CHL','fontsize',fs)
set(gca,'fontsize',fs)

subplot(1,3,2)
fs = 15;
set(gcf,'color','white')
hist(SMmaplog10(SMmaplog10~=-Inf),nbins)
title('SM','fontsize',fs)
set(gca,'fontsize',fs)

subplot(1,3,3)
fs = 15;
set(gcf,'color','white')
hist(CDOMmaplog10(CDOMmaplog10~=-Inf),nbins)
title('CDOM','fontsize',fs)
set(gca,'fontsize',fs)
%% Histogram of concentrations linear scale
nbins = 50;
figure
subplot(1,3,1)
fs = 15;
set(gcf,'color','white')
hist(XResults(:,1),nbins)
title('CHL','fontsize',fs)
set(gca,'fontsize',fs)

subplot(1,3,2)
fs = 15;
set(gcf,'color','white')
hist(XResults(:,2),nbins)
title('SM','fontsize',fs)
set(gca,'fontsize',fs)

subplot(1,3,3)
fs = 15;
set(gcf,'color','white')
hist(XResults(:,3),nbins)
title('CDOM','fontsize',fs)
set(gca,'fontsize',fs)

%% Comparison between two LUTs from Aaron and one curve from HL5.1 in tropos
% % % 
% % % r_02502514 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/r_025025014.txt');
% % % figure
% % % plot(LUT1(10,:))
% % % hold on
% % % plot(LUT1000(10,:),'k')
% % % % plot(r_02502514(:,2),'r')

%% Figure for Dr. John - 09/23/13 Linear Scale
figure
set(gcf,'color','white')
subplot(2,2,1)
imagesc(impos)
title('RGB image ','fontsize',fs)
set(gca,'fontsize',fs)
axis equal
axis off

subplot(2,2,2)
fs = 15;
% set(gcf,'color','white')
imagesc(CHLmap)
title('CHL map ','fontsize',fs)
set(gca,'fontsize',fs)
colorbar
axis equal
axis off

subplot(2,2,3)
fs = 15;
% set(gcf,'color','white')
imagesc(SMmap)
title('SM map ','fontsize',fs)
set(gca,'fontsize',fs)
colorbar
axis equal
axis off

subplot(2,2,4)
fs = 15;
% set(gcf,'color','white')
imagesc(CDOMmap)
title('CDOM map ','fontsize',fs)
set(gca,'fontsize',fs)
colorbar
axis equal
axis off

%% Figure for Dr. John - 09/23/13 Log Scale
figure
fs = 20;
set(gcf,'color','white')
subplot(2,2,1)
imagesc(impos)
title('RGB image ','fontsize',fs)
set(gca,'fontsize',fs)
axis equal
axis off

subplot(2,2,2)
imagesc(CHLmaplog10)
title('CHL map ','fontsize',fs)
set(gca,'fontsize',fs)
colorbar
axis equal
axis off

subplot(2,2,3)
imagesc(SMmaplog10)
title('SM map ','fontsize',fs)
set(gca,'fontsize',fs)
colorbar
axis equal
axis off

subplot(2,2,4)
imagesc(CDOMmaplog10)
title('CDOM map ','fontsize',fs)
set(gca,'fontsize',fs)
colorbar
axis equal
axis off

%% Comparison one curve with new IOP

j = (LUTconc(:,1) == 2.5 & LUTconc(:,2) == 0 & LUTconc(:,3) == 0);
%0.300000000000000   0.250000000000000   0.135000000000000
figure
fs = 15;
set(gcf,'color','white')
plot(L8bands*1000,LUT(j,:))
hold on
plot(wavelength,newR025,'r')
title('Reflectance - HydroLight','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
legend('Old: 2.5 0.0 0.0 ','New 0.25 0.25 0.25')
set(gca,'fontsize',fs)
%% Figure for IGARSS14 abstract - Cranberry
ROIcranELMbased = [
0.016134;
0.020274;
0.045366;
0.027074;
0.008125];

ROIcranELM = [
0.017207;
0.022009;
0.053848;
0.032086;
0.008883];

figure
ylimit = [0 0.06];
fs = 15;
set(gcf,'color','white')
plot(L8bands(1:5)*1000,ROIcranELMbased)
hold on
plot(L8bands(1:5)*1000,ROIcranELM,'k')
title('Cranberry Pond ROI','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
legend('ELM-based ','ELM')
set(gca,'fontsize',fs)
ylim(ylimit)

%% Figure for IGARSS14 abstract - Long Pond
ROIlonELMbased = [
0.014342;
0.017705;
0.032769;
0.021847;
0.009542];

ROIlonELM = [
0.014612;
0.018567;
0.038103;
0.025818;
0.010349];

figure
ylimit = [0 0.06];
fs = 15;
set(gcf,'color','white')
plot(L8bands(1:5)*1000,ROIlonELMbased)
hold on
plot(L8bands(1:5)*1000,ROIlonELM,'k')
title('Long Pond ROI','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
legend('ELM-based ','ELM')
set(gca,'fontsize',fs)
ylim(ylimit)

%% Figure for IGARSS14 abstract - Near Shore Lake
ROInearELMbased = [
0.017452;
0.022167;
0.021507;
0.005976;
0.001954];

ROInearELM = [
0.019115;
0.024545;
0.024026;
0.006784;
0.0025];

figure
ylimit = [0 0.06];
fs = 15;
set(gcf,'color','white')
plot(L8bands(1:5)*1000,ROInearELMbased)
hold on
plot(L8bands(1:5)*1000,ROInearELM,'k')
title('Near shore lake','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
legend('ELM-based ','ELM')
set(gca,'fontsize',fs)
ylim(ylimit)

%% Figure for IGARSS14 abstract - Off Shore Lake
ROIoffELMbased = [
0.027888;
0.035044;
0.028504;
0.007039;
0.00107];

ROIoffELM = [
0.034223;
0.041796;
0.032772;
0.008058;
0.001585];

figure
ylimit = [0 0.06];
fs = 15;
set(gcf,'color','white')
plot(L8bands(1:5)*1000,ROIoffELMbased)
hold on
plot(L8bands(1:5)*1000,ROIoffELM,'k')
title('Off shore lake','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
legend('ELM-based ','ELM')
set(gca,'fontsize',fs)
ylim(ylimit)

%% Figure for IGARSS14 abstract - All

figure
ylimit = [0 0.06];
fs = 15;
set(gcf,'color','white')
plot(L8bands(1:5)*1000,ROIcranELMbased,'-bx')
hold on
plot(L8bands(1:5)*1000,ROIcranELM,'--bx')
plot(L8bands(1:5)*1000,ROIlonELMbased,'-ro')
plot(L8bands(1:5)*1000,ROIlonELM,'--ro')
plot(L8bands(1:5)*1000,ROInearELMbased,'-mv')
plot(L8bands(1:5)*1000,ROInearELM,'--mv')
plot(L8bands(1:5)*1000,ROIoffELMbased,'-k+')
plot(L8bands(1:5)*1000,ROIoffELM,'--k+')


% title('Off shore lake','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
legend('Cranberry ','Cranberry','Long','Long','Nearshore',...
    'Nearshore','Offshore','Offshore')
set(gca,'fontsize',fs)
ylim(ylimit)

%% Stat for IGARSS14 abstract
disp('reflectance units')
diffELM = [ROIcranELMbased-ROIcranELM ROIlonELMbased-ROIlonELM ...
    ROInearELMbased-ROInearELM ROIoffELMbased-ROIoffELM];
diffELM = 100*abs(diffELM); % in reflectance units
max(max(diffELM,[],2))
min(min(diffELM,[],2))  

disp('percentage difference')
diffELM = [(ROIcranELMbased-ROIcranELM)./max([ROIcranELM ROIcranELMbased],[],2) ...
    (ROIlonELMbased-ROIlonELM)./max([ROIlonELM ROIlonELMbased],[],2) ...
    (ROInearELMbased-ROInearELM)./max([ROInearELM ROInearELMbased],[],2) ...
    (ROIoffELMbased-ROIoffELM)./max([ROIoffELM ROIoffELMbased],[],2)];
diffELM = 100*abs(diffELM); % in reflectance units
max(max(diffELM,[],2))
min(min(diffELM,[],2)) 
%% SNR Landsat-8
SNRL8Ltyp = [
130;
130;
100;
90;
90;
100;
100];

SNRL7Ltyp = [
72;
49;
30;
27;
19;
14];

SNRL8 = [
439;
402;
256;
207;
156;
134;
111];

figure
fs = 15;
set(gcf,'color','white')
plot(L8bands*1000,SNRL8Ltyp,'-bx')
hold on
plot(L8bands*1000,SNRL8,'-ro')
plot(L8bands(2:end)*1000,SNRL7Ltyp,'-kv')
% plot(L8bands(1:5)*1000,ROIcranELM,'--bx')
% plot(L8bands(1:5)*1000,ROIlonELMbased,'-ro')
% plot(L8bands(1:5)*1000,ROIlonELM,'--ro')
% plot(L8bands(1:5)*1000,ROInearELMbased,'-mv')
% plot(L8bands(1:5)*1000,ROInearELM,'--mv')
% plot(L8bands(1:5)*1000,ROIoffELMbased,'-k+')
% plot(L8bands(1:5)*1000,ROIoffELM,'--k+')


% title('Off shore lake','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('SNR','fontsize',fs)
legend('L8 SNR at L typical','L8 SNR','L7 SNR')
set(gca,'fontsize',fs)

