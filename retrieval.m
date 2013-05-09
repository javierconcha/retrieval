% Retrieve the concentration from a image using a LUT from HL
cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/retrieval/
%% L5 image
im0115big = imread('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L5images/LT50160302000115AAA02/LT50160302000115AAA02radtif.tif');
im0115bigmask = imread('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L5images/LT50160302000115AAA02/LT50160302000115AAA02watermask.tif');
im0115bigmask(im0115bigmask>0)=1;

L5bands = [0.485,0.560,0.660,0.830,1.650,2.220];
%% for displaying
im0115R = im0115big(:,:,3);
im0115G = im0115big(:,:,2);
im0115B = im0115big(:,:,1);

im0115R = im0115R/max(im0115R(:));% normalize for maximum
im0115G = im0115G/max(im0115G(:));
im0115B = im0115B/max(im0115B(:));

im0115RGB(:,:,1) = histeq(im0115R);% normalize for maximum
im0115RGB(:,:,2) = histeq(im0115G);
im0115RGB(:,:,3) = histeq(im0115B);

figure
fs = 15;
set(gcf,'color','white')
imagesc(double(im0115RGB))
axis equal
hold on

% plot(2778,2123,'g.')
% plot(2814,2123,'g.')
% plot(2814,2159,'g.')
% plot(2778,2159,'g.')

xA = 2500;
xB = 2600;
yA = 2000;
yB = 2100;


plot(xA,yA,'g.')
plot(xB,yA,'g.')
plot(xB,yB,'g.')
plot(xA,yB,'g.')
%% water mask for big image
figure
imagesc(im0115bigmask)
colormap(gray)
axis equal

%% ROI middle of the lake
% RR = im0115big(2778:2814,2123:2159,:);
RR = im0115big(xA:xB,yA:yB,:);
ROImiddleOfTheLake = reshape(RR,[size(RR,1)*size(RR,2) size(RR,3)]);


figure
fs = 15;
set(gcf,'color','white')
plot(L5bands,ROImiddleOfTheLake')
title('ROI middle of the lake','fontsize',fs)
xlabel('wavelength [ \mum]','fontsize',fs)
ylabel('Radiance [W/m^2/sr/um]','fontsize',fs)
set(gca,'fontsize',fs)
%% Deglint
figure
fs = 15;
set(gcf,'color','white')
plot(ROImiddleOfTheLake(:,4),ROImiddleOfTheLake(:,1),'.')
title('ROI middle of the lake','fontsize',fs)
xlabel('NIR','fontsize',fs)
ylabel('Band 1','fontsize',fs)
set(gca,'fontsize',fs)
axis equal



figure
fs = 15;
set(gcf,'color','white')
plot(ROImiddleOfTheLake(:,4),ROImiddleOfTheLake(:,2),'.')
title('ROI middle of the lake','fontsize',fs)
xlabel('NIR Band','fontsize',fs)
ylabel('Band 2','fontsize',fs)
set(gca,'fontsize',fs)
axis equal

figure
fs = 15;
set(gcf,'color','white')
plot(ROImiddleOfTheLake(:,4),ROImiddleOfTheLake(:,3),'.')
title('ROI middle of the lake','fontsize',fs)
xlabel('NIR','fontsize',fs)
ylabel('Band 3','fontsize',fs)
set(gca,'fontsize',fs)
axis equal

figure
fs = 15;
set(gcf,'color','white')
plot(ROImiddleOfTheLake(:,4),ROImiddleOfTheLake(:,5),'.')
title('ROI middle of the lake','fontsize',fs)
xlabel('NIR','fontsize',fs)
ylabel('Band 5','fontsize',fs)
set(gca,'fontsize',fs)
axis equal

figure
fs = 15;
set(gcf,'color','white')
plot(ROImiddleOfTheLake(:,4),ROImiddleOfTheLake(:,6),'.')
title('ROI middle of the lake','fontsize',fs)
xlabel('NIR','fontsize',fs)
ylabel('Band 7','fontsize',fs)
set(gca,'fontsize',fs)
axis equal


[bi1,minNIR] = glint_removal(ROImiddleOfTheLake(:,1),ROImiddleOfTheLake(:,4));
[bi2,~] = glint_removal(ROImiddleOfTheLake(:,2),ROImiddleOfTheLake(:,4));
[bi3,~] = glint_removal(ROImiddleOfTheLake(:,3),ROImiddleOfTheLake(:,4));
[bi5,~] = glint_removal(ROImiddleOfTheLake(:,5),ROImiddleOfTheLake(:,4));
[bi7,~] = glint_removal(ROImiddleOfTheLake(:,6),ROImiddleOfTheLake(:,4));

bi1
bi2
bi3
bi5
bi7
minNIR


%% L5 water pixels 
imnew = reshape(im0115big,[size(im0115big,1)*size(im0115big,2) size(im0115big,3)]);
masknew = reshape(im0115bigmask,[size(im0115bigmask,1)*size(im0115bigmask,2) size(im0115bigmask,3)]);

waterpixelsL5 = imnew(masknew==1,:);
waterpixelsL5 = double(waterpixelsL5);

% %% display All water pixels
% 
% figure
% fs = 15;
% set(gcf,'color','white')
% plot(L5bands,waterpixelsL5(1:1000000,:))
% title('Reflectance water L5 image','fontsize',fs)
% xlabel('wavelength [\mu m]','fontsize',fs)
% ylabel('radiance [W/m^2/sr/um]','fontsize',fs)
% set(gca,'fontsize',fs)


%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% read tif
% filename = '/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L5images/LT50160302000115AAA02/LT50160302000115AAA02_ON130408resampledtif.twf';
% info = imfinfo(filename);
% t = Tiff(filename, 'r'); 
% % t.getTag('GeoKeyDirectoryTag')
% %%
% % Get the dimensions of the image to calculate coordinates.
% numRows = t.getTag('ImageLength');
% numCols = t.getTag('ImageWidth');
% % Get the ID number of the tile containing the coordinates.
% tileNum = t.computeStrip([numRows numCols],1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% L5 image cropped
% im0115crop = imread('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L5images/LT50160302000115AAA02/LT50160302000115AAA02_ONresampled.tif');
% im0115crop = imread('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L5images/LT50160302000115AAA02/LT50160302000115AAA02_ON130408resampledtif.tif');

% image with good ELM
im0115crop = imread('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L5images/LT50160302000115AAA02/LT50160302000115AAA02_ON130417resampledtif.tif');

im0115cropmask = imread('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L5images/LT50160302000115AAA02/LT50160302000115AAA02_ONresampled_mask008007tif.tif');
im0115cropmask(im0115cropmask>0)=1;

%% for displaying
im0115cropRGB(:,:,1)=imadjust(im0115crop(:,:,3));
im0115cropRGB(:,:,2)=imadjust(im0115crop(:,:,2));
im0115cropRGB(:,:,3)=imadjust(im0115crop(:,:,1));


impos = double(im0115cropRGB);
% impos(impos<0)=0;% only positive values

maskRGB(:,:,1)=double(im0115cropmask);
maskRGB(:,:,2)=double(im0115cropmask);
maskRGB(:,:,3)=double(im0115cropmask);

impos = impos.*maskRGB;


figure
imagesc(impos)
axis equal

%% mask
figure
imshow(imadjust(im0115cropmask))

%% water pixels. Convert each band in columns.

imnew = reshape(im0115crop,[size(im0115crop,1)*size(im0115crop,2) size(im0115crop,3)]);
masknew = reshape(im0115cropmask,[size(im0115cropmask,1)*size(im0115cropmask,2) size(im0115cropmask,3)]);

waterpixels = imnew(masknew==1,:);
waterpixels = double(waterpixels);
%% negative values
im = double(im0115crop);
imneg = zeros(size(im));
imneg(im<0)=im(im<0);% only negatives

imnegmask = zeros(size(im0115crop));% for displaying
imnegmask(im<0)=1;
imnegmask(im>=0)=0;
imnegmask = imnegmask+0.5*repmat(~double(im0115cropmask),[1 1 size(im0115crop,3)]); % for the land appear gray

figure
subplot(2,3,1)
fs = 15;
set(gcf,'color','white')
imshow(imnegmask(:,:,1))
title('band 1','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,3,2)
fs = 15;
set(gcf,'color','white')
imshow(imnegmask(:,:,2))
title('band 2','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,3,3)
fs = 15;
set(gcf,'color','white')
imshow(imnegmask(:,:,3))
title('band 3','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,3,4)
fs = 15;
set(gcf,'color','white')
imshow(imnegmask(:,:,4))
title('band 4','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,3,5)
fs = 15;
set(gcf,'color','white')
imshow(imnegmask(:,:,5))
title('band 5','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,3,6)
fs = 15;
set(gcf,'color','white')
imshow(imnegmask(:,:,6))
title('band 7','fontsize',fs)
set(gca,'fontsize',fs)

p=mtit('Negative Values',...
 	     'fontsize',fs+1,'xoff',0,'yoff',.025);
     
imnegb7 = imneg(:,:,6);  % to see negative values stats
mean(imnegb7(:))
std(imnegb7(:))
%% Stats water pixels
% % meanwp = mean(waterpixels,1);
% % stdwp = std(waterpixels,1);
% % maxwp = max(waterpixels,[],1);
% % minwp = min(waterpixels,[],1);
% % 
% % 
% % figure
% % fs = 15;
% % set(gcf,'color','white')
% % plot(L5bands,meanwp,'k')
% % title('Reflectance water L5 image','fontsize',fs)
% % xlabel('wavelength [\mu m]','fontsize',fs)
% % ylabel('reflectance','fontsize',fs)
% % set(gca,'fontsize',fs)
% % hold on
% % plot(L5bands,meanwp+stdwp,'g')
% % plot(L5bands,meanwp-stdwp,'g')
% % plot(L5bands,maxwp,'r')
% % plot(L5bands,minwp,'r')
% % xlim([min(L5bands) max(L5bands)])

%% display All water pixels

figure
fs = 15;
set(gcf,'color','white')
plot(L5bands,waterpixels(waterpixels(:,5)<0.08,:)')
title('Reflectance water L5 image','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
ylim([0 0.18])
%% Display negative values
waterpixels_neg = waterpixels(waterpixels(:,6)<0,:);

figure
fs = 15;
set(gcf,'color','white')
plot(L5bands,waterpixels_neg')
title('Reflectance water L5 image','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
ylim([0 0.18])


%% LUTs from HydroLight
LUT1temp = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/ReformedLUT.txt');
LUT2temp = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/ReformedLUT1.txt');
LUTconc = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/ReformedLUTConc.txt');
LUTconc = LUTconc';%ex: 1000x3


wl120 = LUT1temp(:,1);
wl140 = LUT2temp(:,1);

clear LUT1temp LUT2temp

LUT1 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/ReformedLUT1speclyb.txt');
LUT2 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/ReformedLUT2speclyb.txt');
LUT1 = LUT1';% ex: 1000x6
LUT2 = LUT2';% ex: 1000x6

% before 04/21/13
rr = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/LUTjavierL5.txt');
LUT3 = rr(:,2:end)';
LUTconc3 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/concentration_list.txt');

% LUT1000 = LUT3(LUTconc(:,2)<=24,:);
% LUT1000conc = LUTconc(LUTconc(:,2)<=24,:);
%



% new 04/21/13
rr = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/LUTL5130421.txt');
LUT1200 = rr(:,2:end)';
LUT1200conc = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/concentration_list130421.txt');

% new 05/08/13
rr = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/LUTL5130508.txt');
LUT1568 = rr(:,2:end)';
LUT1568conc = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/concentration_list130508.txt');



% LUT1000 = LUT1200(LUT1200conc(:,1)<=100,:);
% LUT1000conc = LUT1200conc(LUT1200conc(:,1)<=100,:);

cond = LUT1200conc(:,1)<=100;
LUT1000 = LUT1200(cond,:);
LUT1000conc = LUT1200conc(cond,:);


%%
figure
fs = 15;
set(gcf,'color','white')
plot(L5bands,LUT1568)
title('Reflectance LUT1000 - HydroLight','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
% % 
% % figure
% % fs = 15;
% % set(gcf,'color','white')
% % plot(L5bands,LUT2)
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
% plot(L5bands,waterdarkestred')
% hold on
% plot(L5bands,waterpixelspermag','c')
% plot(L5bands,LUTdarkestred','r')
% plot(L5bands,LUT1000(1,:)','k')
% 
% legend('L5 per b3','L5 per mag','HL:0.25;0.25;14','HL:0.25;0.25;0.25')
% title('Darkest pixel - L5 vs HL','fontsize',fs)
% xlabel('wavelength [\mu m]','fontsize',fs)
% ylabel('reflectance','fontsize',fs)
% set(gca,'fontsize',fs)
% % ylim([0 0.2])

%% 
disp('--------------------------------------------------------------------------')
disp('Running Optimization Routine')
    XResults = opt(waterpixels(:,1:4),LUT1000(:,1:4),LUT1000conc);
disp('Optimization Routine finished Successfully')
%% Test the optimization algorhythm
disp('--------------------------------------------------------------------------')
disp('Running Optimization Routine')
    XResultstest = opt(LUT2,LUT2,LUTconc);
disp('Optimization Routine finished Successfully')
%% Test the optimization algorhythm
LUTconc = LUT1000conc;
CDOMconc = unique(LUTconc(:,3))'
SMconc   = unique(LUTconc(:,2))'
CHLconc  = unique(LUTconc(:,1))'
disp('--------------------------------------------------------------------------')
disp('Running Optimization Routine')
    XResultstest = opt(LUT1000,LUT1000,LUTconc);
disp('Optimization Routine finished Successfully')
%% E_RMS
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
%% Display concentrations - data vs retrieved
figure
fs = 15;
set(gcf,'color','white')
plot(LUTconc(:,1),'r')
title('CHL Real vs retrieved','fontsize',fs)
xlabel('real','fontsize',fs)
ylabel('retrieved','fontsize',fs)
set(gca,'fontsize',fs)
hold on
plot(XResultstest(:,1))
legend('real','retrieved')

figure
fs = 15;
set(gcf,'color','white')
plot(LUTconc(:,2),'r')
title('SM Real vs retrieved','fontsize',fs)
xlabel('real','fontsize',fs)
ylabel('retrieved','fontsize',fs)
set(gca,'fontsize',fs)
hold on
plot(XResultstest(:,2))
legend('real','retrieved')

figure
fs = 15;
set(gcf,'color','white')
plot(LUTconc(:,3),'r')
title('CDOM Real vs retrieved','fontsize',fs)
xlabel('real','fontsize',fs)
ylabel('retrieved','fontsize',fs)
set(gca,'fontsize',fs)
hold on
plot(XResultstest(:,3))
legend('real','retrieved')

%% Mapping Concentrations

ConcRet = zeros(size(masknew),3);
ConcRet(masknew==1,:) = XResults; % Concentration Retrieved

CHLmap  = reshape(ConcRet(:,1),[size(im0115cropmask,1) size(im0115cropmask,2) size(im0115cropmask,3)]);
SMmap   = reshape(ConcRet(:,2),[size(im0115cropmask,1) size(im0115cropmask,2) size(im0115cropmask,3)]);
CDOMmap = reshape(ConcRet(:,3),[size(im0115cropmask,1) size(im0115cropmask,2) size(im0115cropmask,3)]);

figure
fs = 15;
set(gcf,'color','white')
imagesc(CHLmap)
title('CHL map','fontsize',fs)
set(gca,'fontsize',fs)
colorbar

figure
fs = 15;
set(gcf,'color','white')
imagesc(SMmap)
title('SM map','fontsize',fs)
set(gca,'fontsize',fs)
colorbar

figure
fs = 15;
set(gcf,'color','white')
imagesc(CDOMmap)
title('CDOM map','fontsize',fs)
set(gca,'fontsize',fs)
colorbar

%% Comparison between two LUTs from Aaron and one curve from HL5.1 in tropos

r_02502514 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/r_025025014.txt');
figure
plot(LUT1(10,:))
hold on
plot(LUT1000(10,:),'k')
% plot(r_02502514(:,2),'r')

