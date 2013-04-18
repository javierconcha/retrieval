cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/retrieval/
%% L5 image
im0115big = imread('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L5images/LT50160302000115AAA02/LT50160302000115AAA02.tif');
im0115bigmask = imread('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L5images/LT50160302000115AAA02/LT50160302000115AAA02watermask.tif');
im0115bigmask(im0115bigmask>0)=1;

L5bands = [0.485,0.560,0.660,0.830,1.650,2.220];
%% for displaying
im0115RGB(:,:,1)=imadjust(im0115big(:,:,3));%/max(max(max(im0115big(:,:,3))));% normalize for maximum
im0115RGB(:,:,2)=imadjust(im0115big(:,:,2));%/max(max(max(im0115big(:,:,2))));
im0115RGB(:,:,3)=imadjust(im0115big(:,:,1));%/max(max(max(im0115big(:,:,1))));

figure
imagesc(im0115RGB,[2 100])
axis equal

figure
imagesc(im0115bigmask)
colormap(gray)
axis equal

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% L5 image cropped
% im0115crop = imread('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L5images/LT50160302000115AAA02/LT50160302000115AAA02_ONresampled.tif');
im0115crop = imread('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L5images/LT50160302000115AAA02/LT50160302000115AAA02_ON130408resampledtif.tif');



im0115cropmask = imread('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L5images/LT50160302000115AAA02/LT50160302000115AAA02_ONresampled_masktif.tif');
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
plot(L5bands,waterpixels')
title('Reflectance water L5 image','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
% ylim([0 0.2])
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
% ylim([0 0.2])


%% LUTs from Aaron
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

rr = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/LUTjavierL5.txt');
LUT3 = rr(:,2:end)';
LUTconc3 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/concentration_list.txt');

LUT1000 = LUT3(LUTconc(:,2)<=24,:);
LUT1000conc = LUTconc(LUTconc(:,2)<=24,:);

figure
fs = 15;
set(gcf,'color','white')
plot(L5bands,LUT1000)
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
%% Darkest in the red
% Darkest in red band
[~,darkestindex] = min(waterpixels(:,3));
waterdarkestred = waterpixels(darkestindex,:);

[~,LUTdarkestindex] = min(LUT1000(:,3));
LUTdarkestred = LUT1000(LUTdarkestindex,:);
LUT1000conc(LUTdarkestindex,:)

% Darkest in Magnitude
wpnorm = sqrt(sum(waterpixels.^2,2)); % for the waterpixels
[normmin,index] = min(wpnorm);
waterpixelspermag = waterpixels(index,:);

LUT1000norm = sqrt(sum(LUT1000.^2,2)); % for the LUT
[normmin,index] = min(LUT1000norm);
LUT1000conc(index,:)





figure
fs = 15;
set(gcf,'color','white')
plot(L5bands,waterdarkestred')
hold on
plot(L5bands,waterpixelspermag','c')
plot(L5bands,LUTdarkestred','r')
plot(L5bands,LUT1000(1,:)','k')

legend('L5 per b3','L5 per mag','HL:0.25;0.25;14','HL:0.25;0.25;0.25')
title('Darkest pixel - L5 vs HL','fontsize',fs)
xlabel('wavelength [\mu m]','fontsize',fs)
ylabel('reflectance','fontsize',fs)
set(gca,'fontsize',fs)
% ylim([0 0.2])

%% 
disp('--------------------------------------------------------------------------')
disp('Running Optimization Routine')
    XResults = opt(waterpixels,LUT1000,LUTconc);
disp('Optimization Routine finished Successfully')
%% Test the optimization algorhythm
disp('--------------------------------------------------------------------------')
disp('Running Optimization Routine')
    XResultstest = opt(LUT2,LUT2,LUTconc);
disp('Optimization Routine finished Successfully')
%% Test the optimization algorhythm
disp('--------------------------------------------------------------------------')
disp('Running Optimization Routine')
    XResultstest = opt(LUT1000,LUT1000,LUTconc);
disp('Optimization Routine finished Successfully')
%% E_RMS
E_Chl = sqrt(sum((XResultstest(:,1)-LUTconc(:,1)).^2)/size(XResultstest,1))
E_Chl*100/68
E_SM = sqrt(sum((XResultstest(:,2)-LUTconc(:,2)).^2)/size(XResultstest,1))
E_SM*100/24
E_CDOM = sqrt(sum((XResultstest(:,3)-LUTconc(:,3)).^2)/size(XResultstest,1))
E_CDOM*100/14
%% Display data vs retrieved

figure
fs = 15;
set(gcf,'color','white')
plot(LUTconc(:,1),XResultstest(:,1),'.')
title('CHL Real vs retrieved','fontsize',fs)
xlabel('real','fontsize',fs)
ylabel('retrieved','fontsize',fs)
set(gca,'fontsize',fs)

figure
fs = 15;
set(gcf,'color','white')
plot(LUTconc(:,2),XResultstest(:,2),'.')
title('SM Real vs retrieved','fontsize',fs)
xlabel('real','fontsize',fs)
ylabel('retrieved','fontsize',fs)
set(gca,'fontsize',fs)

figure
fs = 15;
set(gcf,'color','white')
plot(LUTconc(:,3),XResultstest(:,3),'.')
title('CDOM Real vs retrieved','fontsize',fs)
xlabel('real','fontsize',fs)
ylabel('retrieved','fontsize',fs)
set(gca,'fontsize',fs)
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

%% 

ROImiddlelakemean = [40.575298 20.721369 10.058275 3.807572 0 0];

ROImiddlelakestdv = [0.750423 0.904410 0.590065 0.265090 0 0];

darkpxmiddlelake = ROImiddlelakemean-2*ROImiddlelakestdv;
