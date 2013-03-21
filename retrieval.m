cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/retrieval/

%% L5 image
im0115 = imread('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L5images/LT50160302000115AAA02/LT50160302000115AAA02_ONresampled.tif');
im0115mask = imread('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/L5images/LT50160302000115AAA02/LT50160302000115AAA02_ONresampled_masktif.tif');
im0115mask(im0115mask>0)=1;

format long;
%% for displaying
im0115RGB(:,:,1)=im0115(:,:,3)/max(max(max(im0115(:,:,3))));% normalize for maximum
im0115RGB(:,:,2)=im0115(:,:,2)/max(max(max(im0115(:,:,2))));
im0115RGB(:,:,3)=im0115(:,:,1)/max(max(max(im0115(:,:,1))));

minRGB = double(min(min(min(im0115RGB))));
maxRGB = double(max(max(max(im0115RGB))));

impos = double(im0115RGB);% only positive values
impos(impos<0)=0;

maskRGB(:,:,1)=double(im0115mask);
maskRGB(:,:,2)=double(im0115mask);
maskRGB(:,:,3)=double(im0115mask);

impos = impos.*maskRGB;
impos(:,:,1) = imadjust(impos(:,:,1));
impos(:,:,2) = imadjust(impos(:,:,2));
impos(:,:,3) = imadjust(impos(:,:,3));


figure
imshow(impos)

%% mask
figure
imshow(imadjust(im0115mask))

%%
L5bands = [0.485,0.560,0.660,0.830,1.650,2.220];

imnew = reshape(im0115,[size(im0115,1)*size(im0115,2) size(im0115,3)]);
masknew = reshape(im0115mask,[size(im0115mask,1)*size(im0115mask,2) size(im0115mask,3)]);

waterpixels = imnew(masknew==1,:);
waterpixels = double(waterpixels);
%% negative values
im = double(im0115);% only positive negatives
imneg = zeros(size(im0115));
imneg(im<0)=1;
imneg(im>=0)=0;
imneg = imneg+0.5*repmat(~double(im0115mask),[1 1 size(im0115,3)]);

figure
subplot(2,3,1)
fs = 15;
set(gcf,'color','white')
imshow(imneg(:,:,1))
title('band 1','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,3,2)
fs = 15;
set(gcf,'color','white')
imshow(imneg(:,:,2))
title('band 2','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,3,3)
fs = 15;
set(gcf,'color','white')
imshow(imneg(:,:,3))
title('band 3','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,3,4)
fs = 15;
set(gcf,'color','white')
imshow(imneg(:,:,4))
title('band 4','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,3,5)
fs = 15;
set(gcf,'color','white')
imshow(imneg(:,:,5))
title('band 5','fontsize',fs)
set(gca,'fontsize',fs)

subplot(2,3,6)
fs = 15;
set(gcf,'color','white')
imshow(imneg(:,:,6))
title('band 7','fontsize',fs)
set(gca,'fontsize',fs)

p=mtit('Negative Values',...
 	     'fontsize',fs+1,'xoff',0,'yoff',.025);
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

%% LUTs from Aaron
LUT1temp = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/ReformedLUT.txt');
LUT2temp = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/ReformedLUT1.txt');
LUTconc = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/ReformedLUTConc.txt');
LUTconc = LUTconc';

wl120 = LUT1temp(:,1);
wl140 = LUT2temp(:,1);

clear LUT1temp LUT2temp

LUT1 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/ReformedLUT1speclyb.txt');
LUT2 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/HLinout/ReformedLUT2speclyb.txt');
LUT1 = LUT1';
LUT2 = LUT2';


% % figure
% % fs = 15;
% % set(gcf,'color','white')
% % plot(L5bands,LUT1)
% % title('Reflectance LUT1 - HydroLight','fontsize',fs)
% % xlabel('wavelength [\mu m]','fontsize',fs)
% % ylabel('reflectance','fontsize',fs)
% % set(gca,'fontsize',fs)
% % 
% % figure
% % fs = 15;
% % set(gcf,'color','white')
% % plot(L5bands,LUT2)
% % title('Reflectance LUT2 - HydroLight','fontsize',fs)
% % xlabel('wavelength [\mu m]','fontsize',fs)
% % ylabel('reflectance','fontsize',fs)
% % set(gca,'fontsize',fs)

%% 
disp('--------------------------------------------------------------------------')
disp('Running Optimization Routine')
    XResults = opt(waterpixels,LUT1,LUTconc);
disp('Optimization Routine finished Successfully')
%% Test the optimization algorhythm
disp('--------------------------------------------------------------------------')
disp('Running Optimization Routine')
    XResultstest = opt(LUT1,LUT1,LUTconc);
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

CHLmap  = reshape(ConcRet(:,1),[size(im0115mask,1) size(im0115mask,2) size(im0115mask,3)]);
SMmap   = reshape(ConcRet(:,2),[size(im0115mask,1) size(im0115mask,2) size(im0115mask,3)]);
CDOMmap = reshape(ConcRet(:,3),[size(im0115mask,1) size(im0115mask,2) size(im0115mask,3)]);

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
plot(LUT2(10,:),'k')
plot(r_02502514(:,2),'r')




