cd /Users/javier/Desktop/Javier/PHD_RIT/LDCM/retrieval
%% Optical Densities
% TSS01 = load('./TSSAB01.ASC');
% TSS02 = load('./TSS02.ASC');

TSS01 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/SM130503/SM110501.ASC');
TSS02 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/InputOutput/SM130503/SM110502.ASC');

wavelength = TSS01(:,1);
ODfilt = TSS01(:,2)-TSS01(1,2); % direct from the filter w/o biass
ODfilt = smooth(ODfilt);

ODsusp = 0.378*ODfilt + 0.523*ODfilt.^2; % scattering correction

Ap_pigm = TSS02(:,2)-TSS02(1,2); % after pigment extraction w/o biass
Ap_pigm = smooth(Ap_pigm);

ODp_pigm = 0.378*Ap_pigm + 0.523*Ap_pigm.^2; % scattering correction


Apigm = ODsusp - ODp_pigm;



figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,ODfilt,'--b')
hold on
plot(wavelength,ODsusp,'b')

plot(wavelength,Ap_pigm,'--','Color',[.7 .5 0])
plot(wavelength,ODp_pigm,'Color',[.7 .5 0])

plot(wavelength,Apigm,'Color',[0 .5 0])

plot([wavelength(1) wavelength(end)],[0 0],'k')
title('Absorbance','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorbance','fontsize',fs)
set(gca,'fontsize',fs)
legend('ODfilt','ODsusp','Ap_pigm','ODp_pigm','Apigm')

%% SM and Chl Absorption Coefficients
V = 0.000050; % m^3 : 130 mL
A = 0.021; % m: 21 mm


a_SM = 2.303*ODp_pigm/(V/A);
a_SM = smooth(a_SM);

a_Chl = 2.303*Apigm/(V/A);
a_Chl = smooth(a_Chl);

figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,a_SM)
hold on
plot(wavelength,a_Chl,'Color',[0 .5 0])
title('Absorption Coefficient','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient m^{-1}','fontsize',fs)
set(gca,'fontsize',fs)
legend('SM','Chl')

%% CDOM absorption coefficient

CDOM1 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/CDOM130613/CD0101.ASC');
CDOM2 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/CDOM130613/CD1201.ASC');
CDOM3 = load('/Users/javier/Desktop/Javier/PHD_RIT/LDCM/WaterQualityProtocols/CDOM130613/CD1202.ASC');

wavelength = CDOM1(:,1);
ODcell1 = CDOM1(:,2)-CDOM1(1,2); % direct from the filter w/o biass
ODcell2 = CDOM2(:,2)-CDOM2(1,2); % direct from the filter w/o biass
ODcell3 = CDOM3(:,2)-CDOM3(1,2); % direct from the filter w/o biass

figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,ODcell1,'--b')
hold on
plot(wavelength,ODcell2,'--r')
plot(wavelength,ODcell3,'--r')
title('CDOM Absorbance','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorbance','fontsize',fs)
set(gca,'fontsize',fs)
legend('0101','1201','1202')
xlimits = get(gca,'xlim');
plot(xlimits,[0 0],'k')

L =0.1; % [m] (Cell length = 10cm)
a_CDOM1 = 2.303*ODcell1/L;
a_CDOM2 = 2.303*ODcell2/L;
a_CDOM3 = 2.303*ODcell3/L;

figure
fs = 15;
set(gcf,'color','white')
plot(wavelength,a_CDOM1,'b')
hold on
plot(wavelength,a_CDOM2,'r')
plot(wavelength,a_CDOM3,'r')
title('CDOM Absorption Coefficient','fontsize',fs)
xlabel('wavelength [nm]','fontsize',fs)
ylabel('Absorption Coefficient [m^{-1}]','fontsize',fs)
set(gca,'fontsize',fs)
legend('0101','1201','1202')
xlimits = get(gca,'xlim');
plot(xlimits,[0 0],'k')