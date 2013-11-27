function [a_SM,a_Chl,wavelength] = ParticleExtraction(Particles_curve,AfterMeth_curve,Vol)

% Particles_curve = SM010701;
% AfterMeth_curve = SM010702;

wavelength = Particles_curve(:,1);
ODfilt = Particles_curve(:,2)-Particles_curve(Particles_curve(:,1)==850,2); % direct from the filter w/o biass
ODfilt = smooth(ODfilt);

ODsusp = 0.378*ODfilt + 0.523*ODfilt.^2; % scattering correction

Ap_pigm = AfterMeth_curve(:,2)-AfterMeth_curve(AfterMeth_curve(:,1)==850,2); % after pigment extraction w/o biass
Ap_pigm = smooth(Ap_pigm);

ODp_pigm = 0.378*Ap_pigm + 0.523*Ap_pigm.^2; % scattering correction


Apigm = ODsusp - ODp_pigm;

%% SM and Chl Absorption Coefficients
V = 1.0000e-06*Vol; % m^3
r = (0.021)/2;
A = pi*r^2; % m: 21 mm


a_SM = 2.303*ODp_pigm/(V/A);
a_SM = smooth(a_SM);

a_Chl = 2.303*Apigm/(V/A);
a_Chl = smooth(a_Chl);


