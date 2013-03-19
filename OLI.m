%%%%% Program to find the constituent concentrations based on reflectance
%%%%% values
%
% Created on 07/31/2011
%

clear all;
format long;
%% READ IN LUT/SAMPLES/OLI RESPONSE FUNCTION
SpectralLUT=load('C:\Users\adgpci\Desktop\BobbyCR\ReformedLUT.txt');
% SpectralLUT=load('C:\Users\adgpci\Desktop\BobbyCR\LUT_OffNadir10.txt');
SpectralSamples=load('C:\Users\adgpci\Desktop\BobbyCR\ReformedSamples.txt');
% SpectralSamples=load('C:\Users\adgpci\Desktop\BobbyCR\Samples_OffNadir10.txt');

SpectralResponse=load('C:\Users\adgpci\Desktop\BobbyCR\SensorResponses\OLI_VNIR_Response.txt');

%% DEFINE CONSTITUENT CONCENTRATIONS
LUT_Conc=(load('C:\Users\adgpci\Desktop\BobbyCR\ReformedLUTConc.txt'))';
Samples_Conc=(load('C:\Users\adgpci\Desktop\BobbyCR\ReformedSamplesConc.txt'))';
%% INTERPOLATE AND SPECTRALLY SAMPLE LUT, SAMPLE DATA, AND ATMOSPHERE
% PROPAGATE SIGNALS TO TOA
% Atmosphere=load('C:\Users\adgpci\Desktop\BobbyCR\Atmospheres\tape7_SenZen_Minus225.txt');
% Atmosphere=load('C:\Users\adgpci\Desktop\SenZen\Nadir\tape7.txt');
% Atmosphere=load('C:\Users\adgpci\Desktop\BobbyCR\Atmospheres\Standard\tape7_2bounces.txt');
Atmosphere=load('C:\Users\adgpci\Desktop\MostRecentTape7Urban_Minus20.txt');
Atmosphere(:,3)=Atmosphere(:,3)/100.0;
% AtmosphereUsedForRetrieval=load('C:\Users\adgpci\Desktop\BobbyCR\Atmospheres\tape7_SenZen_Minus225.txt');
% AtmosphereUsedForRetrieval=load('C:\Users\adgpci\Desktop\SenZen\Nadir\tape7.txt');
% AtmosphereUsedForRetrieval=load('C:\Users\adgpci\Desktop\BobbyCR\Atmospheres\Standard\tape7_2bounces.txt');
AtmosphereUsedForRetrieval=load('C:\Users\adgpci\Desktop\MostRecentTape7Urban_Minus20.txt');
AtmosphereUsedForRetrieval(:,3)=AtmosphereUsedForRetrieval(:,3)/100.0;

InterpolatedSpectralLUT=interp1(SpectralLUT(:,1),SpectralLUT(:,2:1001),SpectralResponse(:,1));
InterpolatedSpectralLUT(find(isnan(InterpolatedSpectralLUT)))=0;            %  THIS DEALS WITH NaNs IN THE DATA

InterpolatedAtmosphere=interp1(AtmosphereUsedForRetrieval(:,1),AtmosphereUsedForRetrieval(:,2:3),SpectralResponse(:,1));
InterpolatedAtmosphere(find(isnan(InterpolatedAtmosphere)))=0;            %  THIS DEALS WITH NaNs IN THE DATA

InterpolatedActualAtmosphere=interp1(Atmosphere(:,1),Atmosphere(:,2:3),SpectralResponse(:,1));
InterpolatedActualAtmosphere(find(isnan(InterpolatedActualAtmosphere)))=0;            %  THIS DEALS WITH NaNs IN THE DATA

TOA_Samples=SpectralSamples(:,2:2001).*repmat(Atmosphere(:,2),1,2000)+repmat(Atmosphere(:,3),1,2000);

InterpolatedSpectralSamples=interp1(SpectralSamples(1:120,1),TOA_Samples,SpectralResponse(:,1));
InterpolatedSpectralSamples(find(isnan(InterpolatedSpectralSamples)))=0;    %  THIS DEALS WITH NaNs IN THE DATA

% SPECTRALLY SAMPLE
for n = 2:6
    SpecSampledLUT(:,n-1) = (sum(InterpolatedSpectralLUT.*repmat(SpectralResponse(:,n),1,1000),1)./sum(repmat(SpectralResponse(:,n),1,1000),1));
    SpecSampledData(:,n-1) = (sum(InterpolatedSpectralSamples.*repmat(SpectralResponse(:,n),1,2000),1)./sum(repmat(SpectralResponse(:,n),1,2000),1));
    SpecSampledAtmosphere(:,n-1) = (sum(InterpolatedAtmosphere.*repmat(SpectralResponse(:,n),1,2),1)./sum(repmat(SpectralResponse(:,n),1,2),1));
    SpecSampledActualAtmosphere(:,n-1) = (sum(InterpolatedActualAtmosphere.*repmat(SpectralResponse(:,n),1,2),1)./sum(repmat(SpectralResponse(:,n),1,2),1));
end

RetrievedSS_Data=(SpecSampledData-repmat(SpecSampledAtmosphere(2,:),2000,1))./repmat(SpecSampledAtmosphere(1,:),2000,1);

%% ADD QUANTIZATION
StepSize=[0.000135498 0.000141846 0.000132813 0.000112793 .0000686035];                             % W/m^2/ster/nm

StepSizeArray=repmat(StepSize,2000,1);
QuantizedData = round(SpecSampledData./StepSizeArray).*StepSizeArray;
RetrievedQuantized_Data=(QuantizedData-repmat(SpecSampledAtmosphere(2,:),2000,1))./repmat(SpecSampledAtmosphere(1,:),2000,1);

%% ADD NOISE
% SNR=[130 130 100 90 90];
SNR=[178.75 237.5 192.5 155.25 144];
% SNR=[227.5 344.5 285 220.5 198];
ROI=1;
for n = 1:2000
    test=repmat(SpecSampledData(n,:),ROI,1);
    Noise=test./repmat(SNR,ROI,1);                                %  SNR=S/N   -->  N=S/SNR
    RandomNumbers=randn(ROI,5);
    NoiseData(n,:)=mean((test+Noise.*RandomNumbers),1);
end
NoiseQuantizedData = round(NoiseData./StepSizeArray).*StepSizeArray;
RetrievedNoiseActualData=(NoiseQuantizedData-repmat(SpecSampledActualAtmosphere(2,:),2000,1))./repmat(SpecSampledActualAtmosphere(1,:),2000,1);
RetrievedNoise_Data=(NoiseQuantizedData-repmat(SpecSampledAtmosphere(2,:),2000,1))./repmat(SpecSampledAtmosphere(1,:),2000,1);



%% OPTIMIZATION STUFF
[f2,f1,f3] = ndgrid([0.25 0.50 0.75 1.0 2.0 4.0 7.0 10.0 12.0 14.0],[0.25 0.50 1.0 2.0 4.0 8.0 10.0 14.0 20.0 24.0],[0.25 0.50 1.0 3.0 5.0 7.0 12.0 24.0 46.0 68.0]);
options = optimset('Display','on','LevenbergMarquardt','on');
options.TolFun =[0.0000000001]

RetrievedSamples=Samples_Conc;
SampleData=RetrievedNoise_Data;
for i = 1:size(RetrievedSamples,1)
%     f = myfun_mod(x0,X1,X2,X3,Ys,Yn)
%     plot(Samples_Data(i,:), 'Color', 'red')
%     hold on
    a=sum((SpecSampledLUT-repmat(SampleData(i,:),1000,1)).^2,2);
    index=find(a==min(a));
    x0=LUT_Conc(index,:);
    if x0(1)==68, x0(1)=46; end
    if x0(2)==24, x0(2)=20; end
    if x0(3)==14, x0(3)=12; end
    RetrievedSamples(i,:) = lsqnonlin(@myfun_mod,x0,[0.25;0.25;0.25],[68.0;24.0;14.0],options,f1,f2,f3,SpecSampledLUT,SampleData(i,:));
    i
end
%% ERROR REPRESENTATION
figure(1)
plot(RetrievedSamples(:,1), Samples_Conc(:,1), 'ro')
hold on;

figure(2)
plot(RetrievedSamples(:,2), Samples_Conc(:,2), 'go')
hold on;

figure(3)
plot(RetrievedSamples(:,3), Samples_Conc(:,3), 'bo')
hold on;

RetrievalErrors=mean(abs(RetrievedSamples-Samples_Conc))
Samples_Conc(find((Samples_Conc==0)))=0.0001;

AbsoluteRetrievalErrors=mean((abs(RetrievedSamples-Samples_Conc))./Samples_Conc)
PercentRetrievalErrors=[RetrievalErrors(:,1)/68.0,RetrievalErrors(:,2)/24.0,RetrievalErrors(:,3)/14.0]

