% Optimization Routine
function [XResults,residual] = opt(Ytest,LUT,LUTconc,InputLabel)
format short;

% Xtest: water pixels concentration from the image; Dim: 2000x3
% Ytest: water pixels reflectance from the image; Dim: 2000x8
% LUT: LUT

% global visual
% global visual2
%% Reading in the LUT
Y = LUT;

%% Optimization
% [f2,f1,f3] = ndgrid(CDOMconc,SMconc,CHLconc); % CDOM SM CHL
options = optimset('Display','off','Tolfun',1e-10);

% matlabpool open 4 % for using paralel computing

XResults = zeros(size(Ytest,1),3);
residual = zeros(size(Ytest));
%%
tic
for i = 1:size(Ytest,1)
%     if visual==1
% figure(30)
% clf
% 
% i
% if i==172
%    disp('i=172') 
% end

%     figure(68)
%     clf
%     ylim([0 0.05])
%     end 
%     
%     if visual2==1
%     figure(69)
%     xlim([0 68])
%     ylim([0 24])
%     zlim([0 14])
%     end
    

    % select x0
    a=sum((Y-repmat(Ytest(i,:),size(Y,1),1)).^2,2);
    [~,index]=min(a);
    x0=LUTconc(index,:);
    % Select from the LUT with same input. From ponds OR lake inputs
    CDconc  = unique(LUTconc(strcmp(InputLabel,InputLabel(index)),3))';
    SMconc  = unique(LUTconc(strcmp(InputLabel,InputLabel(index)),2))';
    CHconc  = unique(LUTconc(strcmp(InputLabel,InputLabel(index)),1))';
    
    LUTconcUsed = LUTconc(strcmp(InputLabel,InputLabel(index)),:);
    
    % for the extremes (myfun_mod.m error otherwise)
    if x0(1)==max(CHconc) , x0(1)=CHconc(end-1); end
    if x0(2)==max(SMconc) , x0(2)=SMconc(end-1); end
    if x0(3)==max(CDconc) , x0(3)=CDconc(end-1); end
    
    % Y: LUT
    % Ytest: TestSamples
    % f1,f2,f3: Grid for the LUT
    [XResults(i,:),~,residual(i,:)]= ...
        lsqnonlin(@MyTrilinearInterp,x0,...
            [min(CHconc);min(SMconc);min(CDconc)],...
            [max(CHconc);max(SMconc);max(CDconc)],...
            options,Y,Ytest(i,:),LUTconcUsed);
    
end
disp('Elapsed time is (min):')
disp(toc/60)

% matlabpool close % for using paralel computing
beep
pause(0.1)
beep
pause(0.1)
beep
pause(0.1)
beep

