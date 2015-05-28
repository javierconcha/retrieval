% Optimization Routine
function [XResults,residual] = opt(Ytest,LUT,LUTconc)
format long;

% Xtest: water pixels concentration from the image; Dim: 2000x3
% Ytest: water pixels reflectance from the image; Dim: 2000x8
% LUT: LUT

% global visual
% global visual2
%% Reading in the LUT
Y = LUT;

%% Optimization
CDOMconc = unique(LUTconc(:,3))';
SMconc   = unique(LUTconc(:,2))';
CHLconc  = unique(LUTconc(:,1))';
% [f2,f1,f3] = ndgrid(CDOMconc,SMconc,CHLconc); % CDOM SM CHL
options = optimset('Display','off','Tolfun',1e-10);

matlabpool open 4 % for using paralel computing

XResults=zeros(size(Ytest,1),3);
%%
tic
parfor i = 1:size(Ytest,1)
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
    
    % for the extremes (myfun_mod.m error otherwise)
    if x0(1)==max(LUTconc(:,1)), x0(1)=CHLconc(end-1); end
    if x0(2)==max(LUTconc(:,2)), x0(2)=SMconc(end-1); end
    if x0(3)==max(LUTconc(:,3)), x0(3)=CDOMconc(end-1); end
    
    % Y: LUT
    % Ytest: TestSamples
    % f1,f2,f3: Grid for the LUT
    [XResults(i,:),~,residual(i,:)]= ...
        lsqnonlin(@MyTrilinearInterp,x0,...
            [min(LUTconc(:,1));min(LUTconc(:,2));min(LUTconc(:,3))],...
            [max(LUTconc(:,1));max(LUTconc(:,2));max(LUTconc(:,3))],...
            options,LUT,Ytest(i,:),LUTconc);
    
end
disp('Elapsed time is (min):')
disp(toc/60)

matlabpool close % for using paralel computing
beep
pause(0.1)
beep
pause(0.1)
beep
pause(0.1)
beep

