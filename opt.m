% Optimization Routine
function XResults = opt(Ytest,LUT,LUTconc)
format long;
 
% Xtest: water pixels concentration from the image; Dim: 2000x3
% Ytest: water pixels reflectance from the image; Dim: 2000x8
% LUT: LUT


%% Reading in the LUT
Y = LUT;

%% Optimization
[f2,f1,f3] = ndgrid(unique(LUTconc(:,3))',...
    unique(LUTconc(:,2))',...
    unique(LUTconc(:,1))'); % CDOM SM CHL
options = optimset('Display','off','Tolfun',1e-10);

matlabpool open 4

XResults=zeros(size(Ytest,1),3);
%%

% h = waitbar(0,'Please wait...');
tic
parfor i = 1:size(Ytest,1)
    
a=sum((Y-repmat(Ytest(i,:),size(Ytest,1),1)).^2,2);
index=find(a==min(a));
x0=LUTconc(index,:);
if x0(1)==max(LUTconc(:,1)), x0(1)=46; end
if x0(2)==max(LUTconc(:,2)), x0(2)=20; end
if x0(3)==max(LUTconc(:,3)), x0(3)=12; end

% Y: LUT
% Ytest: TestSamples
% f1,f1,f3: Grid for the LUT
    XResults(i,:) = lsqnonlin(@myfun_mod,x0,...
        [min(LUTconc(:,1));min(LUTconc(:,2));min(LUTconc(:,3))],...
        [max(LUTconc(:,1));max(LUTconc(:,2));max(LUTconc(:,3))],...
        options,f1,f2,f3,Y,Ytest(i,:));
%     waitbar(i / size(Ytest,1))

    
end
disp('Elapsed time is (min):')
disp(toc/60)
  
  matlabpool close
% close(h)
