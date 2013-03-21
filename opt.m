% Optimization Routine
function XResults = opt(Ytest,LUT,LUTconc)
format long;
 
% Xtest: water pixels concentration from the image; Dim: 2000x3
% Ytest: water pixels reflectance from the image; Dim: 2000x8
% LUT: LUT


%% Reading in the LUT
Y = LUT;

%% Optimization

[f2,f1,f3] = ndgrid([0.25 0.50 0.75 1.0 2.0 4.0 7.0 10.0 12.0 14.0],...
    [0.25 0.50 1.0 2.0 4.0 8.0 10.0 14.0 20.0 24.0],...
    [0.25 0.50 1.0 3.0 5.0 7.0 12.0 24.0 46.0 68.0]); % CDOM SM CHL
options = optimset('Display','off','Tolfun',1e-10);

matlabpool open 4


x0=[.49 .48 .47]; % Starting point

XResults=zeros(size(Ytest,1),3);
%%

% h = waitbar(0,'Please wait...');
tic
parfor i = 1:size(Ytest,1)
    
a=sum((Y-repmat(Ytest(i,:),1000,1)).^2,2);
index=find(a==min(a));
x0=LUTconc(index,:);
if x0(1)==68, x0(1)=46; end
if x0(2)==24, x0(2)=20; end
if x0(3)==14, x0(3)=12; end

% Y: LUT
% Ytest: TestSamples
% f1,f1,f3: Grid for the LUT
    XResults(i,:) = lsqnonlin(@myfun_mod,x0,[0.25;0.25;0.25],[68;24;14],...
        options,f1,f2,f3,Y,Ytest(i,:));
%     waitbar(i / size(Ytest,1))

    
end
disp('Elapsed time is (min):')
disp(toc/60)
  
  matlabpool close
% close(h)
