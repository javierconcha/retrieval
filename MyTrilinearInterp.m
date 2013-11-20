function f = MyTrilinearInterp(x0,LUT,Ytest,LUTconc)
% Trilinear interpolation
% By Javier A. Concha
% 05-09-13
% Trilinear interpolation for non uniform and monotonig grid
%% Intialization
% Initial Concentration

nx = x0(1); 
ny = x0(2);
nz = x0(3);

% Concentrations per components
xx = unique(LUTconc(:,1)); % CHL
yy = unique(LUTconc(:,2)); % SM
zz = unique(LUTconc(:,3)); % CDOM
%% To find index in between xx and weight for each dim
index_up = find(xx >= nx,1);

if isempty(index_up)
    index_up = length(xx);
end

index_low = index_up-1;

if index_low ~= 0
    xl = xx(index_low); % lower concentration
else xl = 0;
end
xu = xx(index_up);      % upper concentration

wx = (nx-xl)/(xu-xl);
if isnan(wx)
    wx = 1;
end
%% To find index in between yy and weight for each dim
indey_up = find(yy >= ny,1);

if isempty(indey_up)
    indey_up = length(yy);
end

indey_low = indey_up-1;

if indey_low ~= 0
    yl = yy(indey_low);
else yl = 0;
end
yu = yy(indey_up);


wy = (ny-yl)/(yu-yl);
if isnan(wy)
    wy = 1;
end
%% To find index in between zz and weight for each dim
indez_up = find(zz >= nz,1);

if isempty(indez_up)
    indez_up = length(zz);
end

indez_low = indez_up-1;

if indez_low ~= 0
    zl = zz(indez_low);
else zl = 0;
end
zu = zz(indez_up);

wz = (nz-zl)/(zu-zl);
if isnan(wz)
    wz = 1;
end
%% Look up the values of the 8 points surrounding the cube
if (index_low ~= 0 && indey_low ~= 0 && indez_low ~= 0)
    V000 = LUT(LUTconc(:,1)==xl & ...
        LUTconc(:,2)==yl & ...
        LUTconc(:,3)==zl,:);
else
    V000 = zeros(1,size(LUT,2));
end

if (index_low ~= 0 && indey_low ~= 0)
    V001 = LUT(LUTconc(:,1)==xl & ...
        LUTconc(:,2)==yl & ...
        LUTconc(:,3)==zu,:);
else
    V001 = zeros(1,size(LUT,2));
end

if (index_low ~= 0 && indez_low ~= 0)
    V010 = LUT(LUTconc(:,1)==xl & ...
        LUTconc(:,2)==yu & ...
        LUTconc(:,3)==zl,:);
else
    V010 = zeros(1,size(LUT,2));
end

if index_low ~= 0
    V011 = LUT(LUTconc(:,1)==xl & ...
        LUTconc(:,2)==yu & ...
        LUTconc(:,3)==zu,:);
else
    V011 = zeros(1,size(LUT,2));
end

if (indey_low ~= 0 && indez_low ~= 0)
    V100 = LUT(LUTconc(:,1)==xu & ...
        LUTconc(:,2)==yl & ...
        LUTconc(:,3)==zl,:);
else
    V100 = zeros(1,size(LUT,2));
end

if indey_low ~= 0
    V101 = LUT(LUTconc(:,1)==xu & ...
        LUTconc(:,2)==yl & ...
        LUTconc(:,3)==zu,:);
else
    V101 = zeros(1,size(LUT,2));
end

if indez_low ~= 0
    V110 = LUT(LUTconc(:,1)==xu & ...
        LUTconc(:,2)==yu & ...
        LUTconc(:,3)==zl,:);
else
    V110 = zeros(1,size(LUT,2));
end

V111 = LUT(LUTconc(:,1)==xu & ...
    LUTconc(:,2)==yu & ...
    LUTconc(:,3)==zu,:);

Vxyz =  ...
    V000*(1-wx)*(1-wy)*(1-wz) +...
    V001*(1-wx)*(1-wy)*wz + ...
    V010*(1-wx)*wy*(1-wz) + ...
    V011*(1-wx)*wy*wz + ...
    V100*wx*(1-wy)*(1-wz) + ...
    V101*wx*(1-wy)*wz + ...
    V110*wx*wy*(1-wz) + ...
    V111*wx*wy*wz;

f = Ytest - Vxyz;


% figure(30)
% plot(Ytest,'r')
% hold on
% plot(Vxyz,'g')
% plot(V000)
% plot(V001)
% plot(V010)
% plot(V011)
% plot(V100)
% plot(V101)
% plot(V110)
% plot(V111)
