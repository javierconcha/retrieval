function f = mytriinterp(x0,LUT,Ytest,LUTconc)
%%
nx = x0(1);
ny = x0(2);
nz = x0(3);
%% To find index in between xx and weight for each dim
index_up = find(xx >= nx,1);
index_low = index_up-1;

if index_low ~= 0
    xl = xx(index_low);
else xl = 0;
end
xu = xx(index_up);

wx = (nx-xl)/(xu-xl);
%% To find index in between yy and weight for each dim
indey_up = find(yy >= ny,1);
indey_low = indey_up-1;

if indey_low ~= 0
    yl = yy(indey_low);
else yl = 0;
end
yu = yy(indey_up);


wy = (ny-yl)/(yu-yl);
%% To find index in between zz and weight for each dim
indez_up = find(zz >= nz,1);
indez_low = indez_up-1;

if indez_low ~= 0
    zl = zz(indez_low);
else zl = 0;
end
zu = zz(indez_up);

wz = (nz-zl)/(zu-zl);
%% Look up the values of the 8 points surrounding the cube
if index_low ~= 0 && indey_low ~= 0 && indez_low ~= 0
    V000 = LUT(LUTconc(:,1)==xl & ...
        LUTconc(:,2)==yl & ...
        LUTconc(:,3)==zl,:);
else
    V000 = zeros(1,size(LUT,2));
end

if index_low ~= 0 && indey_low ~= 0
    V001 = LUT(LUTconc(:,1)==xl & ...
        LUTconc(:,2)==yl & ...
        LUTconc(:,3)==zu,:);
else
    V001 = zeros(1,size(LUT,2));
end

if index_low ~= 0 && indez_low ~= 0
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

if indey_low ~= 0 && indez_low ~= 0
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
