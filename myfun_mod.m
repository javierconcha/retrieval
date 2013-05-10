function f = myfun_mod(x0,f1,f2,f3,LUT,Ytest,LUTconc)
% It uses tripolar interpolation for monotonic data (uniformly distributed grid) 
% x0: starting point
% f1,f2,f3: Grid for the LUT

LUTconc(1,:)
global F
% global visual
% global visual2

% Starting Point
xi = x0(1); % CHL
yi = x0(2); % SM
zi = x0(3); % CDOM

xx = unique(f3);           % CHL concentrations
yy = unique(f1);           % SM concentrations
zz = unique(f2);           % CDOM concentrations

% Determine the nearest location of xi in f1
[~,i] = sort([xx;xi]);
si(i) = (1:length(i));
si = (si(length(xx)+1:end)-1)';

% Map values in xi to index offset (si) via linear interpolation
si(si<1) = 1;
si(si>length(xx)-1) = length(xx)-1;
si = si + (xi-xx(si))./(xx(si+1)-xx(si));

% Determine the nearest location of yi in f2
[~,i] = sort([yy;yi]);
ti(i) = (1:length(i));
ti = (ti(length(yy)+1:end)-1)';

% Map values in yi to index offset (ti) via linear interpolation
ti(ti<1) = 1;
ti(ti>length(yy)-1) = length(yy)-1;
ti = ti + (yi-yy(ti))./(yy(ti+1)-yy(ti));

% Determine the nearest location of zi in f3
[~,i] = sort([zz;zi]);
wi(i) = (1:length(i));
wi = (wi(length(zz)+1:end)-1)';

% Map values in zi to index offset (wi) via linear interpolation
wi(wi<1) = 1;
wi(wi>length(zz)-1) = length(zz)-1;
wi = wi + (zi-zz(wi))./(zz(wi+1)-zz(wi));

[x,y,z] = meshgrid(ones(class(f1)):size(f1,2),...
    ones(superiorfloat(f2,f3)):size(f2,1),1:size(f3,3)); % grid of values 1 to 10

xi = si; yi = ti; zi = wi; % xi,si,ti are integers

nrows = size(f1,1);
ncols = size(f1,2);
npages = size(f1,3);
mx = numel(x); my = numel(y); mz = numel(z);
s = 1 + (xi-x(1))/(x(mx)-x(1))*(ncols-1);
t = 1 + (yi-y(1))/(y(my)-y(1))*(nrows-1);
w = 1 + (zi-z(1))/(z(mz)-z(1))*(npages-1);

% Check for out of range values of s and set to 1
sout = find((s<1)|(s>ncols));
if ~isempty(sout), s(sout) = ones(size(sout)); end

% Check for out of range values of t and set to 1
tout = find((t<1)|(t>nrows));
if ~isempty(tout), t(tout) = ones(size(tout)); end

% Check for out of range values of w and set to 1
wout = find((w<1)|(w>npages));
if ~isempty(wout), w(wout) = ones(size(wout)); end

% Matrix element indexing
nw = nrows*ncols;
ndx = floor(w)+floor(t-1)*nrows+floor(s-1)*nw;

% Compute intepolation parameters, check for boundary value.
if isempty(s), d = s; else d = find(s==ncols); end
s(:) = (s - floor(s));
if ~isempty(d), s(d) = s(d)+1; ndx(d) = ndx(d)-nrows; end

% Compute intepolation parameters, check for boundary value.
if isempty(t), d = t; else d = find(t==nrows); end
t(:) = (t - floor(t));
if ~isempty(d), t(d) = t(d)+1; ndx(d) = ndx(d)-1; end

% Compute intepolation parameters, check for boundary value.
if isempty(w), d = w; else d = find(w==npages); end
w(:) = (w - floor(w));
if ~isempty(d), w(d) = w(d)+1; ndx(d) = ndx(d)-nw; end

% size(LUT)
% ndx
% ndx+1
% ndx+nrows
% ndx+(nrows+1)
% ndx+nw
% ndx+1+nw
% ndx+nrows+nw
% ndx+(nrows+1+nw)

F =  (( LUT(ndx,:).*(1-w) +...
        LUT(ndx+1,:).*w ).*(1-t) + ...
      ( LUT(ndx+nrows,:).*(1-w) + ...
        LUT(ndx+(nrows+1),:).*w ).*t).*(1-s) +...
     (( LUT(ndx+nw,:).*(1-w) +...
        LUT(ndx+1+nw,:).*w ).*(1-t) + ...
      ( LUT(ndx+nrows+nw,:).*(1-w) + ...
        LUT(ndx+(nrows+1+nw),:).*w ).*t).*s;
f = Ytest-F;

% if visual
%     figure(68)
%     ylim([0 0.05])
%     plot(Ytest,'r')
%     hold on
%     plot(F,'b')
%     pause(0.01)
%     plot(F,'g','linewidth',1)
%     
%     
% 
% end
% 
% if visual2
%     
%     figure(69)
%     hold on
%     plot3(x0(1),x0(2),x0(3),'r*')
%     pause(0.01)
%     plot3(x0(1),x0(2),x0(3),'*')
%     view(3)
% end

f = f(:);
% imagesc(F)