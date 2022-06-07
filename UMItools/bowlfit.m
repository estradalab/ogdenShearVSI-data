function [fittedbowl, coefficients]=bowlfit(varargin)
% [fittedbowl, coefficients]=bowlfit(Plane, logicalindex);
%coefficients are for [offset, x, y, x2, xy, y2]
Plane=varargin{1};
si=size(Plane);
if nargin==1;
    logicalindex=~isnan(Plane);
else
    logicalindex=varargin{2};
end


xaxis=(0:si(2)-1)-floor(si(2)/2)+(1-mod(si(2),2))/2;
yaxis=(0:si(1)-1)-floor(si(1)/2)+(1-mod(si(2),2))/2;

[XX, YY]=meshgrid(xaxis,yaxis);
X2=XX.^2;
Y2=YY.^2;
XY=XX.*YY;

DM=Plane(logicalindex);
offset=ones(sum(logicalindex(:)),1);
y=YY(logicalindex);
x=XX(logicalindex);
y2=Y2(logicalindex);
x2=X2(logicalindex);
xy=XY(logicalindex);
ZZ=[offset(:), x(:) y(:) x2(:) xy(:) y2(:)];
coefficients=ZZ\DM;
	
%calculate offset bowl
offset=ones(numel(XX),1);
ZZZ=[offset(:), XX(:) YY(:) X2(:) XY(:) Y2(:)];
fittedbowl=ZZZ*coefficients;
fittedbowl=reshape(fittedbowl,si);
