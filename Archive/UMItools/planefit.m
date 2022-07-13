function [fittedplane, gx, gy,offs]=planefit(varargin)
% [fittedplane, gx,gy]=planefit(Plane, logicalindex);
Plane=varargin{1};
si=size(Plane);
if nargin==1
    logicalindex=~isnan(Plane);
else
    logicalindex=varargin{2};
	if ~islogical(logicalindex)
		logicalindex=logical(logicalindex);
	end
end


xaxis=(0:si(2)-1)-floor(si(2)/2)+(1-mod(si(2),2))/2;
yaxis=(0:si(1)-1)-floor(si(1)/2)+(1-mod(si(2),2))/2;

[XX, YY]=meshgrid(xaxis,yaxis);

DM=Plane(logicalindex);
offset=ones(sum(logicalindex(:)),1);
ctilt=YY(logicalindex);
dtilt=XX(logicalindex);
ZZ=[offset(:), ctilt(:) dtilt(:)];
stitchup=ZZ\DM;

gx=stitchup(3);
gy=stitchup(2);
offs=stitchup(1);
	
%calculate offset plane
offset=ones(numel(XX),1);
ZZZ=[offset(:), YY(:) XX(:)];
fittedplane=ZZZ*stitchup;
fittedplane=reshape(fittedplane,si);
