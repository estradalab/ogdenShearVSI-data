function outmatrix = fit2Do3(varargin)
% outmatrix = fit2Do3(inmatrix [,'plot']);
% meaning fit a 2D taylor expansion, to order 3
% inmatrix is a 2D matrix with optional NaN entries, assumed to lie on
% a cartesian underlying grid.
% outmatrix is a third order Taylor expansion in two coordinates x and y,
% centered in the middle of the inmatrix
% the plot option overlays inmatrix and fit, and plots the residue as well

inmatrix=real(varargin{1});
si=size(inmatrix);
logicalindex=~isnan(inmatrix);

%rescale and shift the inmatrix to assure benign fitting
fit_offset=mean(inmatrix(logicalindex(:)));
s_inmatrix=inmatrix-fit_offset;                 %shifted inmatrix
span=std(s_inmatrix(logicalindex(:)));
ss_inmatrix=s_inmatrix/span;                  %shifted and scaled inmatrix

%these axes are not related to any physical axes, they serve to define the
%grid over which the matrix is fitted.
xaxis=(0:si(2)-1)-floor(si(2)/2)+(1-mod(si(2),2))/2;
yaxis=(0:si(1)-1)-floor(si(1)/2)+(1-mod(si(2),2))/2;
xaxis=xaxis/max(xaxis);
yaxis=yaxis/max(yaxis); 

[XX, YY]=meshgrid(xaxis,yaxis);
X2=XX.^2;
Y2=YY.^2;
XY=XX.*YY;
X2Y=X2.*YY;
XY2=XX.*Y2;
X3=X2.*XX;
Y3=Y2.*YY;


% calculate the basis functions (not orthogonal, problem?)
DM=ss_inmatrix(logicalindex);
offset=ones(sum(logicalindex(:)),1);
y=YY(logicalindex);
x=XX(logicalindex);
y2=Y2(logicalindex);
x2=X2(logicalindex);
xy=XY(logicalindex);
x2y=X2Y(logicalindex);
xy2=XY2(logicalindex);
x3=X3(logicalindex);
y3=Y3(logicalindex);


ZZ=[offset(:), x(:) y(:) x2(:) xy(:) y2(:) x2y(:) xy2(:) x3(:) y3(:)];
coefficients=ZZ\DM;

%calculate offset bowl
offset=ones(numel(XX),1);
ZZZ=[offset(:), XX(:) YY(:) X2(:) XY(:) Y2(:) X2Y(:) XY2(:) X3(:) Y3(:)];
outmatrix=ZZZ*coefficients;
outmatrix=reshape(outmatrix,si);

outmatrix=outmatrix*span + fit_offset;

if any(strcmp(varargin,'plot'));
    figure;
    mesh(inmatrix);
    hold on;
    mesh(outmatrix);
end

