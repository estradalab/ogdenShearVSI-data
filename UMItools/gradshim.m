function [outmatrix,rotmatrix,aux]=gradshim(varargin)
%[outmatrix,rotmatrix,aux]=gradshim(inmatrix [,options])
%checks the coarse grained local derivatives to get the tilt
%rejects wrap around
%
%'diffcut', value		%eliminates pixels with large phase jumps, in degrees
%'dirichlet'            %adds zero gradient at image edge, useless
%'mask', maskmatrix     %masks by amplitude> 10% ; optional maskmatrix defines valid data for image
%'gradmask'             %mask out pixels where the intensity gradient is large (high susceptibility curvature)
%'planefit'             %defaults to planar fit to the derivatives

inmatrix=varargin{1};
outmatrix=inmatrix;

bowlfitflag=true; %set to false to suppress 
if any(strcmp(varargin, 'planefit'));
    bowlfitflag=false;
end


%size of the input matrix
siM=size(inmatrix); sigma=1; %don't fuck with this: siM=size(inmatrix); sigma=1;

%option 1
%diff angle cutuff, in degrees
diffanglecutoff=180;
if any(strcmp(varargin,'diffcut'));
    ind=find(strcmp(varargin,'diffcut'));
    diffanglecutoff=varargin{ind+1};
end

%option 2
%mask out pixels where amplitude gradient is large
gradmask=false(size(inmatrix));
if any(strcmp(varargin,'gradmask'));
    [smoothed_density,dummy,grad2]=blur(abs(inmatrix),2);
    gradmask(grad2>0.5*max(grad2(:)))=true;
end

%option 3
%mask out where there is no tissue, either by amplitude, or by using a mask
ind=find(strcmp(varargin,'mask'));
if ~isempty(ind) && nargin<ind && size(varargin{ind+1}) == size(inmatrix);
    maskmat=varargin{ind+1};
    inmatrix(maskmat)=NaN;
else
    inmatrix(isnan(inmatrix))=0;
    magn=abs(blur(inmatrix));        %blur
    magn=magn(:);
    magn=magn/max(magn);
    limvec=(0:0.05:1);
    ct=histc(magn(:),limvec);
    ct(1)=0; %reject empty space and 5% above it
    ct(2)=0; %reject 5% to 10%
    
    %Now ct contains a histogram of voxel intensities, where the first 10% have
    %been set to zero to reject empty space.
    
    %Now find those parts of the distribution where a threshold percentage
    %of the signal resides;
    thr=0.95;
    cumct=cumsum(ct)/sum(ct);                 %cumulative sum
    maxind=find(cumct<thr,1,'last');
    %display([maxind limvec(maxind), size(magn)]);
    inmatrix(magn(:)<0.25*limvec(maxind))=NaN;
end

%Now mask out where the density gradient was large
inmatrix(gradmask(:))=NaN;
aux.mask=~isnan(inmatrix);


%sparseindices=false(siM);
%samplevectory=ceil(sigma/2):sigma:siM(1);
%samplevectorx=ceil(sigma/2):sigma:siM(2);
%sparseindices(samplevectory,samplevectorx)=true;

%sparseM=inmatrix(sparseindices);
%h=numel(samplevectory);
%w=numel(samplevectorx);
%sparseM=reshape(sparseM,h,w);
sparseM=inmatrix;

%get the y slope, and the x slope account for coarse graining sigma
diffy=diff(angle(sparseM),1,1);
diffx=diff(angle(sparseM),1,2);

%generate masks to eliminate points with phase jumps
masky=abs(diffy)<diffanglecutoff/180*pi;
maskx=abs(diffx)<diffanglecutoff/180*pi;
siy=size(diffy);
six=size(diffx);
masky(:,1:six(2))= masky(:,1:six(2)) &  maskx(1:siy(1),:);
maskx(1:siy(1),:)= maskx(1:siy(1),:) &  masky(:,1:six(2));

%mask out bad data points;
diffy(~masky)=NaN;
diffx(~maskx)=NaN;


%get the x slope, account for coarse graining sigma


%option: impose zero gradient at the image edge
if any(strcmp(varargin,'dirichlet'))
    diffy(:,1)=0;
    diffy(:,end)=0;
    diffy(1,:)=0;
    diffy(end,:)=0;
    diffx(:,1)=0;
    diffx(:,end)=0;
    diffx(1,:)=0;
    diffx(end,:)=0;
end

    

%fit planes to the first derivative, in order to get the coefficients
[dummy1, gxx, gxy,offsx]=planefit(diffx);
[dummy2, gyx, gyy,offsy]=planefit(diffy);
[dummy3, diffxparameters]=bowlfit(diffx);
[dummy4, diffyparameters]=bowlfit(diffy);

%generate the output fit
si=size(inmatrix);
xaxis=(0:si(2)-1)-floor(si(2)/2)+(1-mod(si(2),2))/2;
yaxis=(0:si(1)-1)-floor(si(1)/2)+(1-mod(si(2),2))/2;
[X,Y]=meshgrid(xaxis,yaxis);
X2=X.^2;
Y2=Y.^2;
XY=X.*Y;
X2Y=X2.*Y;
XY2=X.*Y2;
Y3=Y2.*Y;
X3=X2.*X;

%generate planes tilted in the x and in the y direction, to correct for
%read out tilt
aux.ytilt=exp(-1i*(offsy*yaxis')*ones(1,numel(xaxis)));
aux.xtilt=exp(-1i*ones(numel(yaxis),1)*(offsx*xaxis));
dummyy=inmatrix.*aux.ytilt; thetay=mean(dummyy(~isnan(dummyy(:))));
dummyx=inmatrix.*aux.xtilt; thetax=mean(dummyx(~isnan(dummyx(:))));

%The field map can be fitted to second or to third order
%second order fits the derivatives to planes --> planefit
%third order fits the derivatives to quadratic bowls --> bowlfit

if ~bowlfitflag;
    %leading order only!
    A = offsx/sigma;
    B = offsy/sigma;
    C = gxx/sigma^2;
    D = gyy/sigma^2;
    E = (gxy +gyx)/2/sigma^2;
    F=0; G=0; H=0; J=0;
else
    
    %coefficients returned by bowlfit are in this order [offset, x, y, x2, xy, y2]
    %first order
    A=diffxparameters(1) /sigma;     %offset x          fx
    B=diffyparameters(1) /sigma;     %offset y			fy
    %second order
    C=diffxparameters(2) /sigma^2;  %gradient xx		fxx
    D=diffyparameters(3) /sigma^2;  %gradient yy		fyy
    E1=diffxparameters(3) /sigma^2; %gradient xy		fxy
    E2=diffyparameters(2)/sigma^2;  %gradient yx		fxy
    E=(E1+E2)/2;							%mean cross gradient
    %third order
    F1=diffxparameters(5)/sigma^3;   %2*xy coefficient in dx  : fxxy
    F2=2*diffyparameters(4)/sigma^3;%x2 coefficient in dy     : fxxy
    F=(F1+F2)/2;
    G1=2*diffxparameters(6)/sigma^3;    %Y2 coefficient in dx    :fyyx
    G2=diffyparameters(5)/sigma^3;       %xy coefficient in dy    :fyyx
    G=(G1+G2)/2;
    H=2*diffxparameters(4)/sigma^3;                %X2 coefficient in dx
    J=2*diffyparameters(6)/sigma^3;                 %Y2 coefficient in dy
end

ZZ = A*X+B*Y+(C/2)*X2+(D/2)*Y2+E*XY+...
    F/2*X2Y + G/2*XY2 + H/6*X3 + J/6*Y3;

aux.bowl=ZZ;
aux.coeff = [A B C D E F G H J];

rotmatrix=exp(-1i*ZZ);
testmatrix=inmatrix.*rotmatrix;
testmatrix(testmatrix==0)=NaN;
anglebins=(-pi:2*pi/365:pi)';
angledistribution=histc(angle(testmatrix(:)),anglebins);
[dummy,index]=max(angledistribution);
phi0=anglebins(index);

rotmatrix=rotmatrix*exp(-1i*phi0);

outmatrix=outmatrix.*rotmatrix;


%Now the matrix has been rotated to near zero
%make it exact

aux.phi0=phi0;
