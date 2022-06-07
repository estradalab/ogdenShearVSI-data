function out=softshim(varargin)
%out=softshim(in)
%autoshimming of image in structure in (.image, .pars)
%options: 
% 'filter'   [,sigma]        : gaussian smoothing of filters the amplitude
% 'phfilter' [, phsigma]   : high pass filter of the phase, sigma is pixel

in=varargin{1};
out=in;                   %the image will be updated

%determine how many slices, how many echoes
if isfield(in,'pars');
    if isfield(in.pars,'ns'); nsl=in.pars.ns; end
    if isfield(in.pars,'slices'); nsl=in.pars.ns; end
    if isfield(in.pars,'ne'); nec=in.pars.ne; end
    if isfield(in.pars,'echoes'); nec=in.pars.ne; end
else
    display('No parameters in input structure')
    return
end

%check if it is multislice single echo, or single slice multi echo
si_image=size(in.image);
if ndims(in.image)==3;
	if nsl==1; in.image=reshape(in.image,[si_image(1) si_image(2) 1 si_image(3)]); end
	if nec==1; in.image=reshape(in.image,[si_image(1) si_image(2) si_image(3) 1]); end
end

%1. see if amplitude filter was specified
ind=strcmp(varargin,'filter');
if any(ind);
	filterflag=true;
	iii-find(ind);
	if nargin>iii && isreal(varargin{iii+1});
		sigma=varargin{iii+1};
	else
		sigma=0.8; %pixels
	end
else
	filterflag=false;
end


%2. see if high-pass phase filter was specified
ind=strcmp(varargin,'phfilter');
if any(ind);
	phfilterflag=true;
	iii=find(ind);
	if nargin>iii && isreal(varargin{iii+1});
		phsigma=varargin{iii+1};
	else
		phsigma=5; %pixels
	end
else
	phfilterflag=false;
end

	
	
%imagetimes=in.pars.te+[0:nec-1]*in.pars.te2-0.001;
%imagetimes=[1:nec]*in.pars.te2;

for ns=1:nsl;
	for ne=1:nec;
		out.image(:,:,ns,ne)=gradshim(in.image(:,:,ns,ne));
		
		if filterflag || phfilterflag;
			ampmap=abs(out.image(:,:,ns,ne))+eps;
			phasemap=out.image(:,:,ns,ne)./ampmap;
		end
		
		if filterflag; %filter the amplitude
			out.image(:,:,ns,ne)=blur(ampmap,sigma).*phasemap;
		end
		
		if phfilterflag; %highpassfilter of the phase
			out.image(:,:,ns,ne)=out.image(:,:,ns,ne)./blur(phasemap,phsigma);
		end
			
	end
end

function [outmatrix,rotmatrix]=gradshim(varargin)
%[outmatrix,rotmatrix]=gradshim(inmatrix [,options])
%checks the coarse grained local derivatives to get the tilt
%rejects wrap around

inmatrix=varargin{1};
outmatrix=inmatrix;
%size of the input matrix
siM=size(inmatrix);
sigma=1;


inmatrix(isnan(inmatrix))=0;
magn=abs(blur(inmatrix,4));        %blur
magn=magn(:);
magn=magn/max(magn);
limvec=(0:0.05:1);
ct=histc(magn(:),limvec);
ct(1)=0; %reject empty space and 5% above it
ct(2)=0; %reject 5% to 10%

%Now find those parts of the distribution where a threshold percentage
%of the signal resides;
thr=0.95;
cumct=cumsum(ct)/sum(ct);                 %cumulative sum
maxind=find(cumct<thr,1,'last');
inmatrix(magn<0.5*limvec(maxind))=NaN;

%sparseindices=false(siM);
%samplevectory=ceil(sigma/2):sigma:siM(1);
%samplevectorx=ceil(sigma/2):sigma:siM(2);
%sparseindices(samplevectory,samplevectorx)=true;

%sparseM=inmatrix(sparseindices);
%h=numel(samplevectory);
%w=numel(samplevectorx);
%sparseM=reshape(sparseM,h,w);
sparseM=inmatrix;

%get the y slope, account for coarse graining sigma
diffy=diff(angle(sparseM),1,1);
diffy(abs(diffy(:))>pi/sigma)=NaN;

%get the x slope, account for coarse graining sigma
diffx=diff(angle(sparseM),1,2);
diffx(abs(diffx(:))>pi/sigma)=NaN;

%fit planes to the first derivative, in order to get the coefficients
%[dummy, gxx, gxy,offsx]=planefit(diffx);
%[dummy, gyx, gyy,offsy]=planefit(diffy);
[dummy, diffxparameters]=bowlfit(diffx);
[dummy, diffyparameters]=bowlfit(diffy);

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


%coefficients returned by bowlfit are in this order [offset, x, y, x2, xy, y2]
%first order
A=diffxparameters(1) /sigma;     %offset x				fx
B=diffyparameters(1) /sigma;     %offset y				fy
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

ZZ = A*X+B*Y+(C/2)*X2+(D/2)*Y2+E*XY+...
	F/2*X2Y + G/2*XY2 + H/6*X3 + J/6*Y3;


rotmatrix=exp(-1i*ZZ);
cor0=angle(inmatrix.*rotmatrix);
cor180=angle(-inmatrix.*rotmatrix);

std0=std(cor0(~isnan(cor0(:))));
std180=std(cor180(~isnan(cor180(:))));

anglebins=pi*[-1:0.1:1];
if std0<std180;
	ct=histc(cor0(:),anglebins);
	[dummy,maxind]=max(ct);
	offsetangle=anglebins(maxind);
else
	ct=histc(cor180(:),anglebins);
	[dummy,maxind]=max(ct);
	offsetangle=anglebins(maxind)-pi;
end

rotmatrix=rotmatrix*exp(-1i*offsetangle);

if any(strcmp(varargin,'mask'));
	inmatrix(isnan(inmatrix))=0;
	rotmatrix(isnan(inmatrix))=0;
	outmatrix=inmatrix.*rotmatrix;
else
	outmatrix=outmatrix.*rotmatrix;
end

%Now the matrix has been rotated to near zero
%make it exact
phi0=mean(angle(outmatrix(~isnan(inmatrix(:)))));
outmatrix=outmatrix*exp(-1i*phi0);

%for the optional absolute value gradient output, generate the gradients, 
%calculate offset bowl
%coefficients returned by bowlfit are in this order [offset, x, y, x2, xy, y2]
offset=ones(numel(X),1);
ZZZ=[offset(:), X(:) Y(:) X2(:) XY(:) Y2(:)];
fittedbowlx=ZZZ*diffxparameters;
fittedbowlx=reshape(fittedbowlx,si);
fittedbowly=ZZZ*diffyparameters;
fittedbowly=reshape(fittedbowly,si);

grad=(fittedbowlx.^2+fittedbowly.^2).^0.5;  %phase gradient (phase change per pixel)




