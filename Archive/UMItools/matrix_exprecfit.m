function out=matrix_exprecfit(varargin)
%out=matrix_exprecfit(incube,taxis,maxamp[,errmaxamp])
%fit to exponential recovery
%incube: nro X npe (x numslices) X numel(taxis) is a stack of images acquired at different times
%the times passed down in the taxis.
%out is structure with fields 'amplitude' and 'tau', result of a pixelwise
%
%the exponential recovery  ' y= a*(1-exp(bt)) '
%can't be done by straight linear regression, so this routine implements an
%iterative linear regression, where an initail amplitude matrix a is
%guessed, then updated iteratively
%log(1-y/a) = log(c) + bt; 

incube=varargin{1};
taxis=varargin{2};


sincube=size(incube);
nofslices=sincube(3);
noftimes=sincube(4);
bincube=zeros(sincube);
noiseincube=zeros(sincube);





function [ampout, tau, R] = matrix_recfit(incube,ampin,taxis)
%this routine fits the log of the signal to a line
%log(c-y) = log(a) + bt;

noiselevel=1;

[tmax,tmaxind]=max(taxis);
%ensure the starting guess is different from the longest TR amplitude
if all(ampin==incube(:,:,tmaxind));
    ampin=incube(:,:,end)*1.1;
end

%size of the incoming data
sincube=size(incube);

%generate a cube of the same size, populated with the initial amplitude
ampincube=zeros(size(sincube));
for jj=1:sincube(3);
    ampincube(:,:,jj)=ampin;
end

oneovers2 = (abs(ampincube-incube)/noiselevel).^2; %scale noise to account for change of variables to log(y)
incube=log(abs(ampincube-incube));

x=ones(sincube);
for jj=1:sincube(3);
    x(:,:,jj)=x(:,:,jj)*taxis(jj);
end
y=incube;
xy=x.*y;
x2=x.*x;

Sx2=sum(x2.*oneovers2,3);
Sy=sum(y.*oneovers2,3);
Sx=sum(x.*oneovers2,3);
Sxy=sum(xy.*oneovers2,3);
Soneovers2=sum(oneovers2,3);

Delta=Soneovers2.*Sx2-Sx.^2;
b= (Soneovers2.*Sxy - Sx.*Sy)./Delta;
a= (Sx2.*Sy-Sx.*Sxy)./Delta;

ampout=exp(a);
tau=-1./b;
R=-b;

%errors in the fitted parameters:


