function [Mout, kfilter,grad2]=blur(varargin)
%[Mout, kfilter]=blur(Min [,kfilter]);
%blurs the image matrix Min using gaussian convolution, 2D fft;
%the second optional input is either a scalar defining the width of the
%gaussian nearest neighbour blurring, in pixels (default 0.8). 
%0 means no blurring, 
%an if optional kfilter matrix has the right size and is nonzero, then the
%input kfilter is used; otherwise kfilter is calculated. The idea is that 
%if you want to blur many images, you can reuse the kfilter calculated the
%first time around.
Min=varargin{1};

if ndims(Min)~=2;
    display('blur routine requires 2D-matrix, or scalar, input.')
    return
end

sigma=0.8; %default

if nargin==2;
	kfilter=varargin{2};
	infilterflag=true; %default
	%1. check if the sizes of filter and inmatrix Min match
	if any(size(Min)~=size(kfilter));
		infilterflag=false;
		if isempty(kfilter);
			sigma=0.8;
		elseif isnumeric(kfilter);
			sigma=kfilter(1);
		end
	else
		%sizes do match, but also check if there are nonzero entries in the infilter
		if ~any(kfilter(:)~=0);
			infilterflag=false;
			sigma=0.8;
		end
	end

else
	infilterflag=false;
	sigma=0.8;
end

%if required, (re-)make the kspace filter
if ~infilterflag;
    sfilter=zeros(size(Min));
    [h w]=size(sfilter);
	
    fw=ceil(3*sigma);
    x=(-fw:fw);
    y=x';
    fmask=zeros(numel(x));
    for jj=1:length(x);
        for ii=1:length(y);
            fmask(ii,jj)=exp(-(x(jj)^2+y(ii)^2)/(2*sigma^2));
        end
    end
    fmask=fmask/sum(fmask(:));
    
    %fmask is a patch that will go into the middle of sfilter;
    sfilter(floor(h/2+1+y),floor(w/2+1+x))=fmask;
    kfilter=conj(fftshift(fft2(fftshift(sfilter))));
end

FFT_Min=fftshift(fft2(fftshift(Min)));
Mout=fftshift(ifft2(fftshift(FFT_Min.*kfilter)));

if isreal(Min);
    Mout=real(Mout);
end
%Mout=fftshift(ifft2((FFT_Min.*kfilter)));

%to avoid confusion
%Mout=abs(Mout);

%calculate the gradient as well;
[fx,fy]=gradient(abs(Mout));
grad2=(fx.^2+fy.^2).^0.5;
%replace the frame of Mout with the input matrix
%framewidth 2 sigma
%Mout(1:2*sigma,:)=Min(1:2*sigma,:);
%Mout((end-2*sigma):end,:)=Min((end-2*sigma):end,:);
%Mout(:,1:2*sigma)=Min(:,1:2*sigma);
%Mout(:,(end-2*sigma):end)=Min(:,(end-2*sigma):end);
