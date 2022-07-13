function [out, kfilter]=lig_blur3d(varargin)
%[out, kfilter]=lig_blur3D(in [,kfilter,options]);
%in is a structure with fields image, kspace, pars as returned by varianms
%out is a structure with the same fields, but, the image has been blurred in 3 directions 
%options:
%'vox', [sz_ro sz_pe sz_pe2]    three component size vector, inmm, defining the virtual voxel
%'abs'                          blur absolute values, not complex numbers
%'grid', [nro npe nz]           three component vector giving the size of the output image
%the ligament acquisition multi-slice rather than full 3D !

%an if optional kfilter matrix has the right size and is nonzero, then the
%input kfilter is used; otherwise kfilter is calculated. The idea is that 
%if you want to blur many images, you can reuse the kfilter calculated the
%first time around.


in=varargin{1};

if ~isstruct(in);
	structflag=false;
	if ndims(varargin{1})==3;
		Min=varargin{1};
	else
		display('first argument must be 3d matrix or structure with field image containing the 3d matrix.')
	end
	
else
	structflag=true;
	Min=in.image;
	out=in;
end

%check if upsampling has been requested
if any(strcmp(varargin,'grid'));   %upsample the 3 d image
    ind=find(strcmp(varargin,'grid'));
    ftsize=varargin{ind+1};
else
    ftsize=size(Min);  %no zerofilling, keep the voxel number honest
end

Mout=zeros(ftsize);



%default blur coupling, applies if no vox is defined
sigmavec=0.8*[1 1 1];

if any(strcmp(varargin,'vox'));
    %virtual voxel size is given in mm
    ind=find(strcmp(varargin,'vox'));
    if nargin<ind+1 || ~isnumeric(varargin{ind+1});
        display('vox option requires a numeric vector input (3 elements, in mm) for the virtual voxel size.');
        return
    else
        si=size(Min);
        voxvec=varargin{ind+1};
        if numel(voxvec)==3;             %a voxel size has been handed down
            %sz is the voxel size handed down;
            sz(1)=in.pars.lro*10/si(1);     %FOV is in cm, convert to mm
            sz(2)=in.pars.lpe*10/si(2);
            sz(3)=in.pars.thk;              %slice thickness is in mm thickness 
            sigmavec=voxvec./sz*0.8/1.7;
        end
    end
end


ftest=exp(-(-10:10).^2/(2*sigmavec(3)^2));
ftest=ftest/sum(ftest);
display(['z-weighting of parent slice = ' num2str(round(ftest(11)*100)/100) ' in lig_blur3d.m...']);


if any(strcmp(varargin,'abs'));
    Min=abs(Min);
end



simage=size(Min);
Min=reshape(Min, [simage(1:3) prod(simage(4:end))]);

infilterflag=false;

if nargin>1 && isnumeric(varargin{2});
	if numel(varargin{2})==1;
        sigmavec=varargin{2}*[1 1 1];
    end
else
	sigma=0.8;
end


%if required, (re-)make the kspace filter
if ~infilterflag;
    siMin=size(Min);                    %size of the input image matrix
    sfilter=zeros(siMin(1:3));               %Fourier space filter
    [h w d]=size(sfilter);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
    fw=ceil(4*sigmavec);     %filter width
    x=(-fw(1):fw(1)); 
    y=(-fw(2):fw(2));
    z=(-fw(3):fw(3));
    fmask=zeros([numel(x),numel(y),numel(z)]);
    for ii=1:length(x);
        for jj=1:length(y);
            for kk=1:length(z);
                fmask(ii,jj,kk)=exp(-x(ii)^2/(2*sigmavec(1)^2) -y(jj)^2/(2*sigmavec(2)^2) -z(kk)^2/(2*sigmavec(3)^2));
            end
        end
    end
    fmask=fmask/sum(fmask(:));
    
    %fmask is a patch that will go into the middle of sfilter;
    sfilter(h/2+1+x,w/2+1+y,d/2+1+z)=fmask;
    kfilter=conj(fftn(sfilter));
end

if ndims(Min)==4;
    nec=siMin(4);
else
    nec=1;
end

%%
%ftsize=[256 48 20];
ftsize=max([ftsize(1:3); siMin(1:3)]);

for ne=1:nec;
    %works, doesn't upsample
    %FFT_Min_padded=fftn(Min(:,:,:,ne));
    %clear Mout; Mout(:,:,:,ne)=fftshift(ifftn(FFT_Min_padded.*kfilter)); figure; imagesc(angle(Mout(:,:,end/2)));
    
    filled_ft=zeros(ftsize(1:3));
    centerpoint=ftsize/2+1;
    
    
    FFT_Min_padded=(fftn(Min(:,:,:,ne)));
    convol=fftshift(FFT_Min_padded.*kfilter);     %convolution kernel
    sicon=size(convol);
    
    indvec1=(-sicon(1)/2:sicon(1)/2-1) + centerpoint(1);
    indvec2=(-sicon(2)/2:sicon(2)/2-1) + centerpoint(2);
    indvec3=(-sicon(3)/2:sicon(3)/2-1) + centerpoint(3);
    
    filled_ft(indvec1,indvec2,indvec3)=convol;
    
    
    Mout(:,:,:,ne)=fftshift(ifftn(fftshift(filled_ft)));
    
   
    %Mout(:,:,:,ne)=ifftshift(ifftn(dummy,ftsize));
    
    
end
%%

if structflag;
	out.image=Mout;
else
	out=Mout;
end

