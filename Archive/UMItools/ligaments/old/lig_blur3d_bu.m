function [out, kfilter]=lig_blur3d(varargin)
%[out, kfilter]=lig_blur3D(in [,kfilter,options]);
%in is a structure with fields image, kspace, pars as returned by varianms
%out is a structure with the same fields, but, the image has been blurred in 3 directions 
%options:
%'vox', [sz_ro sz_pe sz_pe2]    three component size vector, inmm, defining the virtual voxel
%'abs'                          blur absolute values, not complex numbers
%'zf', [nro npe nz]             three component vector giving the size of the output image
%the ligament acquisition multi-slice rather than full 3D !

%an if optional kfilter matrix has the right size and is nonzero, then the
%input kfilter is used; otherwise kfilter is calculated. The idea is that 
%if you want to blur many images, you can reuse the kfilter calculated the
%first time around.


in=varargin{1};

if ~isstruct(in);
    display('abort: First input to lig_blur3D.m must be a structure containing the data.')
    return
else
    Min=in.image;
    out=in;
end

%check if upsampling has been requested
if any(strcmp(varargin,'zf'));   %upsample the 3 d image
    ind=find(strcmp(varargin,'zf'));
    zfvec=varargin{ind+1};
else
    zfvec=[1 1 1];  %no zerofilling, keep the voxel number honest
end



%default blur coupling, applies if no vox is defined
sigmavec=0.7*[1 1 1];

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

if nargin>1 && isnumeric(varargin{2});
	kfilter=varargin{2};
	infilterflag=true; %default
	%1. check if the sizes of filter and inmatrix Min match
    infilterflag=false;
	if isempty(kfilter);
        sigmavec=[0.8 0.8 0.8];
    end
    if size(kfilter)==[1 1];
        sigma=kfilter(1);
    end

else
	infilterflag=false;
	sigma=0.8;
end


%if required, (re-)make the kspace filter
if ~infilterflag;
    siMin=size(Min);            %size of the input image matrix
    sfilter=zeros(siMin);       %Fourier space filter
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
    kfilter=conj(fftn(sfilter,siMin+4));
end

if ndims(Min)==4;
    nec=siMin(4);
else
    nec=1;
end


%Min_padded is 4 entries larger than Min in each dimension
Min_padded=cat(3,Min(:,:,1,:),Min(:,:,1,:),Min,Min(:,:,end,:),Min(:,:,end,:));
Min_padded=cat(2,Min_padded(:,1,:,:),Min_padded(:,1,:,:),Min_padded,Min_padded(:,end,:,:),Min_padded(:,end,:,:));
Min_padded=cat(1,Min_padded(1,:,:,:),Min_padded(1,:,:,:),Min_padded,Min_padded(end,:,:,:),Min_padded(end,:,:,:));



ftsize= size(Min_padded).*zfvec;
Mout_padded=zeros(ftsize);


for ne=1:nec;
    FFT_Min_padded=fftn(Min_padded(:,:,:,ne));
    Mout_padded(:,:,:,ne)=fftshift(ifftn(FFT_Min_padded.*kfilter,ftsize));
    Mout(:,:,:,ne)=Mout_padded(1+zfvec(1)*2:(end-zfvec(1)*2),1+zfvec(2)*2:(end-zfvec(2)*2),1+zfvec(3)*2:(end-zfvec(3)*2),ne);
end

out.image=Mout;

