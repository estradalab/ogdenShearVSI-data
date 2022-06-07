function [out, kfilter]=blur3d(varargin)
%[out, kfilter]=lig_blur3D(in [,kfilter,options]);
% in is a structure with fields image, kspace, pars as returned by varianms
% OR a 3d data set, in which case the units are pixels
% the full width half max FWHM of the gaussian filter can be specified by
% the 'vox' option, either in mm (structure input), or in pixels matrix
% input.
% verified:  Gaussian with FWHM=1 has sigma=0.424
% 
% out is a structure with the same fields, but, the image has been blurred in 3 directions 
% options:
% 'vox', [sz_ro sz_pe sz_pe2]    three component size vector, inmm, defining the virtual voxel
% 'abs'                          blur absolute values, not complex numbers
% 'grid', [nro npe nz]           three component vector giving the size of the output image
% the ligament acquisition multi-slice rather than full 3D !

%an if optional kfilter matrix has the right size and is nonzero, then the
%input kfilter is used; otherwise kfilter is calculated. The idea is that 
%if you want to blur many images, you can reuse the kfilter calculated the
%first time around.


in=varargin{1};

if ~isstruct(in);
	structflag=false;
	if ndims(varargin{1})==3
		Min=varargin{1};
        input2dflag=false;
    elseif ismatrix(varargin{1});
        Min=repmat(varargin{1},1,1,8);
        input2dflag=true;
	end
	
else
	structflag=true;
	Min=in.image;
	out=in;
end

%check if upsampling has been requested
ftsize=size(Min);
if any(strcmp(varargin,'grid'))   %upsample the 3 d image
    ind=find(strcmp(varargin,'grid'));
    ftsize(1:numel(varargin{ind+1}))=varargin{ind+1};
end
Mout=zeros(ftsize);

if any(strcmp(varargin,'vox'))
    %virtual voxel size is given in mm
    ind=find(strcmp(varargin,'vox'));
    if nargin<ind+1 || ~isnumeric(varargin{ind+1});
        display('vox option requires a numeric vector input (3 elements, in mm) for the virtual voxel size.');
        return
    else
        si=size(Min);
        voxvec=varargin{ind+1};
        if numel(voxvec)==3
            
            voxvec(voxvec(:)==0)=0.001;

            
            if structflag  %a structure permitting calculation of voxel size has been handed down
                
                %sz calculates the true voxel size
                sz(1)=in.pars.lro*10/si(1);     %FOV is in cm, convert to mm
                sz(2)=in.pars.lpe*10/si(2);
                % sz(3) depends on 2d or 3D acquisition
                if strcmp(in.acquisitiontype,'2D');
                    sz(3)=in.pars.thk;              %slice thickness is in mm thickness
                elseif strcmp(in.acquisitiontype,'3D');
                    sz(3)=in.pars.lpe2*10/si(3);
                else
                    display('undefined acquisition type of the input structure in blur3d.')
                    return
                end
            else
                sz=ones(1,3);  
            end
            sigmavec=voxvec./sz*0.424;
        else
            sz=ones(1,3);
            sigmavec=ones(1,3)*voxvec(1)./sz*0.424+0.05;
        end
    end
else
    % default blur coupling: sigma=0.424 voxels, gives FWHM = voxel size
    % applies if no 'vox' is defined by the blur input
    sigmavec=0.424*[1 1 1];
end

if any(strcmp(varargin,'abs'))
    Min=abs(Min);
end

simage=size(Min);
Min=reshape(Min, [simage(1:3) prod(simage(4:end))]);

infilterflag=false;

if nargin>1 && isnumeric(varargin{2});
	if numel(varargin{2})==1;
        sigmavec=varargin{2}*[1 1 1];
    end
end


%if required, (re-)make the kspace filter
if ~infilterflag;
    siMin=size(Min);                    %size of the input image matrix
    sfilter=zeros(siMin(1:3));               %Fourier space filter
    [h w d]=size(sfilter);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
    fw=ceil(3*sigmavec);     %filter width
    x=(-fw(1):fw(1)); 
    y=(-fw(2):fw(2));
    z=(-fw(3):fw(3));
    fmask=zeros([numel(x),numel(y),numel(z)]);
    for ii=1:length(x)
        for jj=1:length(y)
            for kk=1:length(z)
                fmask(ii,jj,kk)=exp(-x(ii)^2/(2*sigmavec(1)^2) -y(jj)^2/(2*sigmavec(2)^2) -z(kk)^2/(2*sigmavec(3)^2));
            end
        end
    end
    fmask=fmask/sum(fmask(:));
    display(['z-weighting of parent pixel = ' num2str(round(max(fmask(:)*100))/100) ' in blur3d.m ...']);
    
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
    filled_ft=zeros(ftsize(1:3));
    centerpoint=ftsize/2+1;
    
    
    FFT_Min_padded=(fftn(Min(:,:,:,ne)));         %not zerofilled or padded!
    
    %multiplication of kspace filter and kspace image 
    convol=fftshift(FFT_Min_padded.*kfilter);     %convolution kernel
    sicon=size(convol);
    
    indvec1=(-sicon(1)/2:sicon(1)/2-1) + centerpoint(1);
    indvec2=(-sicon(2)/2:sicon(2)/2-1) + centerpoint(2);
    indvec3=(-sicon(3)/2:sicon(3)/2-1) + centerpoint(3);
    
    %this is the zerofilling step
    filled_ft(indvec1,indvec2,indvec3)=convol;
    
    
    Mout(:,:,:,ne)=fftshift(ifftn(fftshift(filled_ft)));
end

siMout=size(Mout);
siMout= [siMout(1:3) simage(4:end)];
Mout=reshape(Mout,siMout);

%%
if structflag;
	out.image=Mout;
else
	out=Mout;
    if input2dflag
        out=Mout(:,:,1);
    end
end

