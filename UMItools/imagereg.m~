function out=imagereg(varargin)
% out=imagereg(in,options)
% in is a an array of structures, returned by multiple invocations of varianms
% out is the same structure array, but now with co-registered images
% options:
% 'matchradius', fraction, where the fraction defines what part of the image is registered
%
% This routine handles 4 and higher dimensional images by preforming
% the registration transformation on the 1st element of dimensions 4 and up
% The underlying assumption is multi-echo or arraryed acquisitions in which
% the first acquisition has the largest intensity (first echo, no
% diffusion, etc.)

% defaults:
matchradius_value=0.8;
inputlist={'matchradius'};
parseinput(inputlist,varargin);

% Input checking:
% loop to check what type of stuff has been passed down
for jj=1:nargin;
    structflag(jj)  =isstruct(varargin{jj});
    cellflag(jj)    =iscell(varargin{jj});
    realflag(jj)    =isreal(varargin{jj});
    numflag(jj)     =isnumeric(varargin{jj});
    %optional other types
end

% check if a previously defined x and y ROI has been passed down
% defefault: full range
in=varargin{1};
out=in;
si=[size(in(1).image) 1];
xwinvec=1:si(2);
ywinvec=1:si(1);
if any(cellflag);
    dummy=find(cellflag);
    nolimitsflag = true;
    for jj=1:numel(dummy);
        testcell=varargin{dummy(jj)};
        if strcmp(testcell{1},'ROI limits') && nolimitsflag;
            xlim=testcell{2};
            ylim=testcell{3};
            xwinvec=testcell{4};
            ywinvec=testcell{5};
            nolimitsflag = false;
        end
    end
end

%define the size of the image arrays
nux=numel(xwinvec);
nuy=numel(ywinvec);
nuz=si(3);
nud= numel(in);
preimage=zeros([nuy nux nuz nud]);
pre_fft=preimage;


ycoordinates=(((1:nuy)-nuy/2)/(nuy/2))'*ones(1,nux);
xcoordinates=ones(nuy,1)*(((1:nux)-nux/2)/(nux/2));
mask=exp(-(xcoordinates.^2+ycoordinates.^2).^4/matchradius_value^8);

%
for nd=1:nud;
    for nz=1:nuz;
        preimage(:,:,nz,nd)=abs(in(nd).image(ywinvec,xwinvec,nz));
        prefft(:,:,nz,nd)=fft2(preimage(:,:,nz,nd).*mask);
    end
end
postimage=preimage;

% correlate
for nd=2:nud;
    for nz=1:nuz;
        test=ifft2(prefft(:,:,nz,1).*conj(prefft(:,:,nz,nd))./abs(prefft(:,:,nz,1).*conj(prefft(:,:,nz,nd))));
        
        [dummy,maxind]=max(test(:));
        [yind,xind]=ind2sub([nuy nux],maxind);
        if yind>nuy/2; yind=yind-nuy; end;
        if xind>nux/2; xind=xind-nux; end;
        
        xshift_z(nz,nd)=xind-1;
        yshift_z(nz,nd)=yind-1;
    end
    
    %average the shifts of all slices, for this data set;
    xshift=round(mean(xshift_z));
    yshift=round(mean(yshift_z));
    xshiftvec=mod((1:nux)-1 - xshift(nd),nux)+1;
    yshiftvec=mod((1:nuy)-1 - yshift(nd),nuy)+1;
    
    out(nd).image=in(nd).image(yshiftvec,xshiftvec,:,:,:,:);
end


