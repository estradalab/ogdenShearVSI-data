function out=varianms(varargin)
%out=varianms(filename, [option1, option2....]);
%if filename is a file, we are anlyzing fid's
%if filename is a directory, then we are reading the images (fdf's)
%if filename is a 4 digit number, and 'fid' is not an argument, open the image directory
%created by vjnmr, and stack them into a file;
%options:
%'plot'     : plots the slices, additional option 'col' makes a colorplot
%'iso'      : (default as 11/2016) zerofill (if time domain) or interpolate (if images) to make isotropic pixels
%'1k'       : force 1k pixels in each direction
%'512'      : force 512 pixels in each direction '192', '256', '384'
%'zf'       : zero fills 1.5 times in each direction, when analyzing fids
%'zf2'      : zero fills 2 times, in each direction, when analyzing fids
%'zf3'      : zerofills in the third dimension, for 3D acquisitions only
%'fliplr'   : flips the images left-right
%'flipud'	: flips the images up-down
%'flipimage': flips RO
%'filter'   : running gaussian filter, sigma=0.8pixels
%'apod'     : apodization (smooth roll off), details specified by a vector following 'apod':
%             1. no num input: rolls of at 90% kmax in x and y
%             2. scalar input: roll off begins at fraction of kmax in RO direction only, range [0 1];
%             3. two elements: rolls of in the RO and PE direction, at the fraction [fro fpe]
%             4. three elements: rolls of in the RO and PE and PE2 direction for 3d acquisitions, at the fraction [fro fpe fpe2]
%'upsample' : double the matrix in each direction by interpolation (careful, may not work with complex numbers
%'gradlag'  : shift the timedomain data by gradlag (microseconds) to reduce readout roll
%'freq'     : shift the resonance frequency (shift in time domain)
%'norm'     : normalize by the slice thickness
%'pixvec'   : [ pix1 pix2 pix3], or [pix1] or [pix1 pix2] to define image size larger than the acquired data

%check if the filename is image or fid, directory or file
[inputtype,varargin{1}]=check_img_or_fid(varargin{1});

%make it iso by default
%if ~any(strcmp(varargin,'noiso'));
%    varargin=[varargin 'iso'];
%end

varargin=[varargin 'flipimage'];


if strcmp(inputtype{1}, 'timedomain');          %a filename or directory for the fid was handed down
    display('Reconstructing complex images from fid...')
    if strcmp(inputtype{2}, 'directory');
        varargin{1}=[varargin{1} '/fid'];
    end
    out=reconfids(varargin{:});
    if strcmp(out.acquisitiontype,'1D');
        return;
    else
        image=out.image;
        pars=out.pars;
    end
    
else   %an image file or directory is to be read
    display('Reading real images from Varian .fdf files.');
    if strcmp(inputtype{2}, 'directory');
        

        pars=getparameters([varargin{1} '/procpar']);
        
        d=dir(varargin{1});
        names={d.name};
        for jj=1:numel(d);
            nam=names{jj}; %dummy is a cell array, the first entry is the file name
            matchflag(jj)=~isempty(regexp(nam, '\w+\.fdf'));
        end
        numslabsorslices=sum(matchflag);
        nam=names(matchflag);
        for jjj=1:numslabsorslices;
            rdfdf=readfdf([varargin{1} '/' nam{jjj}]);
            if jjj==1;
                image=zeros([size(rdfdf.image) numslabsorslices]);
            end
            si=size(rdfdf.image);
            switch numel(si);
                case 2
                    image(:,:,jjj)=rdfdf.image;     %multislice
                case 3
                    image(:,:,:,jjj)=rdfdf.image;   %multislab
            end
        end
        
        if strcmp(pars.seqfil,'epip');
            image=reshape(image,[si(1) si(2) pars.images pars.ns pars.ne]);
        end
        
        %pars.lpe=pars.span{1};
        %pars.lro=pars.span{2};
    elseif strcmp(inputtype{2}, 'file');
        rdfdf=readfdf(varargin{1});
        image=rdfdf.image;
        %pars=rdfdf.pars;
        %check this later, to see if axes are changed by orientation info
        pars=getparameters('procpar');
        %pars.lpe=pars.span{1};
        %pars.lro=pars.span{2};
    end
    
    %image files have been read, but the pixels may not be square
    si_image=size(image);
    imres(2)=pars.lpe./si_image(2);  %cm/pixel  xaxis
    imres(1)=pars.lro./si_image(1);  %cm/pixel  yaxis
    
    if any(strcmp('iso',varargin));
        %
        highest_resolution=min(imres);
        highest_resolution=highest_resolution(1);
        largest_FOV=max(pars.lpe,pars.lro);
        largest_FOV=largest_FOV(1);
        smallest_FOV=min(pars.lpe,pars.lro);
        
        %display a warning that this data set has anisotropic FOV
        if abs((1-smallest_FOV/largest_FOV))>0.01;
            display('This data set was acquired with a non-square field of view.')
            display('Images are made isotropic and square in varianms.m');
            display('out.pars.lpe and out.pars.lro were updated.');
        end
        
        npts=round(largest_FOV/highest_resolution);
        y_tempaxis=(1:si_image(1)+2)*imres(1);
        y_tempaxis(1)=0;
        y_tempaxis(end)=100;                        %100cmFOV
        x_tempaxis=(1:si_image(2)+2)*imres(2);
        x_tempaxis(1)=0;
        x_tempaxis(end)=100;                             %100cmFOV
        [Xtmp,Ytmp]=meshgrid(x_tempaxis,y_tempaxis);
        
        %upsample the dimension with fewer points
        nsl=numel(image)/(si_image(1)*si_image(2));
        newimages=zeros([npts npts si_image(3:end)]);
        [XX,YY]=meshgrid((1:npts)/npts*largest_FOV,(1:npts)/npts*largest_FOV);
        pars.lpe=largest_FOV;
        pars.lro=largest_FOV;
        for ns=1:nsl;
            src=zeros(size(Xtmp));
            src(1:si_image(1),1:si_image(2))=image(:,:,ns);
            newimages(:,:,ns)=interp2(Xtmp,Ytmp,src,XX,YY);
        end
        image=newimages;
    end
    out.kspace=[];
    
    %reshape the image files;
    si_image=size(image);
    
    %images are stored by VNMRJ is format slice###image###echo##.fdf
    %where the image### dimesion could be either the 2D array or epip
    %number, eand echo number is for multi echo
    nslices=pars.ns;
    nechos=pars.ne;
    nimages=prod(si_image)/(prod(si_image(1:2))*nechos*nslices);
    
    image=reshape(image,[si_image(1:2) nechos nimages nslices]);
    
    %if strcmp(pars.seqfil,'epip');  %epi image acquisition, forget about array dimension (is used for segmented acquisitions)
    %    image=squeeze(reshape(image,[si_image(1:2) pars.images pars.ns 1 1 1]));
    %
    %elseif isfield(pars,'arraydim') && pars.arraydim~=pars.nv;
    %    %if isfield(pars,'images') && pars.images==1;
    %        image=squeeze(reshape(image,[si_image(1:2) pars.arraydim/pars.nv pars.ne pars.ns 1 1 1]));
    %    %end
    %else
    %    image=squeeze(reshape(image,[si_image(1:2) pars.images pars.ne pars.ns 1 1 1]));
    %end
    image=squeeze(permute(image,[1 2 5 3 4])); %read phase slice echo# (image# or array#)
end
%End distinction of whether data is derived from FID or from IMG source

%image phase shifted in the phase directions
% 'former shphase'
if strcmp(inputtype{1}, 'timedomain');
    %shift the phase according to the parameters ppe and ppe2
    si_image=size(image);
    
    pro=out.pars.pro;
    ppe=out.pars.ppe;
    ppe2=out.pars.ppe2;
    
    lro=out.pars.lro;
    lpe=out.pars.lpe;
    lpe2=out.pars.lpe2;
    
    
    %testing out, I don't think there is a need for this because the
    %freq encode works via setting frequencies, but am not certain
    %shiftindex_pro= -round(si_image(1)*pro/lro);
    %pindex_pro=1+mod(shiftindex_pro+(0:si_image(1)-1),si_image(1));
    %image=image(pindex_pro,:,:,:,:);
    %display('phaseshifting in 1st dimension...');
    %if any(strcmp(varargin,'shphase2'));
    if ppe~=0;
        shiftindex_ppe= -round(si_image(2)*ppe/lpe);
        pindex_ppe=1+mod(shiftindex_ppe+(0:si_image(2)-1),si_image(2));
        image=image(:,pindex_ppe,:,:,:);
        display('phaseshifting in 2nd dimension...');
    end
    
    if strcmp(out.pars.seqfil,'mge3d') || strcmp(out.pars.seqfil,'ge3d');
        %if any(strcmp(varargin,'shphase3'));
        if ppe2~=0;
            shiftindex_ppe2= -round(si_image(3)*ppe2/lpe2);
            pindex_ppe2=1+mod(shiftindex_ppe2+(0:si_image(3)-1),si_image(3));
            image=image(:,:,pindex_ppe2,:,:);
            display('phaseshifting in 3rd dimension...');
        end
    end
end

%check if image should be flipped;
fliplr_flag=strcmp(varargin,'fliplr');
flipud_flag=strcmp(varargin,'flipud');
if any(fliplr_flag | flipud_flag);		%some flipping has been requested, execute it in the order it was requested
    si_image=size(image);
    numsl=numel(image)/(si_image(1)*si_image(2));
    for ns=1:numsl;
        if any(fliplr_flag);
            image(:,:,ns)=fliplr(image(:,:,ns));
        end
        if any(flipud_flag);
            image(:,:,ns)=flipud(image(:,:,ns));
        end
    end
end

%check if image should be filtered
if any(strcmp(varargin,'filter'));
    si_image=size(image);
    numsl=numel(image)/(si_image(1)*si_image(2));
    kfilter=[];
    ind=find(strcmp(varargin,'filter'));
    if nargin>ind && isreal(varargin{ind+1});
        %a filter length (in pixels) has been defined;
        [dummy,kfilter]=blur(image(:,:,1),varargin{ind+1});
    end
    
    for ns=1:numsl;
        display(['Filtering image ' num2str(ns) ' in varianms.m....'])
        [image(:,:,ns),kfilter]=blur(image(:,:,ns),kfilter);
    end
    
    %redistribute the dimensions
    image=reshape(image,si_image);
    clear smi;
    
    %if real image files were read, then take the absolute value of the
    %image
    if strcmp(inputtype{1}, 'image');
        image=abs(image);
    end
end

%check if image should be upsampled, but only if .img files were read
if strcmp(inputtype{1},'image') && (any(strcmp(varargin,'upsample')) || ...
        any(strcmp(varargin,'1k')) || any(strcmp(varargin,'1K'))  || any(strcmp(varargin,'512')) || any(strcmp(varargin,'256')) ...
        || any(strcmp(varargin,'1024')) || any(strcmp(varargin,'384')) || any(strcmp(varargin,'192')) ...
        || any(strcmp(varargin,'128')) || any(strcmp(varargin,'96')));
    si_image=size(image);
    numsl=numel(image)/(si_image(1)*si_image(2));
    if any(strcmp(varargin,'1k')) || any(strcmp(varargin,'1K')) || any(strcmp(varargin,'1024'));
        interpfactor(1)=1024/si_image(1);
        interpfactor(2)=1024/si_image(2);
        newi=1024;
        newj=1024;
    elseif any(strcmp(varargin,'512'));
        interpfactor(1)=512/si_image(1);
        interpfactor(2)=512/si_image(2);
        newi=512;
        newj=512;
    elseif any(strcmp(varargin,'384'));
        interpfactor(1)=384/si_image(1);
        interpfactor(2)=384/si_image(2);
        newi=384;
        newj=384;
    elseif any(strcmp(varargin,'192'));
        interpfactor(1)=192/si_image(1);
        interpfactor(2)=192/si_image(2);
        newi=192;
        newj=192;
    elseif any(strcmp(varargin,'256'));
        interpfactor(1)=256/si_image(1);
        interpfactor(2)=256/si_image(2);
        newi=256;
        newj=256;
    elseif any(strcmp(varargin,'128'));
        interpfactor(1)=128/si_image(1);
        interpfactor(2)=128/si_image(2);
        newi=128;
        newj=128;
    elseif any(strcmp(varargin,'96'));
        interpfactor(1)=96/si_image(1);
        interpfactor(2)=96/si_image(2);
        newi=96;
        newj=96;
    else
        interpfactor=[2 2];
        newi=interpfactor(1)*si_image(1);
        newj=interpfactor(2)*si_image(2);
    end
    smi=zeros(newi, newj,numsl); %interpolated search matrices
    [XX,YY]=meshgrid((0:si_image(2)),(0:si_image(1))');
    for ns=1:numsl;
        %diplay according to options
        if any(strcmp(varargin,'1k')) || any(strcmp(varargin,'1K')) || any(strcmp(varargin,'1024'));
            display(['Setting image ' num2str(ns) ' to 1024 resolution in varianms....']);
        elseif any(strcmp(varargin,'512'));
            display(['Setting image ' num2str(ns) ' to 512 resolution in varianms....']);
        elseif any(strcmp(varargin,'256'));
            display(['Setting image ' num2str(ns) ' to 256 resolution in varianms....']);
        elseif any(strcmp(varargin,'384'));
            display(['Setting image ' num2str(ns) ' to 384 resolution in varianms....']);
        elseif any(strcmp(varargin,'192'));
            display(['Setting image ' num2str(ns) ' to 192 resolution in varianms....']);
        elseif any(strcmp(varargin,'128'));
            display(['Setting image ' num2str(ns) ' to 128 resolution in varianms....']);
        elseif any(strcmp(varargin,'96'));
            display(['Setting image ' num2str(ns) ' to 96 resolution in varianms....']);
        else
            display(['Upsampling image twofold ' num2str(ns) ' in varianms.m....'])
        end
        ZZ=zeros(size(XX));		%magnitude
        ZZ(2:si_image(1)+1,2:si_image(2)+1)=abs(image(:,:,ns));
        smi(:,:,ns)=interp2(XX,YY,ZZ,(1:newj)/interpfactor(2),(1:newi)'/interpfactor(1),'*cubic');
        
        ZZi=zeros(size(XX));	%phase
        smi_phase=ZZi;
        ZZi(2:si_image(1)+1,2:si_image(2)+1)=angle(image(:,:,ns));
        smi_phase=interp2(XX,YY,ZZi,(1:newj)/interpfactor(2),(1:newi)'/interpfactor(1),'linear');
        
        smi(:,:,ns)=smi(:,:,ns).*exp(1i*smi_phase);
    end
    
    %redistribute the dimensions
    si_smi=size(smi);
    image=reshape(smi,[si_smi(1:2) si_image(3:end)]);
    clear smi;
end

%optional flipfb front to back
flipfb_flag=any(strcmp(varargin,'flipfb'));
if flipfb_flag && ndims(image)>2;
    si=size(image);
    image=image(:,:,end:-1:1,:,:,:,:);
    %pars.pss=pars.pss(end:-1:1);
end


out.image=image;
out.pars=pars;

%make axes for plotting
si=size(image);
if out.pars.phi==0;
    %RO is ydirection, PE is x direction
    out.pars.yaxis= ((-si(1)/2:si(1)/2-1)/si(1)*out.pars.lro+out.pars.pro)*10;
    out.pars.xaxis= ((-si(2)/2:si(2)/2-1)/si(2)*out.pars.lpe+out.pars.ppe)*10;
elseif out.pars.phi == 90;
    out.pars.yaxis=-((1:si(1))/si(1)*out.pars.lpe+out.pars.ppe)*10;
    out.pars.xaxis= ((1:si(2))/si(2)*out.pars.lro+out.pars.pro)*10;
else
    out.pars.yaxis=-((1:si(1))/si(1)*out.pars.lro)*10;
    out.pars.xaxis= ((1:si(2))/si(2)*out.pars.lpe)*10;
end
if ~isempty(strfind(out.pars.seqfil,'ge3d'));
    %the acquisition is 3D, account for the phase shifts in the phase
    %encode directions
    out.pars.zaxis=((si(3)/2-0.5:-1:-si(3)/2+0.5)/si(3)*out.pars.lpe2+out.pars.ppe2)*10;
end
if ~isempty(strfind(out.pars.seqfil,'ms'));
    %the acquisition is 2D multislice, account for the phase shifts in the phase
    %encode directions
    out.pars.zaxis=out.pars.pss*10;
end

% transitional code, shift labeling to axis1,axis2 and axis3, coordinates
% to designate voxel center, and ascending order of axes
out.pars.axis1 = - ((-si(1)/2+0.5:si(1)/2-0.5)/si(1)*out.pars.lro - out.pars.pro) *10; %convert mm to cm
out.pars.axis2 = - ((-si(2)/2+0.5:si(2)/2-0.5)/si(2)*out.pars.lpe + out.pars.ppe) *10; %convert mm to cm
if ~isempty(strfind(out.pars.seqfil,'ge3d'));  %case 3-d acquisition
    
    %axis3 gets a minus sign, this may be due to the VNMRJ messing with
    %animal orientations.
    out.pars.axis3 = - ( (-si(3)/2+0.5:si(3)/2-0.5) /si(3) *out.pars.lpe2 - out.pars.ppe2)*10;
end
if ~isempty(strfind(out.pars.seqfil,'ms'))    %case multislice acquisition
    out.pars.axis3 = out.pars.pss*10;
    
    %if there is more than one slice; flip the 3rd axis to descending order
    if numel(out.pars.axis3)>1
    [out.pars.axis3, sortindex]=sort(out.pars.axis3,2,'descend');
    out.image=out.image(:,:,sortindex,:,:,:,:,:);
    end
    %out.kspace=out.kspace(:,:,sortindex,:,:,:,:,:);
end

%make a surveyplot
if any(strcmp(varargin,'plot'));
    h=slice_surveyplot([{image} varargin]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=reconfids(varargin)
filename=varargin{1};
path=fileparts(filename);   %path of fid file to be read
if ~isempty(path);          %pwd is higher up
    procparfile=[path '/procpar'];
    pathname=path;
    %define the path to vnmrsys, to be able to find other auxiliary files
    currentdir=pwd;
    cd(path);
    ind=regexp(pwd,'vnmrsys');
    dummy=pwd;
    pvnmr=[dummy(1:ind-1) 'vnmrsys/'];
    cd(currentdir);
else
    procparfile='procpar';  %the current working directory contains the fid and procpar
    ind=regexp(pwd,'vnmrsys');
    dummy=pwd;
    pathname=pwd;
    pvnmr=[dummy(1:ind-1) 'vnmrsys/'];
end

%load a list of parameters from procpar, which may/will come in handy when
%generating images
pc=getparameters(procparfile);
pc.pathname=pathname;

%read the data into a matrix, td rows and however many columns
[acquireddata,header]=readit(filename);

%check if the gradient lag in the readout direction is specified
remainshift=0;
kshift=0;
if any(strcmp(varargin,'gradlag'));
    dii=find(strcmp(varargin,'gradlag'));
    if numel(varargin)>dii
        if isfloat(varargin{dii+1});
            graddelay=varargin{dii+1};
        else
            graddelay=30; %microseconds
        end
    else
        graddelay=30; %microseconds
    end
    si=size(acquireddata);
    ntd=si(1);
    dwell=pc.at/ntd;
    
    tdshift=round(graddelay*1e-6/dwell);
    remainshift=graddelay-tdshift*dwell*1e6;
    kshift=remainshift*1e-6/dwell*2*pi;
    
    dummy=zeros(size(acquireddata));
    if tdshift>=0
        dummy(1:ntd-tdshift,:)=acquireddata(1+tdshift:ntd,:);
    else
        tdshift=-tdshift;
        dummy(1+tdshift:ntd,:)=acquireddata(1:ntd-tdshift,:);
    end
    acquireddata=dummy;
    clear dummy;
end

%optional automatic frequencyshift;
if any(strcmp(varargin,'autof'));
    %get the offset frequency from and spuls experiment
    d=dir;
    spulsflag=false(numel(d),1);
    for jj=1:numel(d);
        spulsflag(jj)= ~isempty(strfind(d(jj).name,'spuls'));
    end
    if ~any(spulsflag);
        display('No spuls data in this directory');
        return
    end
    
    f=find(spulsflag,1,'last');
    spulsdata=readit([d(f).name '/fid']);
    spulspar=getparameters([d(f).name '/procpar']);
    spulsdwell=spulspar.at/numel(spulsdata);
    faxis=1/spulsdwell*[numel(spulsdata):-2:-numel(spulsdata)+2]/numel(spulsdata);
    spulsft=fftshift(fft(spulsdata));
    [maxsig,maxind]=max(fftshift(fft(spulsdata)));
    freqoffset=faxis(maxind);
end


%rearrange the data according to what pulse sequence it is
si=size(acquireddata); %td x number of traces
td=si(1);
switchstring=pc.seqfil;

%replace the switch string with a generic switchstring, to allow multiple
%versions of the same sequence (differing by a number or an appendix)
seqclass='dwssgg'; if ~isempty(findstr(switchstring,seqclass)); switchstring=seqclass; end
seqclass='gesege'; if ~isempty(findstr(switchstring,seqclass)); switchstring=seqclass; end
seqclass='fsems'; if ~isempty(findstr(switchstring,seqclass)) && ~strcmp('sems',switchstring); switchstring=seqclass; end
seqclass='stems'; if ~isempty(findstr(switchstring,seqclass)); switchstring=seqclass; end;
seqclass='tagcine'; if ~isempty(findstr(switchstring,seqclass)); switchstring=seqclass; end
seqclass='ge3d'; if ~isempty(findstr(switchstring,seqclass)); switchstring=seqclass; end
seqclass='spuls'; if ~isempty(findstr(switchstring,seqclass)); switchstring=seqclass; end
seqclass='profile1d'; if ~isempty(findstr(switchstring,seqclass)); switchstring=seqclass; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
acquisitiontype='2D';       %this is the default type
%any 3D sequences in the switch seegment below
%need to set this to '3D';
navigatorflag=false;        %default no navigator

switch lower(switchstring)
    case {'sems'; 'semsdw'; 'stems';'stemsdw';'gems';'gemsvc'; ...
            'sems_unequal_diff_grads'}; %%%%%%%%%%%%%%%%%%%%%% single echo/FID multislice
        rem=numel(acquireddata)/si(1)/pc.ns/pc.nv;
        
        %time phase slice block/array
        if strcmp(pc.seqcon,'nssnn') || strcmp(pc.seqcon,'ncsnn')
            acquireddata=reshape(acquireddata,[si(1) pc.ns*rem pc.nv]);
            acquireddata=permute(acquireddata,[1 3 2 4]);
        elseif (numel(pc.ti) == rem); %straight inversion recovery
            acquireddata=reshape(acquireddata,[si(1) pc.ns rem pc.nv]);  %seqcon=nccnn
            acquireddata=permute(acquireddata,[1 4 2 3]);
        elseif (numel(pc.te) == rem); %multi te
            acquireddata=reshape(acquireddata,[si(1) pc.ns pc.nv rem]);  %seqcon=nccnn
            acquireddata=permute(acquireddata,[1 3 2 4]);
        else
            acquireddata=reshape(acquireddata,[si(1) pc.ns pc.nv rem]);  %seqcon=nccnn
            acquireddata=permute(acquireddata,[1 3 2 4]);
        end
    case {'gems_blk'; 'gems_x'; 'gems_x3'}
        repeats=numel(acquireddata)/prod([si(1) pc.nsblock pc.nv pc.nblk]); % this is a movie acquisition with repeats
        acquireddata=reshape(acquireddata,[si(1) pc.nsblock pc.nv pc.nblk repeats]);
        acquireddata=permute(acquireddata,[1 3 2 4 5]);
        acquireddata=reshape(acquireddata,[si(1) pc.nv pc.nblk*pc.nsblock repeats]);
    case {'lig_stemems';'lig_stemems_nobutterfly'}
        rem=numel(acquireddata)/(si(1)*pc.nv*pc.ns);
        switch pc.tabscheme
            case 0
                acquireddata=reshape(acquireddata,[si(1),2,pc.ns,pc.nv/2,rem]);
                if pc.epiflag==1 %flip order of time domain acquisition
                    acquireddata(:,2,:,:,:,:)=acquireddata(end:-1:1,2,:,:,:,:);
                end
                acquireddata=permute(acquireddata,[1 2 4 3 5]);                 %stack PE dimension together
                acquireddata=reshape(acquireddata,[si(1),pc.nv,pc.ns,rem]);
            case 1
                acquireddata=reshape(acquireddata,[si(1),2,pc.ns,pc.nv/2,rem]);
                if pc.epiflag==1  %flip order of time domain acquisition
                    acquireddata(:,2,:,:,:,:)=acquireddata(end:-1:1,2,:,:,:,:);
                end
                acquireddata=acquireddata(:,[2 1],:,:,:);                       %flip for reversed PE direction in the 2 echo epi acquisition
                acquireddata=permute(acquireddata,[1 2 4 3 5]);                 %stack PE dimension together
                acquireddata=reshape(acquireddata,[si(1),pc.nv,pc.ns,rem]);         
        end
    case {'lig_stemems0721';'lig_stemems0722';'lig_stemems_dev';'lig_stemems1001b'}
         rem=numel(acquireddata)/(si(1)*pc.nv*pc.ns);
         acquireddata=reshape(acquireddata,[si(1),2,pc.ns,pc.nv/2,rem]);
         acquireddata=permute(acquireddata,[1 2 4 3 5]);                 %stack PE dimension together
         acquireddata=reshape(acquireddata,[si(1),pc.nv,pc.ns,rem]);
         
         if numel(pc.peshiftflag)>1
             acquireddata=reshape(acquireddata,[si(1),pc.nv,pc.ns,numel(pc.dro),numel(pc.peshiftflag)]);
             tempacq=squeeze(acquireddata(:,:,:,1,:));
             tempacq(:,2:(end),:,:)=tempacq(:,2:(end),:,:)+squeeze(acquireddata(:,1:(end-1),:,2,:));
             acquireddata=squeeze(tempacq);
             si=size(acquireddata);
         end
         
         
        
    case {'mems'; 'mgems';'mems_us';'mgems_us';'dwsemgems_us'; ...
            'dwgs1';'dwssgg'};%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% multiecho multislice (T2-weighting)
        rem=numel(acquireddata)/si(1)/pc.ne/pc.ns/pc.nv; %remainder
        
        %check if the data came in blocks, if it did, that is an arrayed
        %experiment, and varian stacks this into the second dimension!
        
        %check if there are multiple parameters arrayed
        if ~isempty(strfind(pc.array,',')); %there is a comma (punctuation) in the parameter pc.arrary
            number_of_arrayed_variables = numel(strfind(pc.array,','))+1;
            separatorindices=strfind(pc.array,',');
            switch number_of_arrayed_variables;
                case 2;
                    arrvar{1}=pc.array(1:(separatorindices(1)-1));
                    arrvar{2}=pc.array(separatorindices(1)+1:end);
                case 3
                    arrvar{1}=pc.array(1:(separatorindices(1)-1));
                    arrvar{2}=pc.array(separatorindices(1)+1:separatorindices(2)-1);
                    arrvar{3}=pc.array(separatorindices(2)+1:end);
            end
        else
            if ~isempty(pc.array);
                arrvar{1}=pc.array;
            end
        end
        
        if ~isempty(pc.array);
            for jj=1:numel(arrvar);
                remcount(jj)= numel(eval(['pc.' arrvar{jj}]));
            end
            rem=prod(remcount);
            acquireddata=reshape(acquireddata,[si(1) pc.ne pc.ns rem pc.nv]);  %seqcon = ccsnn (echo slice ph1 ph2 ph3)
            acquireddata=permute(acquireddata,[1 5 3 2 4]);
            si=size(acquireddata);
            acquireddata=reshape(acquireddata,[si(1:end-1) remcount(end:-1:1)]);
        else
            rem=1;
            remcount=1;
            acquireddata=reshape(acquireddata,[si(1) pc.ne pc.ns pc.nv rem]);  %seqcon = ccsnn (echo slice ph1 ph2 ph3)
            acquireddata=permute(acquireddata,[1 4 3 2 5]);
        end
        
        %time phase slice echo block/array
        
    case {'fsems'; 'fsemsdw'}  %%%%%%%%%%%%%%%%%%%%%%%%%%% fast spin echo (multiecho, different k each echo) multislice
        if strcmp(pc.navigator,'y');
            navigatorflag=true;
            %Navigator
            %for now it is not clear to me how procpar indexes the
            %navigators, but by inspection I see them at the end
            if strcmp(pc.nav_type,'pairwise');
                num_navs=2;
            else
                num_navs=1;
            end
            
            rem=numel(acquireddata)/si(1)/pc.ns/pc.nv/(pc.etl+num_navs)*pc.etl;
            nsegs=pc.nv/pc.etl;
            acquireddata=reshape(acquireddata,[si(1) (pc.etl + num_navs) pc.ns nsegs rem]);  %seqcon = nccnn
            
            siacq=size(acquireddata);
            nav_td_range=5;
            
            %assess where the navigators are finite
            testnav=acquireddata(:,(siacq(2)-(num_navs-1)):siacq(2),:,:,:,:,:);
            sinav=size(testnav);
            testnav=sum(abs(reshape(testnav,[sinav(1) prod(sinav(2:end))])),2);
            [maxnav,peakindex]=max(testnav);
            navindices=peakindex+[-1 0 1]';
            
            navigators=acquireddata(navindices,(siacq(2)-(num_navs-1)):siacq(2),:,:,:,:,:);
            navigator_segment_average=squeeze(mean(navigators,4));
            
            
            compNavi=squeeze(sum(navigators,1));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % what is indexed in acquired data:                  (td,  etl  ,ns,seg,array,array..)
            acquireddata=                             acquireddata(:,1:pc.etl,:,:,:,:,:);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %slice_specific_cn1=squeeze(sum(compNavi(1,:,:,:,:),3));
            %slice_specific_cn2=squeeze(sum(compNavi(2,:,:,:,:),3));
            
            if any(strcmp(varargin,'nav')); %navigator correction is called
                siacq=size(acquireddata);
                for sl=1:pc.ns;
                    for nsg=1:nsegs;
                        for rm=1:rem;
                            
                            [nav1.loc,nav1.amp]=qint(abs(navigators(:,1,sl,nsg,rm)));
                            [nav2.loc,nav2.amp]=qint(abs(navigators(:,2,sl,nsg,rm)));
                            temp1=prod(navigators(:,1,sl,nsg,rm));
                            temp2=prod(navigators(:,2,sl,nsg,rm));
                            nav1.phasefactor=temp1/abs(temp1);
                            nav2.phasefactor=temp2./abs(temp2);
                            
                            NN(nsg,sl,rm,1)=nav1.amp*nav1.phasefactor;
                            NN(nsg,sl,rm,2)=nav2.amp*nav2.phasefactor;
                            %%%%%%%%%
                            phasecorrection_type='complex_amplitude';
                            switch phasecorrection_type;
                                case 'phase';
                                    acquireddata(:,1:2:pc.etl,sl,nsg,rm)=acquireddata(:,1:2:pc.etl,sl,nsg,rm)/nav1.phasefactor; %odd echoes
                                    acquireddata(:,2:2:pc.etl,sl,nsg,rm)=acquireddata(:,2:2:pc.etl,sl,nsg,rm)/nav2.phasefactor; %even  echoes
                                case 'complex_amplitude';
                                    acquireddata(:,1:2:pc.etl,sl,nsg,rm)=acquireddata(:,1:2:pc.etl,sl,nsg,rm)/nav1.phasefactor/nav1.amp; %odd echoes
                                    acquireddata(:,2:2:pc.etl,sl,nsg,rm)=acquireddata(:,2:2:pc.etl,sl,nsg,rm)/nav2.phasefactor/nav2.amp; %even  echoes
                                case 'amplitude'
                                    acquireddata(:,1:2:pc.etl,sl,nsg,rm)=acquireddata(:,1:2:pc.etl,sl,nsg,rm)/nav1.amp; %odd echoes
                                    acquireddata(:,2:2:pc.etl,sl,nsg,rm)=acquireddata(:,2:2:pc.etl,sl,nsg,rm)/nav2.amp; %even  echoes
                            end
                            %%%%%%%%%_
                        end
                    end
                end
            end
            %put the phase dimensions next to each other
            acquireddata=permute(acquireddata, [1 2 4 3 5 6] );
            
            %stack the phase dimension
            si=[size(acquireddata) 1 1 1];
            acquireddata=reshape(acquireddata, [si(1) si(2)*si(3) si(4:end)]);
            
            out.navigators=squeeze(compNavi);
        else
            %No navigator
            rem=numel(acquireddata)/si(1)/pc.ns/pc.nv;
            acquireddata=reshape(acquireddata,[si(1) pc.etl pc.ns pc.nv/pc.etl rem]);  %seqcon = nccnn
            
            %put the phase dimensions next to each other
            acquireddata=permute(acquireddata, [1 2 4 3 5 6] );
            %stack the phase dimension
            si=[size(acquireddata) 1 1 1];
            acquireddata=reshape(acquireddata, [si(1) si(2)*si(3) si(4:end)]);
        end
        
        
        %Now reorder the k space data
        if isfield(pc,'pelist');
            if pc.nv~=numel(pc.pelist);
                pelistflag=false;
            else
                pelistflag=true;
            end
        else
            pelistflag=false;
        end
        if ~pelistflag;
            pelistfromfile=(importdata([pvnmr 'tablib/' pc.petable], '\t',1));
            pc.pelist=(pelistfromfile.data)';
            pc.pelist=pc.pelist(:);
        else
            %display('FSE order in .pars.pelist');
        end
        
        T2weightflag=false;
        if T2weightflag;
            %undo T2weighting of the echoes
            Tweight=100;
            si=size(pelistfromfile.data);
            taxis=pc.te+(0:si(2)-1)*pc.esp;
            wtmap=(ones(si(1),1)*exp(taxis/Tweight))';
            wtlist=wtmap(:);
            wtlist=wtlist/mean(wtlist);
            dummy=ones(size(acquireddata));
            for ii=1:numel(wtlist);
                dummy(:,ii,:,:,:)=dummy(:,ii,:,:,:)*wtlist(ii);
            end
            acquireddata=acquireddata.*dummy;
        end
        
        acquireddata(:,pc.pelist+pc.nv/2,:,:,:,:)=acquireddata;
        
    case {'gemsll'}; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% look locker
        acquireddata=reshape(acquireddata,[si(1) pc.etl pc.ns numel(pc.ti) pc.nv/pc.etl]);  %seqcon = nccnn
        %put the phase dimensions next to each other
        acquireddata=permute(acquireddata, [1 2 5 3 4]);
        %stack the phase dimension
        si=[size(acquireddata) 1 1 1];
        acquireddata=reshape(acquireddata, [si(1) si(2)*si(3) si(4:end)]);
        %Now reorder the k space data
        if isfield(pc,'pelist');
            if pc.nv~=numel(pc.pelist);
                pelistflag=false;
            else
                pelistflag=true;
            end
        else
            pelistflag=false;
        end
        if ~pelistflag;fsemsdwvc
            dummy=(importdata([pvnmr 'tablib/' pc.petable], '\t',1));
            pc.pelist=(dummy.data)';
            pc.pelist=pc.pelist(:);
        else
            %display('FSE order in .pars.pelist');
        end
        acquireddata(:,pc.pelist+pc.nv/2+1,:,:,:,:)=acquireddata;
        acquireddata=squeeze(acquireddata);
    case {'tagcine'; 'tagcinevc'}
        rem=numel(acquireddata)/si(1)/pc.ne/pc.nv;
        acquireddata=reshape(acquireddata,[si(1) pc.ne pc.nv rem 1 1]);
        acquireddata=permute(acquireddata, [1 3 2 4 5]);
    case {'ge3d'}
        acquisitiontype='3D';
        if numel(pc.te)==1;
            acquireddata=squeeze(reshape(acquireddata,[pc.np/2, pc.ne, pc.nv, pc.nv2]));
            
            if ndims(acquireddata)>3; %multidim or multiecho acquisition
                dd=ndims(acquireddata);
                acquireddata=permute(acquireddata, [1 3:dd 2]);
            end
        else
            acquireddata=squeeze(reshape(acquireddata,[pc.np/2, pc.nv, numel(pc.te), pc.nv2]));
            acquireddata=permute(acquireddata,[1 2 4 3]);
        end

    case {'spuls'}
        acquisitiontype='1D';
        %get the offset frequency from and spuls experiment
        spulsdata=acquireddata;
        spulspar=pc;
        spulsdwell=spulspar.at/numel(spulsdata);
        faxis=1/spulsdwell*[numel(spulsdata):-2:-numel(spulsdata)+2]/numel(spulsdata);
        spulsft=fftshift(fft(spulsdata));
        
        out.acquisitiontype=acquisitiontype;
        out.frequencyaxis=faxis';
        out.ftdata=spulsft;
        return
    case {'profile1d'}
        acquisitiontype='1D';
        profile=fftshift(fft(acquireddata),1); %fft along the first dimesnion only
        si=size(profile);
        out.acquisitiontype=acquisitiontype;
        
        roaxis=pc.pro*10 + pc.lro*10*(-si(1)/2:si(1)/2-1)/si(1);
        
        out.roaxis=roaxis(:);
        out.profile1D=profile;
        out.pars=pc;
        return
        
end
out.acquisitiontype=acquisitiontype;

%test: reverse phase domain data; result: don't do it
%si=size(acquireddata);
%acquireddata=acquireddata(:,si(2):-1:1,:,:,:);

%frequency shift if desired
%check if a frequency offset [,'freq', frequency] was handed down in varargin
if any(strcmp(varargin,'freq'));
    dii=find(strcmp(varargin,'freq'));
    if numel(varargin)>dii
        if isfloat(varargin{dii+1});
            freqoffset=varargin{dii+1};
        else
            %find the frequency offset using the spuls data, if it exists
            freqoffset=0; %microseconds
        end
    end
    si=size(acquireddata);
    ntd=si(1);
    dwell=pc.at/ntd;
    freqoffsetmatrix = exp(1i*freqoffset*(1:ntd)')*ones([1 si(2)]);
    restd=prod(si)/(si(1)*si(2));
    acquireddata=reshape(acquireddata, [si(1) si(2) restd]);
    for jj=1:restd;
        acquireddata(:,:,jj)=acquireddata(:,:,jj).*freqoffsetmatrix;
    end
    acquireddata=reshape(acquireddata,si);
end

%generate the images, with zerofill if requested
si=[size(acquireddata) 1];
nimages=prod(si(3:end));
%Default: no zero filling
mrows=td;
ncols=pc.nv;

%forced zero fills to fixed number of pixels if requested *****************
if any(strcmp(varargin,'128'));
    %this will be used if there is additional iso or zf argument in
    %varargin
    mrows=128;
    ncols=128;
end
if any(strcmp(varargin,'256'));
    %this will be used if there is additional iso or zf argument in
    %varargin
    mrows=256;
    ncols=256;
end
if any(strcmp(varargin,'384'));
    %this will be used if there is additional iso or zf argument in
    %varargin
    mrows=384;
    ncols=384;
end
if any(strcmp(varargin,'192'));
    %this will be used if there is additional iso or zf argument in
    %varargin
    mrows=192;
    ncols=192;
end
if any(strcmp(varargin,'512'));
    %this will be used if there is additional iso or zf argument in
    %varargin
    mrows=512;
    ncols=512;
end
if any(strcmp(varargin,'1024'));
    %this will be used if there is additional iso or zf argument in
    %varargin
    mrows=1024;
    ncols=1024;
end

%PIXVEC 
si=[size(acquireddata) 1 1];
cutvec=si(1:3);
if any(strcmp(varargin, 'pixvec'));
    index=find(strcmp(varargin, 'pixvec'));
    if isnumeric(varargin{index+1});
        pixvec=varargin{index+1};
        mrows=si(1);    %default before re-assignment by pixvec
        ncols=si(2);    %default before re-assignment by pixvec
        for jj=1:numel(pixvec);
            switch jj;
                case 1
                    mrows=pixvec(1);
                case 2
                    ncols=pixvec(2);
            end
        end
    else
        display('pixvec must be numeric');
    end
    
    % either cut data. or zerofill
    si=size(acquireddata);
    for jj=1:numel(pixvec);
        
        if pixvec(jj)<si(jj);                   %cut out large k data
            offset=(si(jj)-pixvec(jj))/2;
            indvec=offset+(1:pixvec(jj));
            switch jj
            case 1
                acquireddata=acquireddata(indvec,:,:,:,:,:);
            case 2
                acquireddata=acquireddata(:,indvec,:,:,:,:);
            case 3
                acquireddata=acquireddata(:,:,indvec,:,:,:);
            end
        end
    end
end



%**************************************************************************

%zerofilling either by multiplier zf or zf2, or to generate isotropic pixels
if any(strcmp( 'zf',varargin)) || any(strcmp( 'zf2',varargin))  || any(strcmp( 'iso',varargin));  %zerofill called for
    
    %generate the images, with zerofill if requested
    si=[size(acquireddata) 1];
    imres(1)=pc.lro./si(1);  %cm/pixel  yaxis (readout axis)
    imres(2)=pc.lpe./si(2);  %cm/pixel  xaxis (phase encode axis)
    
    
    %1. the raw length in the time domain and pe direction is given by
    mrows=td;
    ncols=pc.nv;
    %make the pixels isotropic
    if any(strcmp( 'iso',varargin));
        %default iso, no filling factor but isotropic pixels
        fff=1;
        
        if imres(1)<imres(2);                   %highest resolution in the readout direction;
            ncols=2*round(pc.lpe/imres(1)/2);   %set to isotropic pixels, but no extra zerofill
            
            if any(strcmp(varargin,'256'));     %upsample to 256 requested
                mrows=256;
                ncols=2*round(pc.lpe/(pc.lro/256)/2);
            end
            if any(strcmp(varargin,'128'));     %upsample to 256 requested
                mrows=128;
                ncols=2*round(pc.lpe/(pc.lro/128)/2);
            end
            if any(strcmp(varargin,'384'));     %upsample to 384 requested
                mrows=384;
                ncols=2*round(pc.lpe/(pc.lro/384)/2);
            end
            if any(strcmp(varargin,'192'));     %upsample to 384 requested
                mrows=192;
                ncols=2*round(pc.lpe/(pc.lro/192)/2);
            end
            if any(strcmp(varargin,'512'));     %upsample to 512 requested
                mrows=512;
                ncols=2*round(pc.lpe/(pc.lro/512)/2);
            end
            if any(strcmp(varargin,'1024'));    %upsample to 1024 requested
                mrows=1024;
                ncols=2*round(pc.lpe/(pc.lro/1024)/2);
            end
        else
            mrows=2*round(pc.lro/imres(2)/2);   %set to isotropic pixels, but no extra zerofill
            
            if any(strcmp(varargin,'256'));     %upsample to 256 requested
                ncols=256;
                mrows=2*round(pc.lro/(pc.lpe/256)/2);
            end
            if any(strcmp(varargin,'128'));     %upsample to 128 requested
                ncols=128;
                mrows=2*round(pc.lro/(pc.lpe/128)/2);
            end
            if any(strcmp(varargin,'384'));     %upsample to 384 requested
                ncols=384;
                mrows=2*round(pc.lro/(pc.lpe/384)/2);
            end
            if any(strcmp(varargin,'192'));     %upsample to 192 requested
                ncols=192;
                mrows=2*round(pc.lro/(pc.lpe/192)/2);
            end
            if any(strcmp(varargin,'512'));     %upsample to 512 requested
                ncols=512;
                mrows=2*round(pc.lro/(pc.lpe/512)/2);
            end
            if any(strcmp(varargin,'1024'));    %upsample to 1024 requested
                ncols=1024;
                mrows=2*round(pc.lro/(pc.lpe/1024)/2);
            end
        end
    end
    
    %zerofill factor 1.5 in each direction
    if any(strcmp( 'zf',varargin));
        fff=1.5;
    end
    %zerofill by factor 2 in each direction
    if any(strcmp( 'zf2',varargin));
        fff=2;
    end
    mrows=2*round(mrows*fff/2);
    ncols=2*round(ncols*fff/2);
end



if any(strcmp(varargin,'apod')) && strcmp(acquisitiontype,'2D');
    %apodization filter is requested
    display('varianms: apodization option was set');
    
    %1. check if a numeric argument follows apod
    apodindex=find(strcmp(varargin,'apod'));
    if nargin==apodindex || ~isnumeric(varargin{apodindex+1}); %'apod' is last option or next option not numeric
        si=size(acquireddata);
        k2=si(1:2)/2;       %number of (ky, kx) /2
        kc=0.9;             %cutoff k
        w=1-kc;             %width of the apod roll off
        
        ky=[-k2(1):k2(1)-1]'./k2(1);
        kx=[-k2(2):k2(2)-1]./k2(2);
        kabs=sqrt((ky.^2)* ones(1,si(2)) + ones(si(1),1)*(kx.^2) );
        ktophat=kabs<kc;
        kring =kabs>kc+w;
        W=cos(pi*(kabs-kc)/(2*w)).^2;
        W(ktophat(:))=1;
        W(kring(:))=0;
    elseif nargin>apodindex && isnumeric(varargin{apodindex+1})&& numel(varargin{apodindex+1})==1
        %'apod' is followed by a single number, apodize time domain only
        si=size(acquireddata);
        %kc=0.9;             %cutoff k
        kc=varargin{apodindex+1};
        k2=si(1:2)/2;       %number of (ky, kx) /2
        
        w=1-kc;             %width of the apod roll off
        
        ky=[-k2(1):k2(1)-1]'./k2(1);
        kabs=abs(ky)* ones(1,si(2));
        ktophat=kabs<kc;
        kring =kabs>kc+w;
        W=cos(pi*(kabs-kc)/(2*w)).^2;
        W(ktophat(:))=1;
        W(kring(:))=0;
    elseif nargin>apodindex && isnumeric(varargin{apodindex+1}) && numel(varargin{apodindex+1})==4
        %selective set a patch in kspace to zero, useful for killing
        %stimulated echo signatures
        si=size(acquireddata);
        W=ones(si(1:2));
        zpl=varargin{apodindex+1};
        W(zpl(1):zpl(2),zpl(3):zpl(4))=0;
    elseif nargin>apodindex && isnumeric(varargin{apodindex+1}) && numel(varargin{apodindex+1})==2
        %define a central oval in kspace and kill the rest
        si=size(acquireddata);
        k2=si(1:2)/2;                   %number of (ky, kx) /2
        
        if all(varargin{apodindex+1} > 1) %apod cutoff given in terms of number of k-space steps away from k=0;
            kcut=varargin{apodindex+1};   
        elseif all(varargin{apodindex+1} <= 1)
            kcut=k2.*varargin{apodindex+1}; %apod cutoff given as a fraction of kspace
        end
        
            
       % w=0.1;                          %width of the apod roll off
        
        ky=(-k2(1):k2(1)-1)'/kcut(1);    %normalized ky
        kx=(-k2(2):k2(2)-1)./kcut(2);    %normalized kx
        %kabs=sqrt((ky.^2)* ones(1,si(2)) + ones(si(1),1)*(kx.^2) );
        ktophat=abs(ones(si(1),1)*kx)<1 & abs(ky*ones(1,si(2)))<1;
        W=zeros(si(1:2));
        W(ktophat(:))=1;
    end
    
    if ndims(W)==2;
        for jj=1:nimages;
            acquireddata(:,:,jj)=acquireddata(:,:,jj).*W;
        end
    end
end


                 
               

%initialize the phase warp matrix to account for an off-center echo
YY=exp(1i*(-mrows/2:mrows/2-1)'*ones(1,ncols)/mrows*kshift);

%Fourier transform with zerofilling;
if strcmp(acquisitiontype,'2D');
    %reorder the slices to be ascending
    [sliceposition,orderindex]=sort(pc.pss);
    if any(orderindex ~= (1:numel(orderindex)));
        display('varianms.m ordered interleaved slices...');
        display('varianms.m updated the slice positions in out.pars.pss ...');
        display('varianms.m ordered interleaved k-space slices ...');
        pc.pss=pc.pss(orderindex);
        acquireddata=acquireddata(:,:,orderindex,:,:,:);
    end
    
    
    siacq=size(acquireddata);
    ioff=mrows/2-siacq(1)/2;
    joff=ncols/2-siacq(2)/2;
    pic=zeros([mrows ncols si(3:end)]);
    for jj=1:nimages;
        paddedslice=zeros(mrows,ncols);
        paddedslice(ioff+(1:siacq(1)), joff+(1:siacq(2)))=acquireddata(:,:,jj);
        %experimental: virtual reduction of bandwidth by a factor 2
        if any(strcmp('BW',varargin));
            windowsize=2;
            paddedslice=filter(ones(1,windowsize)/windowsize,1,paddedslice);
        end
        pic(:,:,jj)=fftshift(fft2(fftshift(conj(paddedslice))));
        
        if any(strcmp(varargin,'flipimage'));
            pic(:,:,jj)=pic(end:-1:1,:,jj);
        end
        
        if kshift~=0;
            pic(:,:,jj)=pic(:,:,jj).*YY;
        end
    end
    
    
    %reshape the images into the original size
    si=[size(acquireddata) 1 1 1];
    sipic=size(pic);
    longsizevector=[sipic(1:2) si(3:end)];  %pads the size vector with a couple of ones
    
    pic=reshape(pic, longsizevector);
    
    % if the 2D acquisition is a contiguous stack of slices and a 3d pixvec
    % was specified, then one could up the resolution in the stack
    % direction by zerofilling
    if exist('pixvec','var');
        if numel(pixvec) == 3
            
            newpic=zeros([pixvec longsizevector(4:end)]);
            higherDimNumber=prod(longsizevector(4:end));
            
            for jj=1:higherDimNumber %loop over dims 4 and up
                
                dummy=fft(pic(:,:,:,jj),[],3); %zerofill in the third dimentsion
                
                newpic(:,:,:,jj)=ifft(dummy,pixvec(3),3);
            end
            pic=reshape(newpic,[pixvec longsizevector(4:end)]);
            
            % update pc.pss
            zax=pc.pss;
            z1_lim=zax(1)-sign(zax(2)-zax(1))*(pc.thk/10)/2;
            zend_lim=zax(end)+sign(zax(2)-zax(1))*(pc.thk/10)/2;
            deltaz=(zend_lim-z1_lim)/pixvec(3);
            pc.pss=z1_lim +(0.5:pixvec(3))*deltaz;
            
        end
    end
    
    
elseif strcmp(acquisitiontype,'3D');
    
    siacq=size(acquireddata);
    ioff=mrows/2-siacq(1)/2;
    joff=ncols/2-siacq(2)/2;
    
    siacq=[size(acquireddata) 1 1 1];
    
    %zerofill in third dimansion if 'zf3' passed down%%%
    if any(strcmp(varargin, 'zf3'));
        nz=1.5*siacq(3);
        nzoff=round(siacq(3)/4);
        pic=zeros([mrows ncols nz siacq(4:end)]);
        display('Zero filling 3rd dimension.');
    else
        if exist('pixvec','var');
            mrows=max(pixvec(1),siacq(1));        %the requested output size
            ncols=max(pixvec(2),siacq(2));
            
            ioff=max(0,round((pixvec(1)-siacq(1))/2));
            joff=max(0,round((pixvec(2)-siacq(2))/2));
            
            if(pixvec(3)>siacq(3));
                nzoff=max(0,round(pixvec(3)-siacq(3))/2);
                nz=pixvec(3);
            else
                nz=siacq(3);
                nzoff=0;
                display('requested pixvec(3) is less than data acquisition size(3)');
            end
        else
            nz=siacq(3);
            nzoff=0;
        end
        
        pic=zeros([mrows ncols nz siacq(4:end)]);
    end
    
    
    si=size(acquireddata);
    KT=true(si(1:3));
    if any(strcmp(varargin,'apod'));
        %apodization filter is requested
        
        
        %1. check if a numeric argument follows apod
        apodindex=find(strcmp(varargin,'apod'));
        
        if nargin>apodindex && isnumeric(varargin{apodindex+1}) && numel(varargin{apodindex+1})<=3;
            display('varianms: apodization option was set on a 3D data set');
            
            %define a central oval in kspace and kill the rest
            si=size(acquireddata);
            k2=si(1:3)/2;                           %number of (ky, kx) /2
            kcut=[varargin{apodindex+1} 1 1 1];     %cutoff fraction ky kx kz
            kcut=kcut(1:3);
           
            
            ky=[-k2(1):k2(1)-1]'./k2(1)/kcut(1);    %normalized ky
            kx=[-k2(2):k2(2)-1]./k2(2)/kcut(2);     %normalized kx
            kz=[-k2(3):k2(3)-1]./k2(3)/kcut(3);     %normalized kz
            %kabs=sqrt((ky.^2)* ones(1,si(2)) + ones(si(1),1)*(kx.^2) );
            ktophat=abs(ones(si(1),1)*kx)<1 & abs(ky*ones(1,si(2)))<1;
            
            KT=repmat(ktophat,[ 1 1 numel(kz)]);
            blankedkz=abs(kz)>1;
            KT(:,:,blankedkz)=false;                %KT is a 3D mask
        end
    end
    
    
    
    
    for jj=1:siacq(4);
        paddedblock=zeros(mrows,ncols,nz);
        paddedblock(ioff+(1:siacq(1)), joff+(1:siacq(2)),nzoff+(1:siacq(3)))=1*KT.*squeeze(acquireddata(:,:,:,jj));
        %pic(:,:,:,jj)=fftshift(fftn(fftshift(conj(paddedblock))));
        
        pic(:,:,:,jj)=fftshift(fftn(fftshift(conj(paddedblock))));
        if any(strcmp(varargin,'flipimage'));
            pic(:,:,:,jj)=pic(end:-1:1,:,:,jj);
        end
    end
    
    %reshape the images into the original size%%%%%%%%%%%
    si=[size(acquireddata) 1 1 1];
    sipic=size(pic);
    pic=reshape(pic, [sipic(1:2) nz si(4:end)]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end




%normalize the pixel intensity to the voxel size if requested
if any(strcmp(varargin,'norm'));
    display('normalized images by slice thickness');
    pic=pic/pc.thk;
end

out.image=pic;
out.pars=pc;
out.kspace=acquireddata;  %end reconfids

function [acquireddata,header]=readit(filename)
%open the fid
file_id=fopen(filename, 'r', 'b');          %always big endian format, even for Linux

%Start reading the file header to get the first 6 elements
[hdr, count] = fread(file_id, 6, 'int32'); %count is number of

header.nblocks = hdr(1);        %data blocks in a file (second dimension, if arrayed acquisition)
header.ntraces = hdr(2);        %number of traces
header.complex_td= hdr(3)/2;    %real and imaginary points (first dimension = time);
ebytes = hdr(4);                % 4 bytes per point	(double precision real, double precision imag)
tbytes = hdr(5);                %np*ebytes
bbytes = hdr(6);                %ntraces*tbytes+nhead*sizeof(datablockhead=28bytes)-total

%Read 2 16bit integers to see the data type and the status of the binary
[data_holder, count] = fread(file_id, 2, 'int16');
vers_file_id = data_holder(1);
status = data_holder(2);            %status determines whether data are store as 16 or 32 bit integers or floats
BinStat = fliplr(dec2bin(status));  %determine the binary format from the status info
%xx00xxxx means 16 bit data
%xx10xxxx means 32 bit data
%else means float

[data_holder, count] = fread(file_id, 1, 'int32');
nhead = data_holder(1);

data_holder = zeros(hdr(3)*hdr(2), 1);

%INT16 data type **********************************************************
if ((BinStat(3) == '0') & (BinStat(4) == '0'))
    for nblocks_count=1:header.nblocks
        for nhead_count=1:nhead
            fread(file_id, 14, 'int16');			%read in the block headers (28 bytes)
        end
        [b, count] = fread(file_id, hdr(2)*hdr(3), 'int16');   	%read in the actual data (ntraces*np)
        %put the data in the vector
        data_holder( ((nblocks_count-1)*hdr(2)*hdr(3)+1) : nblocks_count*hdr(2)*hdr(3) ) = b;
    end
    %INT32 data type **********************************************************
elseif ((BinStat(3) == '1') & (BinStat(4) == '0'))
    for nblocks_count=1:header.nblocks
        for nhead_count=1:nhead
            [d,count]=fread(file_id, 14, 'int16');			%read in the block headers (28 bytes)
            bhead=d;
        end
        [b, count] = fread(file_id, hdr(2)*hdr(3), 'int32');   	%read in the actual data (ntraces*np)
        %put the data in the vector
        data_holder( ((nblocks_count-1)*hdr(2)*hdr(3)+1) : nblocks_counthdr(2)*hdr(3) ) = b;
    end
    %FLOAT data type **********************************************************
else
    for nblocks_count=1:header.nblocks
        for nhead_count=1:nhead
            fread(file_id, 14, 'int16');			%read in the block headers (28 bytes)
        end
        [b, count] = fread(file_id, hdr(2)*hdr(3), 'float');   	%read in the actual data (ntraces*np)
        %put the data in the vector
        data_holder( ((nblocks_count-1)*hdr(2)*hdr(3)+1) : nblocks_count*hdr(2)*hdr(3) ) = b;
    end
end

%real imaginary pairs, create complex array
acquireddata=complex(data_holder(1:2:end),data_holder(2:2:end));
acquireddata=reshape(acquireddata,[header.complex_td numel(acquireddata)/header.complex_td]);
function [p,y,a] = qint(yin)
%QINT - quadratic interpolation of three adjacent samples
%
% [p,y,a] = qint(yin)
% yin is length 3
%
% returns the extremum location p, height y, and half-curvature a
% of a parabolic fit through three points.
% Parabola is given by y(x) = a*(x-p)^2+b,
% where y(-1)=ym1, y(0)=y0, y(1)=yp1.
%p = (yp1 - ym1)/(2*(2*y0 - yp1 - ym1));
%y = y0 - 0.25*(ym1-yp1)*p;
%a = 0.5*(ym1 - 2*y0 + yp1);

p = (yin(3) - yin(1))/(2*(2*yin(2) - yin(3) - yin(1)));
y = yin(2) - 0.25*(yin(1)-yin(3))*p;
a = 0.5*(yin(1) - 2*yin(2) + yin(3));
function [inputtype, fileidstring]=check_img_or_fid(fileidstring)
%check if the filename is image or fid, directory or file
if isnumeric(fileidstring); 
    fileidstring=['0000' num2str(fileidstring)]; 
    fileidstring=fileidstring((end-3):end);
end
[parent, fnamin, ext]=fileparts(fileidstring); %#ok<*ASGLU>
switch ext;
    case '.img';
        inputtype={'image'; 'directory'};
    case '.fdf';
        inputtype={'image'; 'file'};
    case '.fid'
        inputtype={'timedomain'; 'directory'};
    case ''
        if strcmp(fnamin, 'fid');
            inputtype={'timedomain'; 'file'};
        else
            if (numel(fnamin)==4 && ~isempty(str2num(fnamin)));
                d=dir;
                istime=false(numel(d),1);
                isimg=false(numel(d),1);
                for jj=1:numel(d);
                    istime(jj)=~isempty(strfind(d(jj).name,fnamin));
                    isimg(jj)=~isempty(strfind(d(jj).name,'.img'));
                end
                isfid=~isimg;
                fid_directory_index=find(istime&isfid); inputtype={'timedomain'; 'directory'};
                %image_directory_ndex=find(istime&isimg); inputtype={'image'; 'directory'};
                
                if ~isempty(fid_directory_index);
                    fileidstring=d(fid_directory_index).name; %#ok<FNDSB>
                    %for epi images force reading image directory
                    if ~isempty(strfind(fileidstring,'epip')); %it is an epip so read the images from the image directory
                        image_directory_index=find(istime&isimg); inputtype={'image'; 'directory'};
                        fileidstring=d(image_directory_index).name; %#ok<FNDSB>
                    end
                else
                    err=MException('directory_not_found:wrong_time_stamp',['No directories with time stamp ' fnamin]);
                    throw(err);
                end
                
            end
        end
end
function h=slice_surveyplot(in)
pic=in{1};
if any(strcmp(in,'col'));
    cmap=jet(128);
else
    cmap=gray(128);
end
%auto choose the layout
si=size(pic);
maxamp=max(pic(:));
meanamp=mean(abs(pic(:)));
switch ndims(pic);
    case 2
        nfig=1;
        nx=1;
        ny=1;
        nimperfigure=1;
    case 3;						%stack of slices or time stamped single slice
        nfig=1;
        nx=ceil(sqrt(si(3)));
        ny=ceil(si(3)/nx);
        nimperfigure=si(3);	%number of images per figure
    case 4
        nfig=1;
        nx=si(3);
        ny=si(4);
        nimperfigure=nx*ny;
    case 5
        nfig=si(5);
        nx=si(3);
        ny=si(4);
        nimperfigure=nx*ny;
end

%plot it
for ff=1:nfig;
    %h(ff)=figure('Name', pc.comment,'NumberTitle', 'off');
    h(ff)=figure;
    colormap(cmap);
    for ii=1:nimperfigure;
        tsubplot(ny,nx,ii);
        imagesc(abs(pic(:,:,ii)));
        shading flat;
        axis image
        set(gca,'clim',[0 3*meanamp],'TickDir', 'in','Xticklabel','','YTicklabel','');
        simage=size(pic);
        ht=textn(0.9,0.1, num2str(ii+(ff-1)*nimperfigure));
        set(ht, 'color', [0 0.8 0]);
    end
end
function out=readfdf(filename)
%read a varian image, return it in a matrix, loosely
%based on ReadVarFDF by Lana Kaiser, Varian
if ~exist(filename, 'file')==2;
    display(['Image file ' filename ' does not exist!']);
    return;
else
    [fid,msg]=fopen(filename,'r');
end

%1. read text lines and define a pars structure
parscanflag=true;
jj=1;
rdln=(fgetl(fid)); %read the first line
while parscanflag;
    jj=jj+1;
    
    %read the line
    rdln=(fgetl(fid));
    
    %remove '*' and '[]' and '"' from the line
    rdln=strrep(rdln,'*', '');
    rdln=strrep(rdln,'[]','');
    rdln=strrep(rdln,'"', '''');
    
    [token,remain]=strtok(rdln);
    remain=strrep(remain,' ','');
    %display(['pars.' remain]);
    eval(['pars.' remain]);
    
    %display(strtok(rdln));
    if strfind(rdln, 'checksum');
        %display(['Last fdf header line at line ' num2str(jj)]);
        parscanflag=false;
    end
    
    if jj>50; parscanflag=false; end
end

%2D or 3D data?
dimensionality=pars.rank;

%bigendian or little endian?
if pars.bigendian==0;
    endian='l';
else
    endian='b';
end

%define the dimensions,read the data
if dimensionality ==2;
    dim(1)=pars.matrix{1};
    dim(2)=pars.matrix{2};
    status = fseek(fid, -dim(1)*dim(2)*pars.bits/8, 'eof');
    image=fread(fid,[dim(1), dim(2)],'float32',endian);
elseif dimensionality==3;
    dim(1)=pars.matrix{1};
    dim(2)=pars.matrix{2};
    dim(3)=pars.matrix{3};
    status = fseek(fid, -dim(1)*dim(2)*dim(3)*pars.bits/8, 'eof');
    image=fread(fid,[dim(1)*dim(2), dim(3)],'float32',endian);
    image=reshape(image,dim);
end

%permute the image data so that it is consistent with direct reconstruction
indvec=[1:ndims(image)];
indvec(1)=2;
indvec(2)=1;
%exchange x and y axes
image=permute(image,indvec);
%invert axis directions
si=(size(image));
image=image(si(1):-1:1,si(2):-1:1,:,:,:);
out.image=image;
out.pars=pars;

%Now determine the size of the data, locate the beginning of the binary
%data, read it in, redimension
fclose(fid);
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
function out=makemask(in)
out=in;
blurvoxmm=[0.4 0.4 0.4];
blurreddata=blur3d(in,'vox',blurvoxmm);
blurredpic=blurreddata.image;

[blurredsignallevel,blurrednoiselevel] = estimate_noiselevel(blurredpic); 

out.mask=abs(blurredpic)>2.5*blurrednoiselevel;

