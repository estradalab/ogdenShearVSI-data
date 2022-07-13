function out=varianms(varargin)
%out=varianms(filename, [option1, option2....]);
%if filename is a file, we are anlyzing fid's
%if filename is a directory, then we are reading the images (fdf's)
%created by vjnmr, and stack them into a file;
%options:
%'plot'     : plots the slices
%'iso'      : zerofill to make isotropic pixels
%'zf'       : zero fills 1.5 times in each direction, when analyzing fids
%'zf2'      : zero fills 2 times, in each direction, when analyzing fids
%'fliplr'   : flips the images left-right
%'flipud'	: flips the images up-down
%'filter'   : running gaussian filter, sigma=0.8pixels
%'upsample' : double the matrix in each direction by interpolation (careful, may not work with complex numbers
%'gradlag'  : shift the timedomain data by gradlag (microseconds) to reduce readout roll
%'freq'     : shift the resonance frequency (shift in time domain)
%'shphase'  : shift the image in the PE direction, [-1 1]

%check if the filename is image or fid, directory or file
inputtype=check_img_or_fid(varargin{1});

%make it iso
if any(strcmp(varargin,'1k')) && ~any(strcmp(varargin,'iso'))
    varargin=[varargin 'iso'];
end

if strcmp(inputtype{1}, 'timedomain');          %a filename or directory for the fid was handed down
    display('Reconstructing complex images from fid...')
    if strcmp(inputtype{2}, 'directory');
        varargin{1}=[varargin{1} '/fid'];
    end
    out=reconfids(varargin{:});
    image=out.image;
    pars=out.pars;
else   %an image file or directory is to be read
    display('Reading real images from Varian .fdf files.');
    if strcmp(inputtype{2}, 'directory');
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
        pars=getparameters([varargin{1} '/procpar']);
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
    imres(1)=pars.lpe./pars.nv;  %cm/pixel
    imres(2)=pars.lro./pars.fn;  %cm/pixel
    
    si_image=size(image);
    
    if any(strcmp('iso',varargin)) && si_image(1)~=si_image(2);
        %upsample the dimension with fewer points
        nsl=numel(image)/(si_image(1)*si_image(2));
        szm=max(si_image(1),si_image(2));
        newimages=zeros([szm szm si_image(3:end)]);
        [XX,YY]=meshgrid(1:szm,1:szm);
        for ns=1:nsl;
            xskip=ceil(szm/si_image(2));
            yskip=ceil(szm/si_image(1));
            newimages(yskip:end,xskip:end,ns)=...
                interp2((1:si_image(2))/si_image(2)*szm, ...
                (1:si_image(1))/si_image(1)*szm,...
                image(:,:,ns),...
                XX(yskip:end,xskip:end),YY(yskip:end,xskip:end));
        end
        image=newimages;
    end
    out.kspace=[];
    
    %reshape the image files;
    si_image=size(image);
    image=squeeze(reshape(image,[si_image(1:2) pars.ne pars.arraydim pars.ns]));
    image=squeeze(permute(image,[1 2 5 3 4 6 7]));
end

%check if image should be shifted in the phase direction
if any(strcmp(varargin,'shphase'));
    ii=find(strcmp(varargin,'shphase'));
    if nargin>ii && isnumeric(varargin{ii+1}) && varargin{ii+1}>=-0.5 && varargin{ii+1}<=0.5;
        
        si_image=size(image);
        shiftindex=-round(si_image(2)*varargin{ii+1});
        pindex=1+mod(shiftindex+(0:si_image(2)-1),si_image(2));
        numsl=numel(image)/(si_image(1)*si_image(2));
        for ns=1:numsl;
            image(:,:,ns)=image(:,pindex,ns);
            
        end
    else
        display('aborted: varianms Input argument shphase must be followed by a number in [-0.5 0.5].');
        return
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
    for ns=1:numsl;
        display(['Filtering image ' num2str(ns) ' in varianms.m....'])
        [image(:,:,ns),kfilter]=blur(image(:,:,ns),kfilter);
    end
    
    %redistribute the dimensions
    image=reshape(image,si_image);
    clear smi;
end

%check if image should be upsampled
if any(strcmp(varargin,'upsample')) || any(strcmp(varargin,'1k')) || any(strcmp(varargin,'1K'));
    si_image=size(image);
    numsl=numel(image)/(si_image(1)*si_image(2));
    if any(strcmp(varargin,'1k')) || any(strcmp(varargin,'1K'));
        interpfactor(1)=1024/si_image(1);
        interpfactor(2)=1024/si_image(2);
        newi=1024;
        newj=1024;
    else
        interpfactor=[2 2];
        newi=interpfactor(1)*si_image(1);
        newj=interpfactor(2)*si_image(2);
    end
    smi=zeros(newi, newj,numsl); %interpolated search matrices
    [XX,YY]=meshgrid((0:si_image(2)),(0:si_image(1))');
    for ns=1:numsl;
        display(['Upsampling image ' num2str(ns) ' in varianms.m....'])
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

out.image=image;
out.pars=pars;

%make a surveyplot
if any(strcmp(varargin,'plot'));
   h=slice_surveyplot([{image} varargin]);
end
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
    if tdshift>=0;
        dummy(1:ntd-tdshift,:)=acquireddata(1+tdshift:ntd,:);
    else
        tdshift=-tdshift;
        dummy(1+tdshift:ntd,:)=acquireddata(1:ntd-tdshift,:);
    end
    acquireddata=dummy;
    clear dummy;
end

%rearrange the data according to what pulse sequence it is
si=size(acquireddata); %td x number of traces
td=si(1);
switch lower(pc.seqfil)
    case {'sems'; 'semsdw'; 'stems';'stemsdw';'gems';'gemsvc'};%%%%%%%%%%%%%%%%%%%%%% single echo/FID multislice
        rem=numel(acquireddata)/si(1)/pc.ns/pc.nv;
        
        %time phase slice block/array
        if strcmp(pc.seqcon,'nssnn');
            acquireddata=reshape(acquireddata,[si(1) pc.ns*rem pc.nv]);
            acquireddata=permute(acquireddata,[1 3 2]);
        else
            acquireddata=reshape(acquireddata,[si(1) pc.ns pc.nv rem]);  %seqcon=nccnn
            acquireddata=permute(acquireddata,[1 3 2 4]);
        end
    case {'mems'; 'mgems';'mems_us';'mgems_us';'dwsemgems_us';'dwgesege_us'; ...
			'dwgesege_us_eq';'dwgesege_us_drodpedsl'}%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% multiecho multislice (T2-weighting)
        rem=numel(acquireddata)/si(1)/pc.ne/pc.ns/pc.nv; %remainder
        
        %check if the data came in blocks, if it did, that is an arrayed
        %experiment, and varian stacks this into the second dimension!
        if rem==header.nblocks; %the FID is blocked, and varian writes blocks into the second dimension
            acquireddata=reshape(acquireddata,[si(1) pc.ne pc.ns pc.nv rem]);  %seqcon = ccsnn (echo slice ph1 ph2 ph3)
            acquireddata=permute(acquireddata,[1 4 3 2 5]);
        else
            acquireddata=reshape(acquireddata,[si(1) pc.ne pc.ns pc.nv rem]);  %seqcon = ccsnn (echo slice ph1 ph2 ph3)
            acquireddata=permute(acquireddata,[1 4 3 2 5]);
        end
        %time phase slice echo block/array
        
    case {'fsems'; 'fsemsdw'}  %%%%%%%%%%%%%%%%%%%%%%%%%%% fast spin echo (multiecho, different k each echo) multislice
        rem=numel(acquireddata)/si(1)/pc.ns/pc.nv;
        acquireddata=reshape(acquireddata,[si(1) pc.etl pc.ns pc.nv/pc.etl rem]);  %seqcon = nccnn
        %put the phase dimensions next to each other
        acquireddata=permute(acquireddata, [1 2 4 3 5 6] );
        %stack the phase dimension
        si=[size(acquireddata) 1 1 1];
        acquireddata=reshape(acquireddata, [si(1) si(2)*si(3) si(4:end)]);
        %Now reorder the k space data
        if ~isfield(pc,'pelist') && isfield(pc,'petable');
            dummy=(importdata([pvnmr 'tablib/' pc.petable], ' ',1));
            pc.pelist=(dummy.data)';
            pc.pelist=pc.pelist(:);
        else
            %display('FSE order in .pars.pelist');
        end
        acquireddata(:,pc.pelist+pc.nv/2,:,:,:,:)=acquireddata;
    case {'tagcine'; 'tagcinevc'};
        rem=numel(acquireddata)/si(1)/pc.ne/pc.nv;
        acquireddata=reshape(acquireddata,[si(1) pc.ne pc.nv rem 1 1]);
        acquireddata=permute(acquireddata, [1 3 2 4 5]);
    otherwise
end

%test: reverse phase domian data; result: don't do it
%si=size(acquireddata);
%acquireddata=acquireddata(:,si(2):-1:1,:,:,:);

%frequency shift if desired
%check if a frquency offset [,'freq', frequency] was handed down in varargin
if any(strcmp(varargin,'freq'));
	dii=find(strcmp(varargin,'freq'));
	if numel(varargin)>dii
		if isfloat(varargin{dii+1});
			freqoffset=varargin{dii+1};
		else
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

%set the size of zero filling, if any
mrows=td;
ncols=pc.nv;
if any(strcmp( 'zf',varargin)) || any(strcmp( 'zf2',varargin))  || any(strcmp( 'iso',varargin));  %zerofill called for
    %1. the raw length in the time donain and pe direction is given by
    mrows=td;
    ncols=pc.nv;
    %keep ISO, no additional zero filling
    if any(strcmp( 'iso',varargin));
        maxsize=max(td,pc.nv);
        mrows=maxsize;
        ncols=maxsize;
        fff=1;
    end
    
    %zerofill factor 1.5 in each direction
    if any(strcmp( 'zf',varargin));
        fff=1.5;
    end
    %zerofill by factor 2 in each direction
    if any(strcmp( 'zf2',varargin));
        fff=2;
    end
    mrows=mrows*fff;
    ncols=ncols*fff;
end

if any(strcmp(varargin,'512'));
    mrows=512;
    ncols=512;
end

if any(strcmp(varargin,'1024'));
    mrows=1024;
    ncols=1024;
end

%initialize the phase warp matrix to account for an off-center echo
YY=exp(1i*(-mrows/2:mrows/2-1)'*ones(1,ncols)/mrows*kshift);

%initialize the picture array and Fourier transform with zerofilling;
pic=zeros([mrows ncols si(3:end)]);
siacq=size(acquireddata);
ioff=mrows/2-siacq(1)/2;
joff=ncols/2-siacq(2)/2;
for jj=1:nimages;
    paddedslice=zeros(mrows,ncols);
    paddedslice(ioff+(1:siacq(1)), joff+(1:siacq(2)))=acquireddata(:,:,jj);
    %experimental: virtual reduction of bandwidth by a factor 2
    if any(strcmp('BW',varargin));
        windowsize=2;
        paddedslice=filter(ones(1,windowsize)/windowsize,1,paddedslice);
    end
    pic(:,:,jj)=fftshift(fft2(fftshift(paddedslice)));
	if kshift~=0;
		pic(:,:,jj)=pic(:,:,jj).*YY;
	end
end


%reshape the images into the original size
si=[size(acquireddata) 1];
sipic=size(pic);
pic=reshape(pic, [sipic(1:2) si(3:end)]);

%reorder the slices to be ascending
[sliceposition,orderindex]=sort(pc.pss);
if any(orderindex ~= (1:numel(orderindex)));
    display('varianms.m ordered interleaved slices...');
    display('varianms.m updated the slice positions in out.pars.pss ...');
    pic=pic(:,:,orderindex,:,:,:);
    pc.pss=pc.pss(orderindex);
end



out.image=pic;
out.pars=pc;
out.kspace=acquireddata;

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


function [Mout, kfilter]=blur(varargin)
%[Mout, kfilter]=blur(Min [,kfilter]);
%blurs the image matrix Min using gaussian convolution, 2D fft
%an if optional kfilter matrix has the right size and is nonzero, then the
%input kfilter is used; otherwise kfilter is calculated. The idea is that
%if you want to blur many images, you can reuse the kfilter calculted the
%first time around.
Min=varargin{1};

if ndims(Min)~=2;
    display('blur routine requires matrix input.')
    return
end


if nargin==2;
    kfilter=varargin{2};
    infilterflag=true; %default
    
    %1. check if the sizes of filter and inmatrix Min match
    if any(size(Min)~=size(kfilter));
        infilterflag=false;
    else
        %sizes do match, but also check if there are nonzero entries in the infilter
        if ~any(kfilter(:)~=0);
            infilterflag=false;
        end
    end
else
    infilterflag=false;
end

%if required, (re-)make the kspace filter
if ~infilterflag;
    sfilter=zeros(size(Min));
    [h w]=size(sfilter);
    
    sigma=0.8;
    fw=6;
    x=[-fw:fw];
    y=x';
    fmask=zeros(numel(x));
    for jj=1:length(x);
        for ii=1:length(y);
            fmask(ii,jj)=exp(-(x(jj)^2+y(ii)^2)/(2*sigma^2));
        end
    end
    fmask=fmask/sum(fmask(:));
    
    %fmask is a patch that will go into the middle of sfilter;
    sfilter(h/2+y,w/2+x)=fmask;
    kfilter=fft2(sfilter);
end

FFT_Min=fft2(fftshift(Min));
Mout=ifft2(FFT_Min.*kfilter);
function inputtype=check_img_or_fid(fileidstring)
%check if the filename is image or fid, directory or file
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
            display('First input to varianms must be directory or file of type: ');
            display('*.img, *.fid, *.fdf, fid !') ;
            return
        end
    otherwise
        display('First input must be a directory or a file of type: ');
        display('*.img, *.fid, *.fdf, fid !') ;
        return
end
function h=slice_surveyplot(in)
pic=in{1};
if any(strcmp(in,'bw'));
    cmap=gray(64);
else
    cmap=jet(64);
end
%auto choose the layout
si=size(pic);
maxamp=max(pic(:));
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
        subplot(ny,nx,ii);
        pcolor(abs(pic(:,:,ii)));
        shading flat;
        axis image
        set(gca,'clim',[0 abs(maxamp)],'TickDir', 'in');
        title(['Image #' num2str(ii+(ff-1)*nimperfigure)]);
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
