function DESTE=lig_deste_3d(varargin)
% displacement-encoding-stimulated echoes
% out=lig_deste(DEtimestamp1,DEtimestamp2,DEtimestamp3,MGEtimestamp_refposition, MGEtimestamp other position
% 'plot', [sliceinexvec]


%1. identify the files with the time stamp
if ~isnumeric(varargin{1});
	error('Timestamps identifying DESTE fid-data must be numeric!')
end

%2. Check how many time stamps were handed down
jj=1;
istimestampflag=true;
while istimestampflag
	if jj<= nargin;
		dummy=varargin{jj};
		if isnumeric(dummy) &&   (numel(num2str(dummy)) ==4 || numel(num2str(dummy)) ==3);
			%this is a time stamp
			timestampvec{jj}=['0000' num2str(dummy)];
            timestampvec{jj}=timestampvec{jj}(end-3:end);
			jj=jj+1;
		else
			istimestampflag =false;
		end
	else
		istimestampflag=false;  %ends parsing if the last input argument has been reached
	end
end


%% check if T2star or ge3d or contiguous sems data are being read, if so that is the reference for gridding and masking
hiresflag=false; %=default, reset if more than 3 inputs...
if numel(timestampvec) >3;
    hiresflag=true;
	d=dir;
    cell_of_names={d(:).name};
    
	jj=0;
	%the remaining timestamps must be mge3d t2*data
	for nt=4:numel(timestampvec);
		jj=jj+1;
		
        currenttimestamp=timestampvec{nt};
        
        hiresfileindex = strncmp(cell_of_names,currenttimestamp,4) & ~cellfun('isempty',strfind(cell_of_names,'.fid'));
        hiresname=d(hiresfileindex).name;
        t2starflag=~isempty(strfind(hiresname,'_mge3d'));    %true if its mge3d, false otherwise
        ge3dflag=~isempty(strfind(hiresname,'_ge3d'));       %true if its lig_ge3d, false otherwise
        sems3dflag=~isempty(strfind(hiresname,'_sems'));     %true if its lig_sems, false otherwise
        gems3dflag=~isempty(strfind(hiresname,'_gems'));       %true if it is lig_gems, flase otherwise
        
        %define the hires output filename
        if t2starflag
            HIRES_matfilename=['HIRES_' num2str(currenttimestamp) '_T2star.mat'];
        elseif ge3dflag;
            HIRES_matfilename=['HIRES_' num2str(currenttimestamp) '_magnitude.mat'];
        elseif sems3dflag;
            HIRES_matfilename=['HIRES_' num2str(currenttimestamp) '_sems_magnitude.mat'];
        elseif gems3dflag;
            HIRES_matfilename=['HIRES_' num2str(currenttimestamp) '_gems_magnitude.mat'];
        end
        
        %% read the high resolution data
        if t2starflag
            HIRES(jj)=lig_T2star(currenttimestamp);
        elseif ge3dflag
            HIRES(jj)=lig_ge_3d(currenttimestamp,'vox',0.4);  %coarse grain it to (350 micron)^3 voxels
        elseif sems3dflag
            HIRES(jj)=lig_sems(currenttimestamp,'vox',0.4);  %coarse grain it to (350 micron)^3 voxels
        elseif gems3dflag;
            HIRES(jj)=lig_sems(currenttimestamp,'vox',0.4);  %coarse grain it to (350 micron)^3 voxels
        end
	end
	outgridvec=size(HIRES(1).magnitude);
    display(['High resolution 3d data size: ' num2str(outgridvec,'% d') ]);
end

clear dummy;

lambda=zeros(1,numel(timestampvec));
encodedirectionvector={};

blurvec=[1.5 1.5 1.5];

for nf=1:3
	d=dir;
	time_flag=false(numel(d),1);
	img_flag=time_flag;
	PH_flag=time_flag;
	fid_flag=time_flag;
	
	timestamp=timestampvec{nf};
	
	for jj=1:numel(dir);
		filename=d(jj).name;
		[fpath,fname,fext]=fileparts(filename);
		
		%get or fake the last three characters of the filename to identfy phase
		%files
		if numel(fname)>2;
			fnamelast3characters=fname((length(fname)-2):length(fname));
		else
			fnamelast3characters='xxx';
		end
		
		time_flag(jj)=~isempty(strfind(fname,timestamp));
		img_flag(jj)=~isempty(strfind(fext,'img'));
		PH_flag(jj)=strcmp(fnamelast3characters,'_PH');
		fid_flag(jj)=~isempty(strfind(fext,'fid'));
	end
	
	fid_directory=d(time_flag & fid_flag).name;
	display(['Reconstructing data from fid, timestamp =' num2str(timestamp)]);
    
    fid_dir{nf}=fid_directory;
    
    %read the phase encoding stimulated echoes
    % VARIANMS VARIANMS VARIANMS
    prelim(nf)=varianms(fid_directory, 'pixvec',[128 32 32],'apod', [58 14],'gradlag',48); %straight up reconstruction from fid, 

    
%     % Legacy: adaptive blurring based on lambda, and a blurfactor2
%     blurfactor=2;
%     blurvec=blurfactor*...                    %default blurfactor is 2, the factor 10 is for conversion to mm
%         [prelim(nf).pars.dro(1) prelim(nf).pars.dpe(1) prelim(nf).pars.dsl(1)]*10*prelim(nf).pars.lambda + blurvec;
%     

    timestampstr=['000' num2str(timestampvec{nf})];
    timestampstr=timestampstr(end-3:end);
    tssvec(nf,:)=timestampstr;
end


if any(strcmp(varargin,'blurvec'));
    ind=find(strcmp(varargin,'blurvec'));
    blurvec=varargin{ind+1};
    bvn=['_blurvec' num2str(blurvec(1)) num2str(blurvec(2)) num2str(blurvec(3))];
else
    bvn='';
end

%% optional error diagnosing control module
cmflag=false;
if cmflag;
	prelim=lig_diagnose_phase_magnitude_consistency(prelim);
end
	
%% Axis check
minax1=inf; minax2=inf;minax3=inf;maxax1=-inf;maxax2=-inf;maxax3=-inf;
for nf=1:3 
    center(nf,1)=prelim(nf).pars.pro;
    center(nf,2)=prelim(nf).pars.ppe;
    center(nf,3)=prelim(nf).pars.pss0;
    FOV(nf,1)=prelim(nf).pars.lro;
    FOV(nf,2)=prelim(nf).pars.lpe;
    FOV(nf,3)=(prelim(nf).pars.ns-1)*prelim(nf).pars.thk/10;
    si(nf,:)=size(prelim(nf).image);
    
    minax1=min(minax1,min(prelim(nf).pars.axis1));
    minax2=min(minax2,min(prelim(nf).pars.axis2));
    minax3=min(minax3,min(prelim(nf).pars.axis3));
    maxax1=max(maxax1,max(prelim(nf).pars.axis1));
    maxax2=max(maxax2,max(prelim(nf).pars.axis2));
    maxax3=max(maxax3,max(prelim(nf).pars.axis3));
end
minFOV=min(FOV);
display(FOV);

% Compare FOV of phase data with the FOV of HIRES data, if they were
% acquired
if hiresflag && ~sems3dflag
    hicen(1)=HIRES(1).params.pro;
    hicen(2)=HIRES(1).params.ppe;
    hicen(3)=HIRES(1).params.ppe2;
    hifov(1)=HIRES(1).params.lro;
    hifov(2)=HIRES(1).params.lpe;
    hifov(3)=HIRES(1).params.lpe2;
    display('HIRES FOV:');
    display(hifov);
elseif hiresflag && sems3dflag
	hicen(1)=HIRES(1).params.pro;
    hicen(2)=HIRES(1).params.ppe;
    hicen(3)=HIRES(1).params.pss0;
    hifov(1)=HIRES(1).params.lro;
    hifov(2)=HIRES(1).params.lpe;
    hifov(3)=(HIRES(1).params.ns-1)*HIRES(1).params.thk/10;
    display('HIRES FOV:');
    display(hifov);
end
    
%find the data set containing the reference axes
for nf=1:3;
    refdata(nf)=all(minFOV==FOV(nf,:));
end
nfref=find(refdata,1);  %the index to the data containing the reference axes
offref=find(~refdata);  %the indeices to data that are off FOV

%cycle through the data files if any field of views are to be adjusted
if any(~refdata);
    for nf=offref;
        timestampstr=['000' num2str(timestampvec{nf})];
        timestampstr=timestampstr(end-3:end);
        tssvec(nf,:)=timestampstr;
        
        display(['adjust resolution of file ' timestampstr]);
        pixvec=2 * round(si(nf,1:3).*FOV(nf,1:3)./minFOV(1:3)/2+eps);
        prelim(nf)=varianms(fid_dir{nf}, 'apod','pixvec',pixvec); %straight up reconstruction from fid
        
        
        % axis1
        nref =   numel(prelim(nfref).pars.axis1);
        nnew =   numel(prelim(nf).pars.axis1);
        if ~(nref==nnew);
            for jj=0:(nnew-nref);
                indexvec=jj+(1:nref);
                mintest(jj+1)=sum((prelim(nfref).pars.axis1-prelim(nf).pars.axis1(indexvec)).^2);
            end
            [minmismatch,index]=min(mintest);
            indexvec=index-1+(1:nref);
            prelim(nf).pars.axis1=prelim(nf).pars.axis1(indexvec);
            prelim(nf).image=prelim(nf).image(indexvec,:,:,:);    
        end
        
        %axis2
        nref =   numel(prelim(nfref).pars.axis2);
        nnew =   numel(prelim(nf).pars.axis2);
        if ~(nref==nnew);
            for jj=0:(nnew-nref);
                indexvec=jj+(1:nref);
                mintest(jj+1)=sum((prelim(nfref).pars.axis2-prelim(nf).pars.axis2(indexvec)).^2);
            end
            [minmismatch,index]=min(mintest);
            indexvec=index-1+(1:nref);
            prelim(nf).pars.axis2=prelim(nf).pars.axis2(indexvec);
            prelim(nf).image=prelim(nf).image(:,indexvec,:,:);     
        end
        
        %axis3
        nref =   numel(prelim(nfref).pars.axis3);
        nnew =   numel(prelim(nf).pars.axis3);
        if ~(nref==nnew);
            for jj=0:(nnew-nref);
                indexvec=jj+(1:nref);
                mintest(jj+1)=sum((prelim(nfref).pars.axis3-prelim(nf).pars.axis3(indexvec)).^2);
            end
            [minmismatch,index]=min(mintest);
            indexvec=index-1+(1:nref);
            prelim(nf).pars.axis3=prelim(nf).pars.axis3(indexvec);
            prelim(nf).image=prelim(nf).image(:,:,indexvec,:);    
        end
    end
end

% End of Axis check

%% Initialize some variables for storage
si=size(prelim(1).image);				%[n1 n2 n3 enc+-] 
orig_phasefactor=zeros([si(1:3) 3]);
orig_compleximage=zeros([si(1:3) 3]);
if hiresflag
    QC.data=zeros([outgridvec 2 3]);					%[n1 n2 n3 enc+- encdir3]
    gridvec=outgridvec;
else
    QC.data=zeros([si 2 3]);
    gridvec=si(1:3);
end

%% Loop over the 3 phase encoded data, spatial averaging and upsampling
for nf=1:3
    % pick the varianms structure and work with that
	out=prelim(nf);

	orig_phasefactor(:,:,:,nf)=out.image(:,:,:,1)./out.image(:,:,:,2);
	orig_phasefactor(:,:,:,nf)=orig_phasefactor(:,:,:,nf)./abs(orig_phasefactor(:,:,:,nf));
	orig_compleximage(:,:,:,nf)=mean(abs(out.image),4).*orig_phasefactor(:,:,:,nf);  %averages over plus minus lambda
	orig_kspace(:,:,:,:,nf)=out.kspace;
	si=size(orig_compleximage(:,:,:,nf));
    
    	
    %% blur 3d
    %get the voxel size of the original stimulated echo data;
    acquiredvoxelsize=abs([out.pars.axis1(1)-out.pars.axis1(2) ...
                            out.pars.axis2(1)-out.pars.axis2(2) ...
                            out.pars.axis3(1)-out.pars.axis3(2)]);
    

    %set the voxelsize, if the data is upsampled it must be at least as
    %large as the original voxel size
    voxelsize=max(0.6,max(acquiredvoxelsize));
	%blurvec= [1 voxelsize voxelsize]; %in mm
    
    %blurvec=[1 0.5 0.75];
	
    %create a temporarystucture with fields image and parss, to blur using mm voxel sizes
	temp.image=orig_compleximage(:,:,:,nf);
	temp.pars=out.pars;
	temp.acquisitiontype='2D';
	
	% blur the STE, before any resampling
	temp=blur3d(temp,'vox', blurvec);
	display(['blur size= ' num2str(blurvec) 'mm']); %#ok<*DISPLAYPROG>
    
	% upsample without blurring
    temp=blur3d(temp,'vox',[0 0 0],'grid',gridvec);
	
	%update the complex image and phase factor to the new grid size
	compleximage(:,:,:,nf)=temp.image;
	
	
	%average the quality control data over the same volume
	qc_temp=blur3d(out,'vox', blurvec);  %blurs out, the displacement encoding STE
	qc_temp=blur3d(qc_temp,'vox',[0 0 0],'grid',gridvec); %upsamples the displacement encoding STE
	
	QC.data(:,:,:,:,nf)=qc_temp.image;
	QC.axis1{nf}=newaxis(out.pars.axis1,gridvec(1))';
	QC.axis2{nf}=newaxis(out.pars.axis2,gridvec(2))';
	QC.axis3{nf}=newaxis(out.pars.axis3,gridvec(3))';
	
	%get the paramters out of the fid directory
	pars=out.pars;
    origpars(nf)=pars;
	
	%identify which of the three displacement encoding directions was
	%active, read, phase or slice, generate a filename for the output file
	
	dd(1)=(pars.dro(1)-pars.dro(2));  %this is either zero or two
	dd(2)=(pars.dpe(1)-pars.dpe(2));  %this is either zero or two
	dd(3)=(pars.dsl(1)-pars.dsl(2));  %this is either zero or two
	
	directionindicator = {'RO'; 'PE'; 'SL'};

    
	for jj= 1:3;
		if dd(jj)~= 0
			lambda(nf)=round(10*pars.lambda/dd(jj) *1000);
			encodedirectionvector{nf}=directionindicator{jj};
			encodedirectionindex(nf)=jj;
			%outfilename=['analyzed_' num2str(timestamp) '_' directionindicator{jj} '_lambda_' num2str(lambda(jj)) 'mu.mat'];
		end
	end
	wigglecomment{nf}=[tssvec(nf,:) ': ' pars.comment];
	
	parspars{nf}=pars;
end


%% SORTING, in case the time stamps are out of order
%  this is where the phase data gets sorted into RO, PE, SL
[dummy, phaseindvec]=sort(encodedirectionindex);   %
compleximage=compleximage(:,:,:,phaseindvec);
orig_compleximage=orig_compleximage(:,:,:,phaseindvec);
lambda=lambda(phaseindvec);
orig_kspace=orig_kspace(:,:,:,:,phaseindvec);
QC.data=QC.data(:,:,:,:,phaseindvec);
QC.comment={'data field: nro nph nslice 2_values_encoding_polarity 3_encoding_directions'; ...
	'axis1 readout axis2 phse axis3 slice'};


si=size(compleximage);
absvalueimage=zeros(si(1:3));
for nf=1:3;
	absvalueimage=absvalueimage+abs(compleximage(:,:,:,nf));
	dumdum1{nf}=wigglecomment{phaseindvec(nf)};
	dumdum2{nf}=timestampvec{phaseindvec(nf)};
	dummy_axis1{nf}=QC.axis1{phaseindvec(nf)};
	dummy_axis2{nf}=QC.axis2{phaseindvec(nf)};
	dummy_axis3{nf}=QC.axis3{phaseindvec(nf)};
end
wigglecomment=dumdum1;
timestampvec=dumdum2;
QC.axis1=dummy_axis1;
QC.axis2=dummy_axis2;
QC.axis3=dummy_axis3;
% *********************************************************


%%
zfis=size(compleximage(:,:,:,1));  % zero filled image size
DESTE.axis1=newaxis(pars.axis1,zfis(1)); %in mm
DESTE.axis2=newaxis(pars.axis2,zfis(2)); %in mm
DESTE.axis3=newaxis(pars.axis3,zfis(3)); %in mm
DESTE.blurvec=blurvec;
DESTE.data=compleximage;

[signallevel,noiselevel] = estimate_noiselevel(absvalueimage);
DESTE.mask=1*(absvalueimage>2*noiselevel);  %numeric mask, ones for good pixels, zeros for insufficient signal
DESTE.QC=QC;



%% Create a high resolution based mask on the grid defined by the stimulated echo DESTE acquisition
%  Then run the movie comparing the maskrd stimulated echo amplitude (left panel) with and the HIRES amplitude (right panel) 
if exist('HIRES','var')
    
    [X1ref,X2ref,X3ref]=ndgrid(HIRES(1).axis1,HIRES(1).axis2,HIRES(1).axis3);
    [X1q,X2q,X3q]=ndgrid(DESTE.axis1,DESTE.axis2,DESTE.axis3);
    DESTE.mask_hr = +interpn(X1ref,X2ref,X3ref,HIRES(1).mask,X1q,X2q,X3q,'nearest',0);  %
    DESTE.mask_hr(DESTE.mask_hr>0.5)=1;
    DESTE.mask_hr(DESTE.mask_hr<=0.5)=0;
	
	dummymask=DESTE.mask_hr.*DESTE.mask;
    
    figure('position',[100         100        1500         600]);
    si=size(DESTE.mask);
    for jj=1:si(1);
        subplot(1,2,1);
        imagesc(DESTE.axis3, DESTE.axis2,log10(squeeze(absvalueimage(jj,:,:).*dummymask(jj,:,:))+1),[0 5]);
        axis square;
        grid on;
        
        subplot(1,2,2);
        imagesc(HIRES(1).axis3,HIRES(1).axis2,log10(squeeze(HIRES(1).magnitude(jj,:,:))),[0 6]);
        axis square;
        grid on;
        pause(0.05);
    end
end

%% now calculate the numeric drivatives for the strain image
si=size(compleximage);
complex_strainimage=zeros([si 3]);
ind1=[2:si(1)-1];
ind2=[2:si(2)-1];
ind3=[2:si(3)-1];

phaseimage=compleximage./abs(compleximage);

% RO encode            lambda  subtraction
complex_strainimage(ind1,:,:,1,1)=(phaseimage(ind1+1,:,:,1)./phaseimage(ind1-1,:,:,1));
complex_strainimage(:,ind2,:,1,2)=(phaseimage(:,ind2+1,:,1)./phaseimage(:,ind2-1,:,1));
complex_strainimage(:,:,ind3,1,3)=(phaseimage(:,:,ind3+1,1)./phaseimage(:,:,ind3-1,1));
% PE encode          lambda  subtraction
complex_strainimage(ind1,:,:,2,1)=(phaseimage(ind1+1,:,:,2)./phaseimage(ind1-1,:,:,2));
complex_strainimage(:,ind2,:,2,2)=(phaseimage(:,ind2+1,:,2)./phaseimage(:,ind2-1,:,2));
complex_strainimage(:,:,ind3,2,3)=(phaseimage(:,:,ind3+1,2)./phaseimage(:,:,ind3-1,2));
% SL strain         lambda  subtraction
complex_strainimage(:,:,ind3,3)=(phaseimage(:,:,ind3+1,3)./phaseimage(:,:,ind3-1,3));

complex_strainimage(ind1,:,:,3,1)=(phaseimage(ind1+1,:,:,3)./phaseimage(ind1-1,:,:,3));
complex_strainimage(:,ind2,:,3,2)=(phaseimage(:,ind2+1,:,3)./phaseimage(:,ind2-1,:,3));
complex_strainimage(:,:,ind3,3,3)=(phaseimage(:,:,ind3+1,3)./phaseimage(:,:,ind3-1,3));

delta(1)=DESTE.axis1(2)-DESTE.axis1(1);
delta(2)=DESTE.axis2(2)-DESTE.axis2(1);
delta(3)=DESTE.axis3(2)-DESTE.axis3(1);

strainimage=angle(complex_strainimage);
for encoding_directions =1:3;
	for gradient_directions =1:3;
		strainimage(:,:,:,encoding_directions,gradient_directions)=...
			strainimage(:,:,:,encoding_directions,gradient_directions)/(2*pi)*lambda(encoding_directions)/1000/(2*delta(gradient_directions));
	end
end

%% Now calculate the Lagrangian strains using Callan's derivs2strains
% 1. smooth the strain image using blur3d, over the same voxel size used to
% smooth the complex data
clear dummy;
dummy1.image=strainimage;
dummy1.pars=out.pars;
dummy1.pars.axis1=DESTE.axis1;
dummy1.pars.axis2=DESTE.axis2;
dummy1.pars.axis3=DESTE.axis3;
dummy1.acquisitiontype=out.acquisitiontype;
dummy1=blur3d(dummy1,'vox',voxelsize*[1 1 1]);
strainimage=dummy1.image;

[Lagrange, Ql, Vl]=derivs2strain(strainimage);

DESTE.strains=strainimage;
DESTE.Lagrange=Lagrange;
DESTE.Q=Ql;
DESTE.V=Vl;
DESTE.lambda=lambda/1000;   %in mm
DESTE.comment=wigglecomment;
DESTE.timest=timestampvec;

%% check if a wiggle or stretch was specified in the comment
wiggleword={'wiggle','stretch'};
for www=1:numel(wiggleword);
    
    for jj=1:3;
        comment=wigglecomment{jj};
        if ~isempty(strfind(comment,wiggleword{www}));
            numbers_index = isstrprop(comment,'digit') | isstrprop(comment,'punct');
            minussignindex=strfind('comment','-');
            if any(minussignindex);
                numbers_index(minussignindex)=true;
            end
            totherightofwiggle=false(1,numel(numbers_index));
            wiggleindex=strfind(comment,wiggleword{www});
            totherightofwiggle((wiggleindex+length(wiggleword{www})):(wiggleindex+length(wiggleword{www})+5))=true;
            dummy=false(size(totherightofwiggle));
            
            dummy(1:numel(numbers_index))=numbers_index;
            numbers_index=dummy;
            numbers_index=numbers_index & totherightofwiggle;
            wiggle(jj)=str2num([comment(numbers_index)]);
        end
    end
   
end

if exist('wiggle','var');
     DESTE.wiggle=wiggle;
end


%% check if a position was specified in the comment
for jj=1:3
	comment=wigglecomment{jj};
	if ~isempty(strfind(comment,'pos'));
		numbers_index = isstrprop(comment,'digit') | isstrprop(comment,'punct');
        minussignindex=strfind('comment','-');
        if any(minussignindex);
            numbers_index(minussignindex)=true;
        end
		totherightofpos=false(1,numel(numbers_index));
		posindex=strfind(comment,'pos');
		totherightofpos((posindex+3):(posindex+8))=true;
        
        totherightofpos=totherightofpos(1:numel(numbers_index));
        
		numbers_index=numbers_index & totherightofpos;
		position{jj}=str2num([comment(numbers_index)]);
	else
		position{jj}=[];
	end
end

DESTE.position=position;
DESTE.note='origdata is constructed from +-lambda, at the true resolution';
DESTE.origdata=orig_compleximage;
DESTE.origkspace=orig_kspace;
DESTE.origpars=origpars;


%% quick view of colored phases
slvec=1:si(3);

figure('position',[50 50 900 900]);

for sl=slvec
    nanmask=DESTE.mask(:,:,sl);
    nanmask(nanmask==0)=NaN;
    
    
    for endir=1:3;
        %			tsubplot(3,4,(endir-1)*4+1,10);
        subplot(1,3,endir);
        imagesc(DESTE.axis2,DESTE.axis1,angle(compleximage(:,:,sl,endir)).*nanmask,[-pi pi]);
        axis image;
        ylabel('1-dir')
        xlabel('2-dir')
        title(['encode dir ' num2str(endir)]);
    end
    
    pause(0.1);
end



%%
if exist('HIRES','var')
    DESTE.HIRES=HIRES(1);
    outfilename=['DESTE_strains_' DESTE.timest{1} '_' DESTE.timest{2} '_' DESTE.timest{3} '_' num2str(HIRES(1).timest) bvn '.mat'];

    
        
    
%     si=size(HIRES(1).magnitude);
%     
%     high_res_strain =zeros([si 3 3]);
%     
%     
%     %%
%     ordercontrol=false;
%     if ordercontrol;
%         hir_absvalueimage=interp3(X2q,X1q,X3q,absvalueimage,X2ref,X1ref,X3ref);
%         
%         figure('position',[400         563        1097         406], 'Name', 'Coronal view');
%         sir=size(hir_absvalueimage);
%         for jj=1:sir(3);
%             subplot(1,2,1);
%             imagesc(DESTE.HIRES.axis2,DESTE.HIRES.axis1, log10(hir_absvalueimage(:,:,jj)),[1 5]);
%             axis square;
%             
%             subplot(1,2,2); imagesc(DESTE.HIRES.axis2,DESTE.HIRES.axis1, log10(DESTE.HIRES.magnitude(:,:,jj)),[1 5]);
%             set(gcf,'Name',num2str(jj));
%             axis square;
%             pause(0.1);
%         end
%         
%     end
    %%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%% INTERPOLATION OF THE STRAIN TO RESOLUTION OF T2STAR DATA (OPTIONAL)
%     hires_strainflag= false;
%     if hires_strainflag;
%         for ii=1:3;
%             for jj=1:3;
%                 hi_res_strain(:,:,:,ii,jj)=interp3(X2q,X1q,X3q,squeeze(DESTE.strains(:,:,:,ii,jj)),X2ref,X1ref,X3ref);
%             end
%         end
%         DESTE.hi_res_strain=hi_res_strain;
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    outfilename=['DESTE_strains_' DESTE.timest{1} '_' DESTE.timest{2} '_' DESTE.timest{3} '.mat'];
end

save(outfilename, '-struct', 'DESTE');

if any(strcmp(varargin,'view3d'));
    view3dgui(cat(2,angle(compleximage(:,:,:,1)),angle(compleximage(:,:,:,2)),angle(compleximage(:,:,:,3))),'mlim',[-pi pi]);
    view3dgui(cat(3,angle(compleximage(:,:,:,1)),angle(compleximage(:,:,:,2)),angle(compleximage(:,:,:,3))),'mlim',[-pi pi]);
    
end

end

%%

function [signallevel, noiselevel]=estimate_noiselevel(varargin)
incube=varargin{1};
%estimate the noise using NxN patches in dimensions 1 & 2 of the dataset
%default 8x8 pixel patch
%default estimate for empty space in the data set = 25%
N=8;                    %the size of a patch of pixels used to coarse-grain the image;
emptyspacefraction=0.25;%estimate of how much empty space there is in the image data, for the noise calculation

si=size(incube);
incube=reshape(incube,[si(1) si(2) numel(incube)/(si(1)*si(2))]);
refplane=incube(:,:,(1 : end));
si=size(refplane);
ni=floor(si(1)/N);
nj=floor(si(2)/N);
nk=floor(si(3)/N);
noisetest=zeros(ni,nj,nk);
amplitudetest=zeros(ni,nj,nk);
for kk=0:nk-1;
    for ii=0:ni-1;
        for jj=0:nj-1;
            ivec=(1:N) + ii*N;
            jvec=(1:N) + jj*N;
            kvec=(1:N) + kk*N;
            patch=refplane(ivec,jvec,kvec);
            noisetest(ii+1,jj+1,kk+1)=sqrt(var(abs(patch(:))));
            amplitudetest(ii+1,jj+1,kk+1)=abs(mean(patch(:)));
        end
    end
end

sortednoiselevel=sort(noisetest(:));
sortedamplitudelevel=sort(amplitudetest(:));

% the noise level is calculated where the signal is zero,
% (or finite for abs value data
% by taking the mean of the first 25% percent of noise (asumed empty)
% and the mean of the first 25% of the signal
noiselevel= mean(sortednoiselevel(1:round(numel(sortednoiselevel)*emptyspacefraction))) + ...
    mean(sortedamplitudelevel(1:round(numel(sortedamplitudelevel)*emptyspacefraction)));


%signal level is calculate at the end of the signal distribution
signallevel=mean(sortedamplitudelevel(round(end*0.95):(end-10)));

if any(strcmp(varargin,'plot'));
    figure;
    plot(sortedamplitudelevel,'r');
    hold on;
    plot(sortednoiselevel,'b');
    set(gca,'Yscale','log');
    legend({'signal'; 'noise'});
    legend boxoff;
    title(['S=' num2str(round(signallevel)) ' / N=' num2str(round(noiselevel))]);
end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outaxis=newaxis(oldaxis,numberofentries)

spacevec=diff(oldaxis);
if any(abs(spacevec/spacevec(1)-1)>1e-10);
	display('wonky axes with unven spacing');
end
spacing=mean(spacevec);
FOV=spacevec(1)*numel(oldaxis);
Limit1=oldaxis(1)-spacing/2;      %this accounts for the voxel size
Limit2=oldaxis(end)+spacing/2;

if abs(abs(FOV)-abs(Limit1-Limit2))/abs(FOV)>1e-10;
	display('Error in axis calculation');
	return
end
newspacing=(Limit2-Limit1)/numberofentries;
outaxis=Limit1+[0.5:1:numberofentries-0.5]*newspacing;

end








