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
		if isnumeric(dummy) &&     (numel(num2str(dummy)) ==4 || numel(num2str(dummy)) ==3);
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


%% check if T2star or ge3d data are being read, if so that is the reference for gridding and masking
hiresflag=false;
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
        t2starflag=~isempty(strfind(hiresname,'mge3d')); %true if its mge3d, false otherwise
        if t2starflag;
            HIRES_matfilename=['HIRES_' num2str(currenttimestamp) '_T2star.mat'];
        else
            HIRES_matfilename=['HIRES_' num2str(currenttimestamp) '_magnitude.mat'];
        end
        
        %% Always recalculate the T2 star data until troubleshooting of varianms is done!
        if false;
			HIRES(jj)=load(HIRES_matfilename);
        else
            if t2starflag;
                HIRES(jj)=lig_T2star(currenttimestamp);
            else
                HIRES(jj)=lig_ge_3d(currenttimestamp);
            end
		end
	end
	outgridvec=size(HIRES(1).magnitude);
    display(['A high resolution mge3d was acquired with td=' num2str(outgridvec(1)) '.' ]);
    
end

clear dummy;

lambda=zeros(1,numel(timestampvec));
encodedirectionvector={};
for nf=1:3;
	
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
	if ~hiresflag;
    prelim(nf)=varianms(fid_directory, 'apod'); %straight up reconstruction from fid
    else
        prelim(nf)=varianms(fid_directory, 'apod','pixvec',outgridvec); %straight up reconstruction from fid
    end
    
    
    timestampstr=['000' num2str(timestampvec{nf})];
    timestampstr=timestampstr(end-3:end);
    tssvec(nf,:)=timestampstr;
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Quality control check to see if all the phase data were acquired on the same grid
for nf=1:3; 
    center(nf,1)=prelim(nf).pars.pro;
    center(nf,2)=prelim(nf).pars.ppe;
    center(nf,3)=prelim(nf).pars.pss0;
    FOV(nf,1)=prelim(nf).pars.lro;
    FOV(nf,2)=prelim(nf).pars.lpe;
    FOV(nf,3)=(prelim(nf).pars.ns-1)*prelim(nf).pars.thk;
    si(nf,:)=size(prelim(nf).image);
end
minFOV=min(FOV);
display(FOV);

%find the data set containing the reference axes
for nf=1:3;
    refdata(nf)=all(minFOV==FOV(nf,:));
end
nfref=find(refdata,1);  %the index to the data containing the reference axes
offref=find(~refdata);  %the indeices to data that are off FOV

%cycle through the data files
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



for nf=1:3;
    out=prelim(nf);
	%get and reconstruct the multislice image from the fid directory
	%the image is 4-d three spatial dimensions and plus minus lambda.

	
	orig_phasefactor(:,:,:,nf)=out.image(:,:,:,1)./out.image(:,:,:,2);
	orig_phasefactor(:,:,:,nf)=orig_phasefactor(:,:,:,nf)./abs(orig_phasefactor(:,:,:,nf));
	orig_compleximage(:,:,:,nf)=mean(abs(out.image),4).*orig_phasefactor(:,:,:,nf);  %averages over plus minus lambda
	si=size(orig_compleximage(:,:,:,nf));
    
    %upsample and smooth the 2D multislice phase data
    if hiresflag;   %if hires data was acquired:
        %1. match the time domain dimensions of 2D multislice to that of hires 3D mge3d
        %2. double resolution in the phase and slice directions
        %gridvec=[outgridvec(1) 2*si(2) 2*si(3)];
        gridvec=outgridvec;
    else            %otherwise
        %1. keep the timedomain as it came up
        %2. double resolution in the phase and slice directions
        gridvec=[2*si(1) 2*si(2) 2*si(3)]; %zerofill in all directions
    end
	
	display(['zerofilled STE image to size [' num2str(gridvec(1)) ' ' num2str(gridvec(2)) ' ' num2str(gridvec(3)) '];']);
	
    %blur 3d
	voxelsize=0.4;
	display(['blur size=' num2str(voxelsize) 'mm']);
	dummy.image=orig_compleximage(:,:,:,nf);
	dummy.pars=out.pars;
	dummy=lig_blur3d(dummy,'vox',voxelsize*[1 1 1],'grid',gridvec);
	
	%update the complex image and phase factor to the new grid size
	compleximage(:,:,:,nf)=dummy.image;
	
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


% SORTING, in case the time stamps are out of order
%this is where the phase data gets sorted into RO, PE, SL
[dummy, phaseindvec ]=sort(encodedirectionindex);   %

compleximage=compleximage(:,:,:,phaseindvec);
lambda=lambda(phaseindvec);

si=size(compleximage);
absvalueimage=zeros(si(1:3));
for nf=1:3;
	absvalueimage=absvalueimage+abs(compleximage(:,:,:,nf));
	dumdum1{nf}=wigglecomment{phaseindvec(nf)};
	dumdum2{nf}=timestampvec{phaseindvec(nf)};
end
wigglecomment=dumdum1;
timestampvec=dumdum2;
%% *********************************************************


%%
zfis=size(compleximage(:,:,:,1));  % zero filled image size
DESTE.axis1=newaxis(pars.axis1,zfis(1)); %in mm
DESTE.axis2=newaxis(pars.axis2,zfis(2)); %in mm
DESTE.axis3=newaxis(pars.axis3,zfis(3)); %in mm
DESTE.data=compleximage;

[signallevel,noiselevel] = estimate_noiselevel(absvalueimage);
DESTE.mask=1*(absvalueimage>3*noiselevel);  %numeric mask, ones for good pixels, zeros for insufficient signal

%% 
if exist('HIRES','var');
     figure('position',[100         100        1500         600]);
     si=size(DESTE.mask);
    for jj=1:si(1);
        subplot(1,2,1);
        imagesc(DESTE.axis3, DESTE.axis2,log10(squeeze(absvalueimage(jj,:,:).*DESTE.mask(jj,:,:))+1),[0 5]);
        axis square;
        grid on;
        subplot(1,2,2);
        imagesc(HIRES(1).axis3,HIRES(1).axis2,log10(squeeze(HIRES(1).magnitude(jj,:,:))),[0 6]);
        axis square;
        grid on;
        pause(0.1);
    end
    
   %% 
    diagnosticplot=true;
    if diagnosticplot;
        figure;
        si=size(absvalueimage);
        subplot(1,2,1);
        imagesc(DESTE.axis3,DESTE.axis1,squeeze(absvalueimage(:,si(2)/2,:)));
        
        si=size(HIRES(1).magnitude);
        subplot(1,2,2);
        imagesc(HIRES(1).axis3,HIRES(1).axis1,squeeze(HIRES(1).magnitude(:,si(2)/2,:)));
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

DESTE.strains=strainimage;
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
            numbers_index=numbers_index & totherightofwiggle;
            wiggle(jj)=str2num([comment(numbers_index)]);
        end
    end
end
DESTE.wiggle=wiggle;

%% check if a position was specified in the comment
for jj=1:3;
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
		numbers_index=numbers_index & totherightofpos;
		position{jj}=str2num([comment(numbers_index)]);
	else
		position{jj}=[];
	end
end

DESTE.position=position;
DESTE.note='orig_compleximage is constructed from +-lambda, at true resolution';
DESTE.origdata=orig_compleximage;
DESTE.origpars=origpars;


%% Plot option
if any(strcmp(varargin,'plot'));
	ind=find(strcmp(varargin,'plot'));
	if nargin>ind && isnumeric(varargin{ind+1});
		slvec=varargin{ind+1};
		
	else
		slvec=1:si(3);
	end
	
	for sl=slvec;
		nanmask=DESTE.mask(:,:,sl);
		nanmask(nanmask==0)=NaN;
		figure('name',num2str(sl),'position',[50 50 650 900]);
		
		for endir=1:3;
			tsubplot(3,4,(endir-1)*4+1,10);
			imagesc(DESTE.axis2,DESTE.axis1,angle(compleximage(:,:,sl,endir)).*nanmask);
			axis image;
			ylabel('1-dir')
			xlabel('2-dir')
			title(['encode dir ' num2str(endir)]);
			
			for graddir=1:3;
				
				tsubplot(3,4,(endir-1)*4 + graddir+1,10);
				imagesc(DESTE.axis2,DESTE.axis1,strainimage(:,:,sl,endir,graddir).*nanmask,[-0.5 0.5]);
				axis image;
				ylabel('1-dir')
				xlabel('2-dir')
				title(['(d' num2str(endir) ',g' num2str(graddir) ')']);
			end
		end
		
		
	end
end

%%
if exist('HIRES','var');
	DESTE.HIRES=HIRES(1);
	outfilename=['DESTE_strains_' DESTE.timest{1} '_' DESTE.timest{2} '_' DESTE.timest{3} '_' num2str(HIRES(1).timest) '.mat'];
	
	
	[X1ref,X2ref,X3ref]=ndgrid(HIRES(1).axis1,HIRES(1).axis2,HIRES(1).axis3);
	
	
	[X1q,X2q,X3q]=ndgrid(DESTE.axis1,DESTE.axis2,DESTE.axis3);
	
	%Vq = interp3(X2ref,X1ref,X3ref,refdata,X2q,X1q,X3q);
	
	si=size(HIRES(1).magnitude);
	
	high_res_strain =zeros([si 3 3]);
	
    
    %%
    ordercontrol=true;

	if ordercontrol;
		hir_absvalueimage=interp3(X2q,X1q,X3q,absvalueimage,X2ref,X1ref,X3ref);
		
		figure('position',[400         563        1097         406], 'Name', 'Coronal view');
        sir=size(hir_absvalueimage);
		for jj=1:sir(3);
			subplot(1,2,1);
			imagesc(DESTE.HIRES.axis2,DESTE.HIRES.axis1, log10(hir_absvalueimage(:,:,jj)),[1 5]);
            axis square;
			
			subplot(1,2,2); imagesc(DESTE.HIRES.axis2,DESTE.HIRES.axis1, log10(DESTE.HIRES.magnitude(:,:,jj)),[1 5]);
			set(gcf,'Name',num2str(jj));
            axis square;
			pause(0.1);
        end
		
	end
	%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% INTERPOLATION OF THE STRAIN TO RESOLUTION OF T2STAR DATA (OPTIONAL)
    hires_strainflag= false;
    if hires_strainflag;
        for ii=1:3;
            for jj=1:3;
                hi_res_strain(:,:,:,ii,jj)=interp3(X2q,X1q,X3q,squeeze(DESTE.strains(:,:,:,ii,jj)),X2ref,X1ref,X3ref);
            end
        end
        DESTE.hi_res_strain=hi_res_strain;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
	outfilename=['DESTE_strains_' DESTE.timest{1} '_' DESTE.timest{2} '_' DESTE.timest{3} '.mat'];
end

save(outfilename, '-struct', 'DESTE');

end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [signallevel, noiselevel]=estimate_noiselevel(incube)
%estimate the noise using NxN patches of the first and last slice
N=8;
refplane=incube(:,:,[1 : end]);
si=size(refplane);
ni=floor(si(1)/N);
nj=floor(si(2)/N);
nk=si(3); %for each slice
noisetest=zeros(ni,nj,nk);
amplitudetest=zeros(ni,nj,nk);
for kk=1:nk;
	for ii=0:ni-1;
		for jj=0:nj-1;
			ivec=(1:N) + ii*N;
			jvec=(1:N) + jj*N;
			kvec=kk;
			patch=refplane(ivec,jvec,kvec);
			noisetest(ii+1,jj+1,kk)=sqrt(var(abs(patch(:))));
			amplitudetest(ii+1,jj+1,kk)=mean(abs(patch(:)));
		end
	end
end

sortednoiselevel=sort(noisetest(:));
sortedamplitudelevel=sort(amplitudetest(:));

%the noise level is calculated where the signal is zero, or finite for abs
%value noise;
noiselevel=sqrt(2)*mean(sortednoiselevel(1:round(numel(sortednoiselevel)/2))) + ...
    mean(sortedamplitudelevel(1:round(numel(sortedamplitudelevel)/2)));
%signal level is calculate at the end of the signal distribution
signallevel=mean(sortedamplitudelevel(round(end-300):(end-10)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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











