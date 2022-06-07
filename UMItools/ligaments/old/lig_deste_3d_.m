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
		if isnumeric(dummy) && numel(num2str(dummy)) ==4;
			%this is a time stamp
			timestampvec{jj}=num2str(dummy);
			jj=jj+1;
		else
			istimestampflag =false;
		end
	else
		istimestampflag=false;
	end
end


% check if T2star data was specified as well
if numel(timestampvec) >3;
	
	jj=0;
	%the remaining timestamps must be mge3d t2*data
	for nt=4:numel(timestampvec);
		jj=jj+1;
		gridinfo=timestampvec{nt};
		T2starmatfilename=['HIRES_' num2str(gridinfo) '_T2star.mat'];
		if exist(T2starmatfilename,'file');
			HIRES(jj)=load(T2starmatfilename);
		else
			HIRES(jj)=lig_T2star(gridinfo);
		end
	end
	outgridvec=size(HIRES(1).magnitude);
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
	
    out=varianms(fid_directory, 'noiso','flipfb'); %straight up reconstruction from fid
	
	
	
	%get and reconstruct the multislice image from the fid directory
	%the image is 4-d three spatial dimensions and plus minus lambda.
	
	orig_phasefactor(:,:,:,nf)=out.image(:,:,:,1)./out.image(:,:,:,2);
	orig_phasefactor(:,:,:,nf)=orig_phasefactor(:,:,:,nf)./abs(orig_phasefactor(:,:,:,nf));
	orig_compleximage(:,:,:,nf)=mean(abs(out.image),4).*orig_phasefactor(:,:,:,nf);  %averages over plus minus lambda
    
    
end

%Now only some of the slices in the 2-d acquisition actuall contain
%ligament, so before processing that further set the noise part of the
%volume to zero using maximum intensity projection

%dummy=sum(abs(orig_compleximage),4);
%projection1=sum(dummy,3);
%projection2=squeeze(sum(dummy,2));

%dummy2=HIRES(1).mask;
%pro1=sum(dummy2,3);
%pro2=squeeze(sum(dummy2,2));

for nf=1:3;
    

	
	si=size(orig_compleximage(:,:,:,nf));
	gridvec=[si(1) 2*si(2) 2*si(3)]; %zerofill in the PE and SL directions
	
	display(['zerofilled STE image to size [' num2str(gridvec(1)) ' ' num2str(gridvec(2)) ' ' num2str(gridvec(3)) '];']);
	%blur 3d
	voxelsize=0.2;
	display(['blur size=' num2str(voxelsize) 'mm']);
	dummy.image=orig_compleximage(:,:,:,nf);
	dummy.pars=out.pars;
	dummy=lig_blur3d(dummy,'vox',voxelsize*[1 1 1],'grid',gridvec);
	
	%update the complex image and phase factor to the new grid size
	compleximage(:,:,:,nf)=dummy.image;
	
	%get the paramters out of the fid directory
	pars=out.pars;
	
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
	wigglecomment{nf}=[timestamp ': ' out.pars.comment];
	
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
if true
     figure('position',[400         563        1097         406]);
    for jj=1:256;
        subplot(1,2,1);
        imagesc(DESTE.axis3, DESTE.axis2,log10(squeeze(absvalueimage(jj,:,:).*DESTE.mask(jj,:,:))+1),[0 5]);
        axis square;
        subplot(1,2,2);
        imagesc(HIRES.axis3,HIRES.axis2,log10(squeeze(HIRES.magnitude(jj,:,:))),[0 6]);
        axis square;
        pause(0.1);
    end
end
%%



%now calculate the numeric drivatives for the strain image
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

%% check if a wiggle was specified in the comment
for jj=1:3;
	comment=wigglecomment{jj};
	if ~isempty(strfind(comment,'wiggle'));
		numbers_index = isstrprop(comment,'digit') | isstrprop(comment,'punct');
		totherightofwiggle=false(1,numel(numbers_index));
		wiggleindex=strfind(comment,'wiggle');
		totherightofwiggle((wiggleindex+3):(wiggleindex+7))=true;
		numbers_index=numbers_index & totherightofwiggle;
		wiggle(jj)=str2num(comment(numbers_index));
	else
		wiggle(jj)=[];
	end
end

DESTE.wiggle=wiggle;

%check if a position was specified in the comment
for jj=1:3;
	comment=wigglecomment{jj};
	if ~isempty(strfind(comment,'pos'));
		numbers_index = isstrprop(comment,'digit') | isstrprop(comment,'punct');
		totherightofpos=false(1,numel(numbers_index));
		posindex=strfind(comment,'pos');
		totherightofpos((posindex+3):end)=true;
		numbers_index=numbers_index & totherightofpos;
		position{jj}=str2num(comment(numbers_index));
	else
		position{jj}=[];
	end
end

DESTE.position=position;
DESTE.pars=pars;


%%
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
if false;
%if exist('HIRES','var');
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
        
		
		figure;
		for jj=30:96;
			subplot(1,2,1);
			imagesc(DESTE.HIRES.axis2,DESTE.HIRES.axis1, log(hir_absvalueimage(:,:,jj)));
			
			subplot(1,2,2); imagesc(DESTE.HIRES.axis2,DESTE.HIRES.axis1, log(DESTE.HIRES.magnitude(:,:,jj)));
			set(gcf,'Name',num2str(jj));
			pause(0.5);
        end
		
	end
	%%
	for ii=1:3;
		for jj=1:3;
			hi_res_strain(:,:,:,ii,jj)=interp3(X2q,X1q,X3q,squeeze(DESTE.strains(:,:,:,ii,jj)),X2ref,X1ref,X3ref);
		end
    end		
	DESTE.hi_res_strain=hi_res_strain;
	%% strain_HIRES
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
noiselevel=sqrt(2)*mean(sortednoiselevel(1:round(numel(sortednoiselevel)/2)));
signallevel=mean(sortedamplitudelevel(round(end*0.9):end));
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











