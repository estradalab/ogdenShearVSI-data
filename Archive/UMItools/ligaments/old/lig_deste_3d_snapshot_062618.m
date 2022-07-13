function phase_struct=lig_deste_3d(varargin)
% displacement-encoding-stimulated echoes
% out=deste(timestamp,[timestamp,timestamp...[,options],
%'grid', [vector of grid size or matlab workspace containing datasize],
% 'mask',  %masks based on estimated noise level
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

%2b. check if there are three time stamps
if numel(timestampvec)~=3;
	display('Expecting three time stamps, one for RO, one for PE and one for SL slice encode.');
	complete_3dflag=false;
else
	complete_3dflag=true;
end


% check if the data should be re-gridded
if any(strcmp(varargin,'grid'));
	ind=find(strcmp(varargin,'grid'));
	gridinfo=varargin{ind+1};
	%
	if numel(gridinfo) ==1 && numel(num2str(dummy)) ==4;
		%a time stamp for a .mat file has been handed down, from which the
		%grid size is read
		load(['analyzed_' num2str(gridinfo) '_T2star.mat'],'datasize');
		outgridvec=datasize;
	elseif numel(gridinfo) == 3;
		%a three element size vector has been handed down
		outgridvec=varargin{ind+1};
	else
		error('The gridding information must either contain the time stamp of an analyzed_xxxx_T2star.mat workspace, or a size vector with 3 elements!')
	end
end

clear dummy;



lambda=zeros(1,numel(timestampvec));
encodedirectionvector={};
for nf=1:numel(timestampvec);
	
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
	%img_directory=d(time_flag & img_flag).name;
	%phi_directory=d(time_flag & PH_flag).name;
	
	
	%get and reconstruct the multislice image from the fid directory
	%the image is 4-d three spatial dimensions and plus minus lambda.
	
	display(['Reconstructing data from fid, timestamp =' num2str(timestamp)]);
	out=varianms(fid_directory, 'noiso','flipfb','shphase'); %straight up reconstruction from fid
	
	phasefactor{nf}=out.image(:,:,:,1)./out.image(:,:,:,2);
	phasefactor{nf}=phasefactor{nf}./abs(phasefactor{nf});
	
	compleximage{nf}=mean(abs(out.image),4).*phasefactor{nf};
	
	original_imagesize=size(compleximage{nf});
	
	blurflag=true;
	if blurflag;
		
		si=size(compleximage{nf});
		gridvec=[si(1) 2*si(2) 2*si(3)]; %zerofill in the PE and SL directions
		
		display(['zerofilled STE image to size [' num2str(gridvec(1)) ' ' num2str(gridvec(2)) ' ' num2str(gridvec(3)) '];']);
		%blur 3d
		voxelsize=0.6;
		display(['blur size=' num2str(voxelsize) 'mm']);
		dummy.image=compleximage{nf};
		dummy.pars=out.pars;
		dummy=lig_blur3d(dummy,'vox',voxelsize*[1 1 1],'grid',gridvec);
		
		%update the complex image and phase factor to the new grid size
		compleximage{nf}=dummy.image;
		phasefactor{nf}=compleximage{nf}./abs(compleximage{nf});
	end
	
	
	
	
	
	%get the paramters out of the fid directory
	pars=getparameters([fid_directory '/procpar']);
	
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
	
end

[dummy, phaseindvec ]=sort(encodedirectionindex);   %

%extract the required parameters for plotting, construct the axes
%all lengths in mm
zfis=size(compleximage{1});  % zero filled image size
zf_factor=zfis./original_imagesize;

output.axis1=newaxis(out.pars.axis1,zfis(1)); %in mm
output.axis2=newaxis(out.pars.axis2,zfis(2)); %in mm
output.axis3=newaxis(out.pars.axis3,zfis(3)); %in mm


% SORTING ****************************************************
%this is where the phase data gets sorted into RO, PE, SL
absvalueimage=zeros(size(compleximage{1}));
for nf=1:numel(timestampvec);
	absvalueimage=absvalueimage+abs(compleximage{nf});
	phaseimage(:,:,:,nf)=phasefactor{phaseindvec(nf)}; 
	dumdum1{nf}=wigglecomment{phaseindvec(nf)};
	dumdum2{nf}=timestampvec{phaseindvec(nf)};
end
lambda=lambda(phaseindvec); output.pars.lambda=lambda;
wigglecomment=dumdum1;
timestampvec=dumdum2;
% *********************************************************


%now calculate the numeric drivatives for the strain image
si=size(phaseimage);
complex_strainimage=zeros(si);
ind1=[2:si(1)-1];
ind2=[2:si(2)-1];
ind3=[2:si(3)-1];

% RO strain
complex_strainimage(ind1,:,:,1)=(phaseimage(ind1+1,:,:,1)./phaseimage(ind1-1,:,:,1));
% PE strain
complex_strainimage(:,ind2,:,2)=(phaseimage(:,ind2+1,:,2)./phaseimage(:,ind2-1,:,2));
% SL strain
complex_strainimage(:,:,ind3,3)=(phaseimage(:,:,ind3+1,3)./phaseimage(:,:,ind3-1,3));


%if requested, regrid the data
if any(strcmp(varargin,'grid'));
	dum1.image=absvalueimage;
	dum2=lig_blur3d(dum1,'grid',outgridvec);
	absvalueimage=abs(dum2.image);
	
	dum1.image=phaseimage;
	dum2=lig_blur3d(dum1,'grid',outgridvec);
	phaseimage=dum2.image;
	
	dum1.image=complex_strainimage;
	dum2=lig_blur3d(dum1,'grid',outgridvec);
	complex_strainimage=dum2.image;
	
	%redefine the axes again, due to re-gridding
end

delta1=output.axis1(2)-output.axis1(1);
delta2=output.axis2(2)-output.axis2(1);
delta3=output.axis3(2)-output.axis3(1);

strainimage=angle(complex_strainimage);
strainimage(:,:,:,1)=strainimage(:,:,:,1)/(2*pi)*lambda(1)/1000/(2*delta1);
strainimage(:,:,:,2)=strainimage(:,:,:,2)/(2*pi)*lambda(2)/1000/(2*delta2);
strainimage(:,:,:,3)=strainimage(:,:,:,3)/(2*pi)*lambda(3)/1000/(2*delta3);



%try making a mask
if any(strcmp(varargin,'mask'));
	
	ind=find(strcmp(varargin,'mask'));
	if nargin>ind && isnumeric(varargin{ind+1})
		gridinfo=varargin{ind+1};
		%
		if numel(gridinfo) ==1 && numel(num2str(gridinfo)) ==4;
			%a time stamp for a .mat file has been handed down, from which the
			%grid size is read
			d=dir;
			for jj=1:numel(dir);
				if ~isempty(strfind(d(jj).name,['analyzed_' num2str(gridinfo)]));
					loadfilename=d(jj).name;
				end
			end
			
			if ~exist('loadfilename','var');
				 lig_T2star(gridinfo);
				 loadfilename=['analyzed_' num2str(gridinfo) '_T2star.mat'];
				 pause(2);
			end
			
			maskreference=load(loadfilename);
			refdata=abs(blur3d(maskreference.magnitude,'vox', 0.6));
			[X1ref,X2ref,X3ref]=ndgrid(maskreference.axis1,maskreference.axis2,maskreference.axis3);
			
			display(['Masking on ' loadfilename]);
			
			[X1q,X2q,X3q]=ndgrid(output.axis1,output.axis2,output.axis3);
			
			Vq = interp3(X2ref,X1ref,X3ref,refdata,X2q,X1q,X3q);
			
			[signallevel,noiselevel] = estimate_noiselevel(absvalueimage);
			maskabs=1*(absvalueimage>2*noiselevel); %numeric mask, ones for good pixels, zeros for insufficient signal
			
			
		else
			error('The masking information must contain the time stamp of an analyzed_xxxx_T2star.mat workspace!')
		end
		%sample the mask on the grid of the phase data
	else
		display('Masking on abs value of STE');
		
		[signallevel,noiselevel] = estimate_noiselevel(absvalueimage);
		maskabs=1*(absvalueimage>2*noiselevel);  %numeric mask, ones for good pixels, zeros for insufficient signal
		
	end
	
else
	'Specify a mask giving a timestamp, or without argument.';
	return
end


si=size(absvalueimage);



output.amplitude=absvalueimage;
output.strainimage=strainimage;
output.phaseimage=phaseimage;
output.mask=maskabs;



if any(strcmp(varargin,'plot'));
	ind=find(strcmp(varargin,'plot'));
	if nargin>ind && isnumeric(varargin{ind+1});
		slvec=varargin{ind+1};
		
	else
		slvec=1:si(3);
	end
	nsl=numel(slvec);
	
	%plot the slices
	si = size(compleximage{1});
	figure('position' , [50 50 95*nsl 1000]);
	cmap=jet;
	cmap=[ [0 0 0]; cmap];
	set(gcf,'colormap',cmap);
	
	jj=0;
	for sl=slvec;
		jj=jj+1;
		
		%make a mask for noiselevel
		
		plotmask_amplitude=maskabs(:,:,sl);
		plotmask_angle=plotmask_amplitude;
		plotmask_angle(plotmask_amplitude==0)=NaN;
		
		tsubplot(4,numel(slvec),jj,3);
		imagesc(output.axis2,output.axis1,log10(absvalueimage(:,:,sl)).*plotmask_amplitude);
		set(gca,'Xticklabel','','Yticklabel','','clim',[2 4]);
		axis image;
		title(num2str(sl));
		
		tsubplot(4,numel(slvec),numel(slvec)+jj,3);
		imagesc(output.axis2,output.axis1,angle(phaseimage(:,:,sl,1)).*plotmask_angle,[-pi pi]);
		set(gca,'Xticklabel','','Yticklabel','');
		axis image;
		title(num2str(sl));
		
		tsubplot(4,numel(slvec),2*numel(slvec)+jj,3);
		imagesc(output.axis2,output.axis1,angle(phaseimage(:,:,sl,2)).*plotmask_angle,[-pi pi]);
		set(gca,'Xticklabel','','Yticklabel','');
		axis image;
		title(num2str(sl));
		
		tsubplot(4,numel(slvec),3*numel(slvec)+jj,3);
		imagesc(output.axis2,output.axis1,angle(phaseimage(:,:,sl,3)).*plotmask_angle,[-pi pi]);
		set(gca,'Xticklabel','','Yticklabel','');
		axis image;
		title(num2str(sl));
	end
	
	
	%plot the slices of deirvatives
	strainrange=0.5;
	si = size(compleximage{1});
	figure('position' , [50 50 95*nsl 1000],'Name',['Strain color = +-' num2str(strainrange)]);
	cmap=jet;
	cmap=[ [0 0 0]; cmap];
	set(gcf,'colormap',cmap);
	
	
	
	jj=0;
	for sl=slvec;
		jj=jj+1;
		
		plotmask_amplitude=maskabs(:,:,sl);
		plotmask_angle=plotmask_amplitude;
		plotmask_angle(plotmask_amplitude==0)=NaN;
		
		tsubplot(4,numel(slvec),jj,3);
		imagesc(output.axis2,output.axis1,log10(absvalueimage(:,:,sl)).*plotmask_amplitude);
		set(gca,'Xticklabel','','Yticklabel','','clim',[2 4]);
		axis image;
		title('amp');
		
		tsubplot(4,numel(slvec),numel(slvec)+jj,3);
		imagesc(output.axis2,output.axis1,(strainimage(:,:,sl,1)).*plotmask_angle,[-1 1]*strainrange);
		set(gca,'Xticklabel','','Yticklabel','');
		axis image;
		title(['ud ' num2str(sl)]);
		
		tsubplot(4,numel(slvec),2*numel(slvec)+jj,3);
		imagesc(output.axis2,output.axis1,(strainimage(:,:,sl,2)).*plotmask_angle,[-1 1]*strainrange);
		set(gca,'Xticklabel','','Yticklabel','');
		axis image;
		title(['lr ' num2str(sl)]);
		
		tsubplot(4,numel(slvec),3*numel(slvec)+jj,3);
		imagesc(output.axis2,output.axis1,(strainimage(:,:,sl,3)).*plotmask_angle,[-1 1]*strainrange);
		set(gca,'Xticklabel','','Yticklabel','');
		axis image;
		title(['io ' num2str(sl)]);
		
		
	end
	
	
end

phase_struct.axis1=output.axis1;
phase_struct.axis2=output.axis2;
phase_struct.axis3=output.axis3;
phase_struct.magimage=absvalueimage;
phase_struct.phaseimage=phaseimage;
phase_struct.mask=maskabs;
phase_struct.comment='Doubled resolution in dims 2&3 by zerofilling, mask via noise estimate of magimage.';
phase_struct.wigglecomments=wigglecomment;
phase_struct.timestamps=timestampvec;
phase_struct.true_lambda=lambda;


% if the mge3 T2 star reference data was used to calculate the mask,
% then regrid the strain data onto the grid of the T2 star data and output
% that as well

if exist('maskreference','var');
	
	[signallevel,noiselevel] = estimate_noiselevel(maskreference.magnitude);
	
	mask=1*(maskreference.magnitude>2*noiselevel);
	
	[X1m,X2m,X3m]=ndgrid(maskreference.axis1,maskreference.axis2,maskreference.axis3); %maskreference axes
	
	[X1s,X2s,X3s]=ndgrid(output.axis1,output.axis2,output.axis3); %strain data axes
	
	si=size(maskreference.magnitude);
	upsampled_strain=zeros([si 3]);
	
	for straindirection=1:3;
		
		%Vq = interp3(X2ref,X1ref,X3ref,refdata,X2q,X1q,X3q);
		
		upsampled_strain(:,:,:,straindirection)=interp3(X2s,X1s,X3s,strainimage(:,:,:,straindirection),X2m,X1m,X3m);
	end
	
	axis1=maskreference.axis1;
	axis2=maskreference.axis2;
	axis3=maskreference.axis3;
	magnitude=maskreference.magnitude;
	strain=upsampled_strain;
	timestamps=timestampvec;
	wigglecomments=wigglecomment;
	true_lambda=lambda;
	
	
	outfilename=['strain_and_magnitude' timestampvec{1} '_' timestampvec{2} '_' timestampvec{3} '_' num2str(gridinfo) '.mat'];
	
	save(outfilename, 'axis1','axis2','axis3','magnitude','strain','mask','phase_struct','timestamps', 'wigglecomments', 'true_lambda');
end



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











