function output=lig_deste_3d(varargin)
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
end

[dummy, phaseindvec ]=sort(encodedirectionindex);

%extract the required parameters for plotting, construct the axes
%all lengths in mm
zfis=size(compleximage{1});  % zero filled image size
zf_factor=zfis./original_imagesize;

axis1spacing=(out.pars.yaxis(2)-out.pars.yaxis(1))/zf_factor(1); %in mm
axis2spacing=(out.pars.xaxis(2)-out.pars.xaxis(1))/zf_factor(2); %in mm
axis3spacing=(out.pars.zaxis(2)-out.pars.zaxis(1))/zf_factor(3); %in mm

output.axis1=out.pars.yaxis(1)+(axis1spacing*((0:zfis(1)-1)-0.5)); %in mm
output.axis2=out.pars.xaxis(1)+(axis2spacing*((0:zfis(2)-1)-0.5));
output.axis3=out.pars.zaxis(1)+(axis3spacing*((0:zfis(3)-1)-0.5));
output.pars.lambda=lambda;




%create a pancake stack of absolute value images
absvalueimage=zeros(size(compleximage{1}));
for nf=1:numel(timestampvec);
	absvalueimage=absvalueimage+abs(compleximage{nf});
	phaseimage(:,:,:,nf)=phasefactor{phaseindvec(nf)};
end


%now calculate the numeric drivatives for the strain image
si=size(phaseimage);
strainimage=zeros(si);
ind1=[2:si(1)-1];
ind2=[2:si(2)-1];
ind3=[2:si(3)-1];

delta1=output.axis1(2)-output.axis1(1);
delta2=output.axis2(2)-output.axis2(1);
delta3=output.axis3(2)-output.axis3(1);

strainimage(ind1,:,:,1)=angle(phaseimage(ind1+1,:,:,1)./phaseimage(ind1-1,:,:,1))/(2*pi)*lambda(1)/1000/(2*delta1);
strainimage(:,ind2,:,2)=angle(phaseimage(:,ind2+1,:,2)./phaseimage(:,ind2-1,:,2))/(2*pi)*lambda(2)/1000/(2*delta2);
strainimage(:,:,ind3,3)=angle(phaseimage(:,:,ind3+1,3)./phaseimage(:,:,ind3-1,3))/(2*pi)*lambda(3)/1000/(2*delta3);


%if requested, regrid the data



%try making a mask
if any(strcmp(varargin,'mask'));
	[signallevel,noiselevel] = estimate_noiselevel(absvalueimage);
	SNR=signallevel/noiselevel;
	maskabs=absvalueimage<2*noiselevel;
	
	absvalueimage(maskabs)=NaN;
	for nf=1:numel(timestampvec);
		dummy=phaseimage(:,:,:,nf);
		dummy(maskabs)=NaN;
		phaseimage(:,:,:,nf)=dummy;
	end
end


si=size(compleximage{1});


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
	figure('position' , [50 50 50*nsl 600]);
	cmap=jet;
	cmap=[ [0 0 0]; cmap];
	set(gcf,'colormap',cmap);
	
	jj=0;
	for sl=slvec;
		jj=jj+1;
		
		%make a mask for noiselevel
		
		tsubplot(4,numel(slvec),jj,3);
		imagesc(output.axis2,output.axis1,log10(absvalueimage(:,:,sl)));
		set(gca,'Xticklabel','','Yticklabel','','clim',[2 4]);
		axis image;
		title(num2str(sl));
		
		tsubplot(4,numel(slvec),numel(slvec)+jj,3);
		imagesc(output.axis2,output.axis1,angle(phaseimage(:,:,sl,1)),[-pi pi]);
		set(gca,'Xticklabel','','Yticklabel','');
		axis image;
		title(num2str(sl));
		
		tsubplot(4,numel(slvec),2*numel(slvec)+jj,3);
		imagesc(output.axis2,output.axis1,angle(phaseimage(:,:,sl,2)),[-pi pi]);
		set(gca,'Xticklabel','','Yticklabel','');
		axis image;
		title(num2str(sl));
		
		tsubplot(4,numel(slvec),3*numel(slvec)+jj,3);
		imagesc(output.axis2,output.axis1,angle(phaseimage(:,:,sl,3)),[-pi pi]);
		set(gca,'Xticklabel','','Yticklabel','');
		axis image;
		title(num2str(sl));
		
		
	end
	
	
	%plot the slices of deirvatives
	strainrange=0.4;
	si = size(compleximage{1});
	figure('position' , [50 50 50*nsl 600],'Name',['Strain color = +-' num2str(strainrange)]);
	cmap=jet;
	cmap=[ [0 0 0]; cmap];
	set(gcf,'colormap',cmap);
	
	
	
	
	
	
	
	if any(strcmp(varargin,'grid'));
		absvalueimage(isnan(absvalueimage))=0;
		strainimage(isnan(strainimage))=0;
		phaseimage(isnan(phaseimage))=1;
		
		dum1.image=absvalueimage;
		dum2=lig_blur3d(dum1,'grid',outgridvec);
		output.amplitude=dum2.image;
		
		dum1.image=strainimage;
		dum2=lig_blur3d(dum1,'grid',outgridvec);
		output.strainimage=dum2.image;
		
		dum1.image=phaseimage;
		dum2=lig_blur3d(dum1,'grid',outgridvec);
		output.phaseimage=dum2.image;
	else
		output.amplitude=absvalueimage;
		output.strainimage=strainimage;
		output.phaseimage=phaseimage;
		
	end
	
	
	
	
	jj=0;
	for sl=slvec;
		jj=jj+1;
		
		tsubplot(4,numel(slvec),jj,3);
		imagesc(output.axis2,output.axis1,log10(absvalueimage(:,:,sl)));
		set(gca,'Xticklabel','','Yticklabel','','clim',[2 4]);
		axis image;
		title('amp');
		
		tsubplot(4,numel(slvec),numel(slvec)+jj,3);
		imagesc(output.axis2,output.axis1,(strainimage(:,:,sl,1)),[-1 1]*strainrange);
		set(gca,'Xticklabel','','Yticklabel','');
		axis image;
		title(['ud ' num2str(sl)]);
		
		tsubplot(4,numel(slvec),2*numel(slvec)+jj,3);
		imagesc(output.axis2,output.axis1,(strainimage(:,:,sl,2)),[-1 1]*strainrange);
		set(gca,'Xticklabel','','Yticklabel','');
		axis image;
		title(['lr ' num2str(sl)]);
		
		tsubplot(4,numel(slvec),3*numel(slvec)+jj,3);
		imagesc(output.axis2,output.axis1,(strainimage(:,:,sl,3)),[-1 1]*strainrange);
		set(gca,'Xticklabel','','Yticklabel','');
		axis image;
		title(['io ' num2str(sl)]);
		
		
	end
	
	
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











