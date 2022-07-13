function output=lig_deste(varargin)
% displacement-encoding-stimulated echoes
% out=deste(timestamp,[timestamp,timestamp..., ['grid', vector of grid size or matlab workspace containing datasize])
% script to read the amplitude and phase of stimulated echo images
% for the ligament/tendon project


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


% check if the data should be re-gridded
if any(strcmp(varargin,'grid'));
	ind=find(strcmp(varargin,'grid'));
	gridinfo=varargin{ind+1};
	%
	if numel(gridinfo) ==1 && numel(num2str(dummy)) ==4;
		%a time stamp for a .mat file has been handed down, from which the
		%grid size is read
		load(['analyzed_' num2str(gridinfo) '_T2star.mat'],'datasize');
		gridvec=datasize;
	elseif numel(gridinfo) == 3;
		%a three element size vector has been handed down
		gridvec=varargin{ind+1};
	else
		error('The gridding information must either contain the time stamp of an analyzed_xxxx_T2star.mat workspace, or a size vector with 3 elements!')
	end
end

clear dummy;


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
	out=varianms(fid_directory, 'noiso','flipfb'); %straight up reconstruction from fid
	
	phasefactor=out.image(:,:,:,1)./out.image(:,:,:,2);
	phasefactor=phasefactor./abs(phasefactor);
	
	compleximage=mean(abs(out.image),4).*phasefactor;
	
	blurflag=true;
	if blurflag;
		if ~any(strcmp(varargin,'grid'));
			if any(strcmp(varargin,'double'));
				gridvec=2*size(compleximage);
			else
				gridvec=size(compleximage);
			end
			
		end
		display(['outputgrid = [' num2str(gridvec(1)) ' ' num2str(gridvec(2)) ' ' num2str(gridvec(3)) '];']);
		%blur 3d
		dummy.image=compleximage;
		dummy.pars=out.pars;
		dummy=lig_blur3d(dummy,'vox',0.6*[1 1 1],'grid',gridvec);
		compleximage=dummy.image;
	end
	
	%get the paramters out of the fid directory
	pars=getparameters([fid_directory '/procpar']);
	
	%identify which of the three displacement encoding directions was
	%active, read, phase or slice, generate a filename for the output file
	dd(1)=(pars.dro(1)-pars.dro(2));
	dd(2)=(pars.dpe(1)-pars.dpe(2));
	dd(3)=(pars.dsl(1)-pars.dsl(2));
	directionindicator = {'RO'; 'PE'; 'SL'};
	lambda=zeros(1,3);
	for jj= 1:3;
		if dd(jj)~= 0
			lambda(jj)=round(10*pars.lambda/dd(jj) *1000);
			outfilename=['analyzed_' num2str(timestamp) '_' directionindicator{jj} '_lambda_' num2str(lambda(jj)) 'mu.mat'];
		end
	end
	
	%construct the complex image
	output.image=compleximage;
	si=size(compleximage);
	lambda_direction=find(dd~=0);
	
	%construct the spatial derivative image in the direction of applied
	%lambda
	si=size(compleximage);
	phasedifference=zeros(si);
	if lambda_direction ==1;
		for sl=1:si(3);
			phasedifference(2:si(1)-1,2:si(2)-1,sl)=compleximage(1:si(1)-2,2:si(2)-1,sl)./compleximage(2:si(1)-1,2:si(2)-1,sl);
		end
	end
	
	output.phasedifference=phasedifference;
	
	
	
	
	%extract the required parameters for plotting, construct the axes
	%all lengths in cm
	[numro numpe numsl]=size(output.image);
	lsl=pars.thk*pars.ns; %the length in the slice direction
	output.pars.xaxis=[0:numpe-1]/numpe*pars.lpe;
	output.pars.yaxis=[0:numro-1]/numro*pars.lro;
	output.pars.zaxis=[0:numsl-1]/numsl*lsl;
	output.pars.lambda=lambda;
	
	%claculate the diffusive attenuation factor b=k^2 Delta
	%standard units are s/mm^2
	output.pars.diffusive_b=(2*pi/(pars.lambda*10))^2*pars.TS;
	if any(strcmp(varargin,'double'));
		outfilename=['dblsz_' outfilename];
	end
	save(outfilename,'output');
	
	if any(strcmp(varargin,'plot'));
		ind=find(strcmp(varargin,'plot'));
		if nargin>ind && isnumeric(varargin{ind+1});
			slvec=varargin{ind+1};
			
			pic=compleximage;
			
			%plot the slices
			si = size(pic);
			figure('position' , [50 50 1200 600],'Name', outfilename);
			ns=si(3);
			
			jj=0;
			for sl=slvec;
				jj=jj+1;
				
				%make a mask for noiselevel
				
				subplot(3,numel(slvec),jj);
				imagesc(output.pars.xaxis,output.pars.yaxis,log10(abs(pic(:,:,sl))));
				set(gca,'Xticklabel','','Yticklabel','','clim',[2 4]);
				axis image;
				title(num2str(sl));
				
				subplot(3,numel(slvec),numel(slvec)+jj);
				imagesc(output.pars.xaxis,output.pars.yaxis,angle(pic(:,:,sl)));
				set(gca,'Xticklabel','','Yticklabel','');
				axis image;
				title(num2str(sl));
				
				subplot(3,numel(slvec),2*numel(slvec)+jj);
				imagesc(output.pars.xaxis,output.pars.yaxis,angle(phasedifference(:,:,sl)),[-1 0]);
				set(gca,'Xticklabel','','Yticklabel','');
				axis image;
				title(num2str(sl));
				
				
			end
		else
			error('Specify a vector of slice indices after the plot comman');
		end
		
	end
	
	
end



