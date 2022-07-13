function out=lig_sems(varargin)
%takes a mge3d data set of a ligament, data id bu time stamp
%options: 'plot' to make a movie
timestamp=varargin{1};

ind=find(strcmp(varargin,'gradlag'));
if ~isempty(ind)
	if nargin>ind
		gradlag=varargin{ind+1};
	else
		gradlag=30;
	end
else
	gradlag=40;
end


%% 1a. read the data from FID
%data=varianms(timestamp,'apod',[0.4 0.5 0.5]);
pixvec=[192 64 64];
data=varianms(timestamp,'pixvec',pixvec,'gradlag',gradlag);  %approximately (0.25mm}^3 for typical ligament acquisition

% 1b. read the data at the max resolution, which may be higher than [192 64 64]
siacq=[data.pars.np/2 data.pars.nv data.pars.ns];
orig_data=varianms(timestamp,'pixvec',max(siacq,pixvec),'gradlag',gradlag);

%% 2. apply a 3d blur to get the voxel size right
% blurring the structure, so voxel size in in mm
if any(strcmp(varargin,'vox'));
    ind=find(strcmp(varargin,'vox'));
    vxs=varargin{ind+1};
    if numel(vxs)~=3;
        vox=vxs(1)*[1 1 1];
    else
        vox=vxs(1:3);
    end
    data=blur3d(data,'vox',vox);
end

pic=data.image;

%% 3. calculate the time axis
%taxis=data.pars.te+(0:data.pars.ne-1)*data.pars.te2;


%% 4. pick the region of interest where the relaxation times and amplitudes
%are to be calculated use the central fifth of the readout direction
si=size(pic);
axis1=data.pars.axis1;  %mm
axis2=data.pars.axis2;  %mm
axis3=data.pars.axis3;  %mm

%%
%outf=matrix_expfit(pic,taxis);

magnitude=abs(pic);

%% make a mask based on the [192 64 64] data set
%blur3d is handed a structure, so the voxvec is the three principal axes of the blur volume in mm
blurvoxmm=[0.4 0.4 0.4];
blurreddata=blur3d(data,'vox',blurvoxmm);
blurredpic=blurreddata.image;
[blurredsignallevel,blurrednoiselevel] = estimate_noiselevel(blurredpic); 
mask=abs(blurredpic)>1/2*(blurredsignallevel-blurrednoiselevel);

% make another mask based on the orig_data dataset
%blur this one on the true resolution
sik=size(data.kspace);
origvox=abs([orig_data.pars.lro*10/sik(1), ...
    orig_data.pars.lpe*10/sik(2), ...
    orig_data.pars.thk]);
absmaskflag=true;
if absmaskflag
    temp=orig_data;
    temp.image=abs(temp.image);
    blurreddata=blur3d(temp,'vox',origvox);
    blurredpic=blurreddata.image;
    [blurredsignallevel,blurrednoiselevel] = estimate_noiselevel(blurredpic);
    orig_data.mask=abs(blurredpic)>1/2*(blurredsignallevel-blurrednoiselevel);
else
    blurreddata=blur3d(orig_data,'vox',1.4*origvox);
    blurredpic=blurreddata.image;
    [blurredsignallevel,blurrednoiselevel] = estimate_noiselevel(blurredpic);
    orig_data.mask=abs(blurredpic)>1/2*(blurredsignallevel-blurrednoiselevel);
end
    


comment=data.pars.comment;
timest = timestamp;

%check if a position was specified in the comment
comment=data.pars.comment;
if ~isempty(strfind(comment,'pos'));
    numbers_index = isstrprop(comment,'digit') | isstrprop(comment,'punct');
    totherightofpos=false(1,numel(numbers_index));
    posindex=strfind(comment,'pos');
    totherightofpos((posindex+3):end)=true;
    numbers_index=numbers_index & totherightofpos;
    position=str2num(comment(numbers_index));
else
    position=[];
end

params=data.pars;

%%
if isnumeric(timestamp);
    matfilename=['HIRES_' num2str(timestamp) '_magnitude.mat'];
else
    matfilename=['HIRES_' timestamp '_magnitude.mat'];
end

[signallevel,noiselevel] = estimate_noiselevel(pic); %#ok<ASGLU>

save(matfilename, 'axis1','axis2','axis3','magnitude','mask','noiselevel','comment','timest','position','params','orig_data');

out.axis1=axis1(:);
out.axis2=axis2(:);
out.axis3=axis3(:);
out.magnitude=magnitude;
out.signallevel=signallevel;
out.noiselevel=noiselevel;
out.mask=mask;
out.maskvoxelsize_mm=blurvoxmm;
out.masknoiselevel=blurrednoiselevel;
out.comment=comment;
out.timest=timest;
out.position=position;
out.params=params;
out.data=pic;
out.orig_data=orig_data;


%%
if any(strcmp(varargin,'plot'));
    figure('Position',[ 200 200 1100 400]);
    for jj=1:numel(axis1);
        
        
        subplot(1,2,1);
        MM=log10(squeeze(magnitude(jj,:,:))');
        imagesc(axis3,axis2,MM,[log10(2*out.noiselevel) log10(signallevel)+1]);
        colorbar;
        title(['magnitude, z=' num2str(axis1(jj),'%0.3f') 'mm']);
        axis image;
        xlabel('x (mm)');
        ylabel('y (mm)');
        
        subplot(1,2,2);
        mm=squeeze(mask(jj,:,:))';
        imagesc(axis3,axis2,mm,[-1 1]);
        title('mask');
        axis image;
        xlabel('x (mm)');
        ylabel('y (mm)');
        
        
        colormap jet;
        pause(0.1);
        
        gifflag=false;
        
        if gifflag
            %make the gif;
            frame=getframe(gcf);
            im=frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            
            % Write to the GIF File
            if jj == 1
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append');
            end
        end
    end
    
end;

