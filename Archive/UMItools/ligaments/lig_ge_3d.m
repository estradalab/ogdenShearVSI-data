function out=lig_ge_3d(varargin)
%takes a mge3d data set of a ligament, data id bu time stamp
%options: 'plot' to make a movie
timestamp=varargin{1};


%% 1. read the data from FID
%data=varianms(timestamp,'apod',[0.4 0.5 0.5]);
data=varianms(timestamp,'pixvec',[192 64 64],'gradlag',48);  %approximately (0.25mm}^3 for typical ligament acquisition

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

%%make a mask
%blur3d is handed a structure, so the voxvec is the three principal axes of the blur volume in mm
blurvoxmm=[0.4 0.4 0.4];
blurreddata=blur3d(data,'vox',blurvoxmm);
blurredpic=blurreddata.image;

[blurredsignallevel,blurrednoiselevel] = estimate_noiselevel(blurredpic); 

mask=abs(blurredpic)>2.5*blurrednoiselevel;
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

save(matfilename, 'axis1','axis2','axis3','magnitude','mask','noiselevel','comment','timest','position','params');

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

