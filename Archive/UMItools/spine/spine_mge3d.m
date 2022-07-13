function mge3d=spine_mge3d(varargin)
% out=spine_mge3d(intimes)
%intimes is either a 4 digit time string, or a cell array of 4-digit time strings


%% default values
pix_value=1;   %zerofill to get a minimum pixel number of 192
vox_value=[0.075 0.075 0.25]; %size of the virtual voxel

inputwords={'plot','startslice','vox','pix','fliplr','flipud','flipfb','register','shphase'}; %optional input words 


parseinput(inputwords,varargin); %creates variables logical xx_flag, numeric xx_value, where xx is an input word
                                  

%%
jj=1;
loopend=false;
for kk=1:nargin;
    if (~any(strcmp(inputwords,varargin{kk})) && ~loopend);
        timelist{jj}=varargin{kk};
        jj=jj+1;
    else
        loopend=true;
    end
end

%input checking that all the listed times exist
d=dir;
names={d(:).name}; %cell array of file and directory names
for jj=1:numel(timelist);
    existflag=false;
    for kk=1:numel(names);
        if ~isempty(strfind(names{kk},timelist{jj}));
            existflag=true;
        end
    end
    if ~existflag;
        display(['No directory for time ' timelist{jj}]);
    end
end


%% assemple a fliplist
fliplist={[]};
if flipud_flag;
    fliplist=[fliplist,{'flipud'}];
end
if fliplr_flag;
    fliplist=[fliplist,{'fliplr'}];
end
if flipfb_flag;
    fliplist=[fliplist,{'flipfb'}];
end
if shphase_flag;
    fliplist=[fliplist,{'shphase'}, {shphase_value}];
end

%%
for jj=1:numel(timelist);
    data(jj)=varianms(timelist{jj}, num2str(pix_value),'zf3'); %doubles resolution in the 2nd PE dimension
    %data(jj)=varianms(timelist{jj}, num2str(pix_value),fliplist{:}); %doubles resolution in the 2nd PE dimension
end

if register_flag;
    data=imagereg(data);
end

%average the data;
mge3d=data(1);
mge3d.image=zeros(size(data(1).image));
for jj=1:numel(data);
    mge3d.image=mge3d.image+data(jj).image;
end

%% check if the the slices need to be re-ordered;
if startslice_flag;
    si=size(mge3d.image);
    mge3d.image=mge3d.image(:,:,[startslice_value:si(3) 1:startslice_value-1],:);
end

%% check if the voxel size has been changed, and blur 3d
%mge3d=blur3d(mge3d,'vox',vox_value);

%% optional plot
if plot_flag;
    summedGE=mge3d;
    summedGE.image=sum(mge3d.image,4);
    plotms(summedGE,'clim', [0 2e5]);
end
