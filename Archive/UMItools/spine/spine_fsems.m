function [fsems,data]=spine_fsems(varargin)
% out=spine_mge3d(intimes)
%intimes is either a 4 digit time string, or a cell array of 4-digit time strings

%% default values
pix_value=384;   %zerofill to get a minimum pixel number of 192

inputwords={'plot','nav','vox','pix','register'}; %optional input words 
parseinput(inputwords,varargin);



%% input checking
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

%%
for jj=1:numel(timelist);
    data(jj)=varianms(timelist{jj}, num2str(pix_value),'filter'); %doubles resolution in the 2nd PE dimension
end

if register_flag;
    data=imagereg(data);
end

%out=imagereg(data);


%average the data;
fsems=data(1);
fsems.image=zeros(size(data(1).image));
for jj=1:numel(data);
    fsems.image=fsems.image+abs(data(jj).image);
end

if plot_flag;
    plotms(fsems);
end
