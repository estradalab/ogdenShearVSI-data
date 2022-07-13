function [data, hdr] = read_img(h, name)

%	[data, hdr] = read_img(h, name)
%   or
%       [data, hdr] = read_img(name)
% 
% 
% Luis hernandez
% last edit 10-24-2008
%
% Loads the data from an analyze format file 'name' containing mutislice image data
% if the file contains multiple TIME points,
% this function returns a two dimensional array of data 
% (each row is an image , each column is a time series).
%
% NIFTI support added: If the suffix is .nii or .nii.gz , it calls read_nii_img
% 
% (c) 2006 Luis Hernandez-Garcia 
% University of Michigan
% report bugs to:  hernan@umich.edu
%
global SPM_scale_factor endian
data=[];


if nargin==1
    % if there is only one argument, that means that it is
    % the file name
    name = h;    
end

% what type of file is it ?
[pth fname suffix vrs] = fileparts(name);

if isempty(suffix)
    if exist([name '.img'], 'file')
        suffix = '.img';
		hname = fullfile(pth, [fname '.hdr']);
    end
    if exist([name '.nii'], 'file')
        suffix = '.nii';
    end
    if exist([name '.gz'], 'file')
        suffix = '.gz';
    end
end

name = fullfile(pth,[fname suffix]);

% now read it appropriately
switch(suffix)
    case '.img'
		hname = fullfile(pth, [fname '.hdr']);
        hdr = read_hdr(hname);
    case '.hdr'
		hname = fullfile(pth, [fname '.hdr']);
        hdr = read_hdr(hname);
 
	case '.gz'
    % these two cases are for NIFTI files
        [data, hdr] =read_nii_img(name);
		hdr = nii2avw_hdr(hdr);
        return
    case '.nii'
        [data, hdr]=read_nii_img(name);
		hdr = nii2avw_hdr(hdr);
        return
    otherwise
        hdr = read_hdr(sprintf('%s.hdr',fname));
        name = sprintf('%s.img',fname);
end


if (abs(SPM_scale_factor < eps)) 
    SPM_scale_factor=1;
end
%SPM_scale_factor


[pFile,messg] = fopen(name, 'r', endian);
if pFile == -1
    fprintf('%s - could not open %s \n', messg, name);
    return;
end

xdim = hdr.xdim;
ydim = hdr.ydim;
zdim = hdr.zdim;
tdim = hdr.tdim;

switch hdr.datatype     
    case 0
        fmt = 'uint8';
    case 2
        fmt = 'uint8';
    case 4
        fmt = 'short';
    case 8
        fmt = 'int';
    case 16
        fmt = 'float';
    case 32
        fmt = 'float';
        xdim = hdr.xdim * 2;
        ydim = hdr.ydim * 2;
    case 64
        fmt = 'int64';    
    otherwise
        fprintf('\nData Type %d Unsupported. Aborting\n',hdr.datatype);
        return
        
end



% Read in data.
%     for t=1:tdim 
%         fprintf('\r Reading time frame ... %d', t);
%     	d = (fread(pFile,[xdim*ydim*zdim], fmt))'; 	
%         data = [data ; d];
%     end
d = (fread(pFile,[xdim*ydim*zdim*tdim], fmt))'; 

% saving memory
%if hdr.datatype==4 || hdr.datatype==8
%    d=int16(d);
%end

if tdim >=2
    d = reshape(d, xdim*ydim*zdim, tdim);
    d=d';
%else
%    d = reshape(d, [xdim ydim zdim]); 
end
fclose(pFile);
data = d * SPM_scale_factor;
%SPM_scale_factor
%whos data
return





