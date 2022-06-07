function [data, hdr] = read_img2(hdr, name)

% [image_data, hdr] = read_img2([hdr,] name)
%
% Luis hernandez
% last edit 2-23-2007
%
% Loads the data from a file 'name' containing mutislice image data
% 		hdr: is a header structure (use read_hdr()).  Not necessary by,
% 		left there to make compatible with old code.
% 		name:  file name (.img)
% The function returns a three dimensional array 
% NIFTI support added: If the suffix is .nii or .nii.gz , it calls read_nii_img
% Note (2/26/08) - converts nifti headers to analyze headers
%
% returns only the first time point as a 3D matrix
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
    name = hdr;    
end

if length(name)> 3
   suffix = name(end-3:end);
   switch(suffix)
   case '.img'
       hdr = read_hdr(sprintf('%s.hdr',name(1:end-4)));
   case '.hdr'
       hdr = read_hdr(name);
       name = sprintf('%s.img',name(1:end-4));
   % these two cases are for NIFTI files
   case 'i.gz'
       [data, hdr] =read_nii_img(name);
	   hdr = nii2avw_hdr(hdr);
       return
   case '.nii'
       [data, hdr]=read_nii_img(name);
	   hdr = nii2avw_hdr(hdr);
       return
   otherwise
       hdr = read_hdr(sprintf('%s.hdr',name));
       name = sprintf('%s.img',name);
   end
else
    hdr = read_hdr(sprintf('%s.hdr',name));
    name = sprintf('%s.img',name);
end


if (abs(SPM_scale_factor < eps)) 
    SPM_scale_factor=1;
end
%SPM_scale_factor


[pFile,messg] = fopen(name, 'r', endian);
if pFile == -1
    fprintf('%s - could not open %s', messg, name);
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
        errormesg(sprintf('Data Type %d Unsupported. Aborting',hdr.datatype));
        return
        
end



% Read in data.
%     for t=1:tdim 
%         fprintf('\r Reading time frame ... %d', t);
%     	d = (fread(pFile,[xdim*ydim*zdim], fmt))'; 	
%         data = [data ; d];
%     end
d = (fread(pFile,[xdim*ydim*zdim*tdim], fmt))'; 	
if tdim >=2
    d = reshape(d, xdim*ydim*zdim, tdim);
    d=d';
else
   d = reshape(d, [xdim ydim zdim]); 
end
fclose(pFile);
data = d * SPM_scale_factor;
%SPM_scale_factor
%whos data
return

