function result = readvol(name)
%
% function result = readvol(varargin)
% Luis hernandez
% last edit 6-1-2001
%
%
% The function returns a 3D array 
% of dimensions determined by the header
%      dim x dim x no_slices


% Let user select filename ...

if nargin==0
   [file path] = uigetfile('*.img','Select Analyze file');
   name = strcat(path,file);
end
    
   sz = size(name);
   hdrname = strcat(name(1:(sz(2) - 4)), '.hdr');
   
   h = read_hdr(hdrname);
   result = read_img2(h, name);

   
return   
