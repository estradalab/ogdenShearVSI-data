function result = makslice()

% Luis Hernandez
% last edit : 1-8-97
%
% This function defines the slice structures as 
% described below
%
%result = struct(...
%'xdim',64,...      % number of pixels in the x direction
%'ydim',64,...      % number of pixels in the y direction
%'slicenum', 1,...  % slice number within the volume
%'n_slices',1,...   % Number of slices in the volume
%'data',zeros(1, 1));% data proper
%

result = struct(...%
'xdim',64,...      % number of pixels in the x direction
'ydim',64,...      % number of pixels in the y direction
'slicenum', 1,...  % slice number within the volume
'n_slices',1,...   % Number of slices in the volume
'pix_mm',4.735, ...% mm per pixel
'slice_mm',5, ...  % mm per slice
'data',zeros(1, 1));% data proper
  
return