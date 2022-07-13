function [vx, vy, vz] = mm2vox( hdr, mxyz )
%function [vx, vy, vz] = mm2vox( hdr, [mx,my,mx])
%
% convert from mm to voxels using the AVW
% header information (including the origin)
% in AVW the origin is in voxels
%
% (c) 2005 Luis Hernandez-Garcia 
% University of Michigan
% report bugs to:  hernan@umich.edu
%
%
mx = mxyz(:,1);
my = mxyz(:,2);
mz = mxyz(:,3);

vx = round((mx )/hdr.xsize)+ hdr.origin(1);
vy = round((my )/hdr.ysize)+ hdr.origin(2);
vz = round((mz)/hdr.zsize) + hdr.origin(3);

return
