function [mx, my, mz] = vox2mm( hdr, vxyz)
% function [mx, my, mz] = vox2mm( hdr, [vx, vy, vz])
%
% convert from voxels to mm using the AVW header information
% accounting for the origin.  in AVW, the origin is in voxels?
%
%
% (c) 2005 Luis Hernandez-Garcia 
% University of Michigan
% report bugs to:  hernan@umich.edu
%

vx = vxyz(:,1);
vy = vxyz(:,2);
vz = vxyz(:,3);

mx = (vx - hdr.origin(1)) * hdr.xsize;
my = (vy - hdr.origin(2)) * hdr.ysize;
mz = (vz - hdr.origin(3)) * hdr.zsize;

return