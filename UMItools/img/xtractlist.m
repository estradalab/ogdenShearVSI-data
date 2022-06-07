
function tcourse = xtractlist(hdr, data, list)
%
% function tcourse = xtractlist(hdr, data, list)
%
%   UNTESTED
%
% extracts a single time series from a list of voxels 
% time series of images
%
% hdr - header structure of one of the images in the time series.
% data - a 2D matrix where each row is an image and each column is a time
%       series
% list - a list of 3D coordinates in rows:  [x1 y1 z1 ; x2 y2 z2 ; ...]
%


fprintf('calculating the appropriate voxel coordinates...\n')
inds = sub2ind([hdr.xdim hdr.ydim hdr.zdim],list(:,1),list(:,2), list(:,3));
num=0;

tcourse=zeros(size(data,1),1);

for n=1:size(list,1)
    tcourse = tcourse + ...
        data( :,inds(n) );
    num = num +1;
end

tcourse = tcourse/num;

return
