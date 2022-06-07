
function tcourse = xtractROI(hdr, data, x, y, z)
%
% function tcourse = xtractROI(hdr, data, x, y, z)
%
% here the data is a 2D matrix with all the time points in rows
% x,y,z are vectors defining an arbitratry ROI
%
% (c) 2005 Luis Hernandez-Garcia 
% University of Michigan
% report bugs to:  hernan@umich.edu
%
%

num = length(x);


tcourse=zeros(size(data,1),1);
% for count=1:num
% 
%     tcourse = tcourse + ...
%         data( :,(z(count)-1)*(hdr.xdim*hdr.ydim) + (y(count)-1)*hdr.xdim + x(count)-1);
% end
% 
%tcourse = tcourse/num;
inds = sub2ind([hdr.xdim, hdr.ydim, hdr.zdim], x , y, z);
raw = data(:,inds);
for count=1:length(tcourse)
    tmp = raw(count,:);
    tcourse(count) = mean(tmp(isfinite(tmp)));
end

fprintf('... %d voxels in ROI ', num)


return
