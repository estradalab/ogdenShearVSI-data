
function tcourses = xtractROI02(hdr, data, x, y, z)
%
% function tcourses = xtractROI(hdr, data, x, y, z)
%
% here the data is a 2D matrix with all the time points in rows
% x,y,z are vectors defining an arbitratry ROI
%
% the output is several columns with a time course for every pixel in the
% ROI
%
% (c) 2010 Luis Hernandez-Garcia 
% University of Michigan
% report bugs to:  hernan@umich.edu
%
%

num = length(x);


tcourses=zeros(size(data,1),num);
for count=1:num

    tcourses(:,count) =  ...
        data( :,(z(count)-1)*(hdr.xdim*hdr.ydim) + (y(count)-1)*hdr.xdim + x(count)-1);
end

fprintf('... %d voxels in ROI ', num)



return
