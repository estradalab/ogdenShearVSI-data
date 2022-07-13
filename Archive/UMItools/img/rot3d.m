function odata = rot3d(indata, mtx)
% function odata = rot3d(indata, mtx)
%
% rotate a 3D data matrix by a transformation matrix and 
% resample it
%

xdim=size(indata,1);
ydim=size(indata,2);
zdim=size(indata,3);

[x1, y1, z1] = meshgrid(1:xdim, 1:ydim, 1:zdim);

x1 = reshape(x1,1, xdim*ydim*zdim);
y1 =  reshape(y1,1, xdim*ydim*zdim);
z1 = reshape(z1,1, xdim*ydim*zdim);

xyz1 = [ x1; y1; z1; ones(1,xdim*ydim*zdim)];

xyz2 = mtx * xyz1;
xyz2 = xyz2(1:3,:);


x2 = xyz2(1,:);
y2 = xyz2(2,:);
z2 = xyz2(3,:);

x2 = reshape(x2, xdim,ydim,zdim);
y2 =  reshape(y2, xdim,ydim,zdim);
z2 = reshape(z2, xdim,ydim,zdim);


[x1, y1, z1] = meshgrid(1:xdim, 1:ydim, 1:zdim);

odata = interp3(x2,y2,z2,indata, x1,y1,z1);

return