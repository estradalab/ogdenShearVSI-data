function M = makeAffine( rx, ry, rz, tx, ty,tz);
% function M = makeAffine( rx, ry, rz, tx, ty, tz);
% 
% Make a transformation matrix from movement parameters
% from here
% http://bishopw.loni.ucla.edu/AIR5/rigidbody.html
% M = eye(4);
% M(1,1) = cos(rz)*cos(ry)+sin(rx)*sin(ry)*sin(rz);
% M(1,2) = sin(rz)*cos(ry)-cos(rz)*sin(rx)*sin(ry);
% M(1,3) = cos(rx)*sin(ry);
% M(1,4) = tx;
% M(2,1) =  -sin(rz)*cos(rx);
% M(2,2) = cos(rx)*cos(rz);
% M(2,3) = sin(rx);
% M(2,4) = ty;
% M(3,1) = sin(rz)*sin(rx)*cos(ry) - cos(rz)*sin(ry);
% M(3,2) = -cos(rz)*sin(rx)*cos(ry) - sin(rz)*sin(ry);
% M(3,3) = cos(rx)*cos(ry);
% M(3,4) = tz;
%   
if nargin==1
    p = rx;
    rx = p(1); ry = p(2); rz = p(3); tx = p(4); ty=p(5);tz=p(6);
end

M = eye(4);
M(1,1) = cos(rz)*cos(ry)+sin(rx)*sin(ry)*sin(rz);
M(1,2) = sin(rz)*cos(ry)-cos(rz)*sin(rx)*sin(ry);
M(1,3) = cos(rx)*sin(ry);
M(1,4) = tx;
M(2,1) =  -sin(rz)*cos(rx);
M(2,2) = cos(rx)*cos(rz);
M(2,3) = sin(rx);
M(2,4) = ty;
M(3,1) = sin(rz)*sin(rx)*cos(ry) - cos(rz)*sin(ry);
M(3,2) = -cos(rz)*sin(rx)*cos(ry) - sin(rz)*sin(ry);
M(3,3) = cos(rx)*cos(ry);
M(3,4) = tz;

return