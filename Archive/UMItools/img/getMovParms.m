function parms = getMovParms (M)
%function [rotx roty rotz tx ty tz ] = getMovParms (M)
% cant' handle rotations greater than pi/2
%

% from http://bishopw.loni.ucla.edu/AIR5/rigidbody.html
rotx = asin(M(2,3));
roty = atan(M(1,3) / M(3,3));
rotz = -atan(M(2,1)/M(2,2));

tx = M(1,4);
ty = M(2,4);
tz = M(3,4);

parms = [rotx, roty, rotz, tx,ty,tz];

return
