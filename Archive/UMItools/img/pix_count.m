function pix_count(data)
%
% This function will count the voxels on the right and on the left
% of a plane selected graphically by the user. 
%
[x y] = ginput(2)
x = fix(x)
y = fix(y)

line( x,y);

midplane = [...
   x(1) y(1) 1;
   x(2) y(2) 1;
   x(2) y(2) 2]

rl = rightleft2( data, midplane);

return




