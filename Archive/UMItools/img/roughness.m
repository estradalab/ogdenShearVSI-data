function rgh = roughness(vol, mask)
% function rgh = roughness(vol)
%
% computes the norm of the gradient over a volume
% that serves as a decent measure of the roughness of an image
%

[dx, dy, dz] = gradient(vol);

if nargin==2
	mask=shrink(mask);
	dx = dx.*mask;
	dy = dy.*mask;
	dz = dz.*mask;
end

rgh = sqrt(dx.^2 + dy.^2 + dz.^2);
rgh = sum(sum(sum(rgh))) / length(rgh(:));


return
