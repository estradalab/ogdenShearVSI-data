function smoother2(root , sz)
% function smoother (root , sz)
%
% uses Matlab's smooth3 function on a time series
%  sz is the number of voxels in the kernel
%

[p root ext vr] = fileparts(root);
	
[data h] = read_img(root);
odata = zeros(size(data));

for t=1:h.tdim
	in = reshape(data(t,:), h.xdim, h.ydim, h.zdim);
	tmp = smooth3(in, 'gaussian', [ sz sz sz ]);
	odata(t,:) = tmp(:);
end

write_img(['s' root '.img'], odata, h);
