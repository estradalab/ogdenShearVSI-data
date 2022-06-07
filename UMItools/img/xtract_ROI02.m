function tseries = xtract_ROI02(rootname, xyz)
% function tseries = xtract_ROI02(rootname, xyz)
% extracts time courses from a time series of AVW images
% this is a version for large files that overwhelm the memory
% reads files on at a time.

hfiles = dir(sprintf('%s*.hdr', rootname));
imgfiles = dir(sprintf('%s*.img', rootname));

Nframes = length(imgfiles);
tseries = zeros(Nframes,1);

hdr = read_hdr(hfiles(1).name);
xdim = hdr.xdim; ydim = hdr.ydim; zdim = hdr.zdim;
inds = sub2ind([xdim ydim zdim], xyz(:,1), xyz(:,2), xyz(:,3) );

for t=1:Nframes
	tmp = read_img(imgfiles(t).name);
	tseries(t) = mean( tmp(inds) );
end

return
    

