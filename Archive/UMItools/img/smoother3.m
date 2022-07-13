function smoother3(fname , sz)
% function smoother (fname , sz)
%
% uses Matlab's smooth3 function on a time series
%  sz is the number of voxels in the kernel
%
% this is the same as smoother2 but in NIFTI format
%

[p r ext]=fileparts(fname);

if ext=='.nii'
    isNIFTI=1;
    [data h] = read_nii_img(fname);
    h = nii2avw_hdr(h);
else
    isNIFTI=0;
    [data h] = read_img(fname);
end

odata = zeros(size(data));
slsz = sz;

data = reshape(data, h.tdim, h.xdim*h.ydim*h.zdim);

for t=1:size(data,1)
    in = reshape(data(t,:), h.xdim, h.ydim, h.zdim) ;
    if h.zdim==1
        in(:,:,2)=in;
        in(:,:,3)=in(:,:,1);
        slsz=1;
    end
    tmp = smooth3(in, 'gaussian', [ sz sz slsz ], 1.4);
    if h.zdim==1
        tmp = tmp(:,:,1);
    end
    odata(t,:) = tmp(:);
end



if isNIFTI
    h = avw2nii_hdr(h);
    write_nii(fullfile(p,['s' r '.nii']), odata, h, 0);
else 
    write_img(fullfile(p,['s' r ext]), odata, h);
end