function nii2avw(name)
% function nii2avw(name)
% 
% converts NIFTI images to AVW (Analyze) format images disregarding
% the affine transformation parameters
%


[dnii hnii] = read_nii_img(name);
avwh = nii2avw_hdr(hnii);

write_img([name '.img'],dnii,avwh);

return