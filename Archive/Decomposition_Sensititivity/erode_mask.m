function [new_mask] = erode_mask(mask)
bw_mask = mask==1;
se = strel('cube',1);
new_bw_mask = imerode(bw_mask,ones(3,3,3));
new_mask = double(new_bw_mask);
new_mask(new_mask==0) = NaN;
end