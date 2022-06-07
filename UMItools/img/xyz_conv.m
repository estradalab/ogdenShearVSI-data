% Convert the XYZ coordinates that are given by SPM when
%   using results to pixel coordinates in fMRI space.
%   The following is the orientation as used by spm.
%       x: left -> right
%       y: post -> ant
%       z: inf  -> sup
%
% fMRI uses the (x,y)=(0,0) as the origin. This is the upper left
%   corner, or the rightmost and most anterior.
%
% When acquiring only part of the brain, we typicall specify the
%   origin as (x,y,z) = (32,32,0). This is in pixels and SPM
%   allows (0,0,0) as the origin, so adjust accordingly.
%   Thus, negative x is toward the left, negative y is toward
%   the posterior. Z is positive since we use 0 as the most inferior.
%
% Thus, a negative y should come out with a value less than 32
%   indicating more anterior.
%   Negative x is to the left, thus giving a value greater than 32.
%
function pix_mm = XYZ_conv(XYZ)
PIX_mm = 4.375;
SLI_mm = 5.00;
X_org = 32;
Y_org = 32;
x_mm = XYZ(1,:)';
y_mm = XYZ(2,:)';
z_mm = XYZ(3,:)';
x_pix = X_org - x_mm/PIX_mm;
y_pix = Y_org - y_mm/PIX_mm;
z_pix = z_mm/SLI_mm;
pix_mm = [x_pix y_pix z_pix];
