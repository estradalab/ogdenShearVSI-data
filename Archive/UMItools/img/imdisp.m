function imdisp(data, scale);     %, slice_number)
%
% function imdisp(data, scale)
%
% Luis Hernandez
% 7-25-2000
%
% Display the current 2D image with the appropriate gray scale
% displays a 2D array. 
%

cla
[x y] = size(data);
pixels = max(x,y);
buffer = zeros(pixels, pixels);

buffer(1:x, 1:y) = data;

image(buffer'*scale);
  
return 


