function mimdisp(hdr,scale, vol_data);     %, slice_number)
%
% function imdisp(hdr, scale, vol_data)
%
% Multi-slice display
%
% Luis Hernandez
% 7-25-2000
%
% Display the current image on the axis scaled with scale
% 
global MAPMAX

matrix_size = ceil(sqrt(hdr.zdim))

xdim = matrix_size * hdr.xdim;
ydim = matrix_size * hdr.ydim;

pixels = max(xdim,ydim);

b=zeros(pixels,pixels);

row = 1;
col = 1;
sl = 1;

while sl <= hdr.zdim

   [sl row col];
   
   xlo = hdr.xdim*(row-1);
   xhi = hdr.xdim*row;
   ylo = hdr.ydim*(col-1);
   yhi = hdr.ydim*col;
   
   b( xlo + 1 : xhi,   ylo + 1 : yhi) = vol_data(:,:,sl);	
   
   col=col+1;
   sl = sl+1;
   
   if col > matrix_size
      col=1;
      row=row+1;
   end
   
end
mx = max(b);
scale = MAPMAX/max(mx);
image(b'*scale);

  
return 






