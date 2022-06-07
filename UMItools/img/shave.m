function shave(rootname, slice)
% function shave(rootname, slice)
%
% removes slices from specified image file
% assumes analyze format
%

h=read_hdr(strcat(rootname,'.hdr'));

switch h.bits      
	case 8
         fmt = 'uint8';
      case 16
         fmt = 'uint16';
      case 32
         fmt = 'uint32';
      otherwise
         errormesg(sprintf('Data Type %d Unsupported. Aborting',hdr.bits));
         return
end


pFile = fopen(strcat(rootname,'.img'), 'r' );
data_in = fread(pFile, h.xdim*h.ydim*h.zdim, fmt);
fclose(pFile);
sz = size(data_in)

data_out = data_in(1:(slice-1)*h.xdim*h.ydim);
data_out = [data_out;  data_in( slice*h.xdim*h.ydim : sz(1)) ];

size(data_out)

pFile = fopen('out.img', 'wb');
fwrite(pFile,data_out,fmt);
fclose(pFile);

h.zdim = h.zdim-1;
write_hdr('out.hdr', h);

return



