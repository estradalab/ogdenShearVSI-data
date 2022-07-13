function zero_slice(rootname, slice)
% function zero_slice(rootname, slice)
%
% sets to zero all intensity velues in the 
% slice from specified image file
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

data_in((slice-1)*h.xdim*h.ydim: slice*h.xdim*h.ydim) = 0;

pFile = fopen(strcat(rootname,'.img'), 'wb');
fwrite(pFile,data_in,fmt);
fclose(pFile);


return



