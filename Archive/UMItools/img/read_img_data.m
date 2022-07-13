function data = read_img_data(hdr, name)

%	data = read_img_data(hdr, name)
%
% Luis hernandez
% last edit 4-5-98
%
% Loads the data from an analyze format file 'name' containing mutislice image data
% this function returns a one dimensional array of data (row vector).
%
% (c) 2005 Luis Hernandez-Garcia 
% University of Michigan
% report bugs to:  hernan@umich.edu
%

global SPM_scale_factor endian


if (abs(SPM_scale_factor < 0.000000000001)) 
	SPM_scale_factor=1;
end
%SPM_scale_factor


   [pFile,messg] = fopen(name, 'r', endian);
   if pFile == -1
      fprintf('%s - could not open %s', messg, name);
      return;
   end
   
   xdim = hdr.xdim;
   ydim = hdr.ydim;
   zdim = hdr.zdim;
   
   switch hdr.datatype     
       case 0
           fmt = 'uint8';
       case 4
           fmt = 'short';
       case 8
           fmt = 'int';
       case 16
           fmt = 'float';
       case 32
           fmt = 'float';
           xdim = hdr.xdim * 2;
           ydim = hdr.ydim * 2;
       case 64
           fmt = 'int64';
       otherwise
           errormesg(sprintf('Data Type %d Unsupported. Aborting',hdr.bits));
           return

   end

	   
   
   % Read in data.
 
	data = (fread(pFile,[xdim*ydim*zdim], fmt))'; 	

	fclose(pFile);
	data = data * SPM_scale_factor;
   
   return





