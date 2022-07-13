function result = get_img_max(hdr, name)
%
%function result = get_img_max(hdr, name)
%
% Luis hernandez
% last edit 8-24-98
%
% Loads a specific slice of data from an analyze format file 'name' 
% containing mutislice image data and finds the intensity of the brightest 
%pixel

   switch hdr.bits      
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

   
   [pFile,messg] = fopen(name, 'r');
   
   if pFile == -1
      disp(messg);   
      return;
   end
   
   d = (fread(pFile,hdr.xdim * hdr.ydim, fmt))' ;	
   
   fclose(pFile);
   
   result = max(d);
   clear d;
   %%%%
   
   return
   


