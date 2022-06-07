function  write_img_data(name, data, hdr)

% Luis hernandez
% last edit 4-6-98
%
% function  write_img_data(name, data, hdr)
%
% Writes the data to an analyze format file 'name' containing mutislice image data 
%
%
% (c) 2005 Luis Hernandez-Garcia 
% University of Michigan
% report bugs to:  hernan@umich.edu
%



   [pFile,messg] = fopen(name, 'wb');
   if pFile == -1
      errormesg(messg);   
      return;
   end
   
      
   switch hdr.datatype     
   case 2
      fmt = 'char';
   case 4
      fmt = 'short';
   case 8
      fmt = 'int';
   case 16
      fmt = 'float';
   case 32
      fmt = 'float';
           
   otherwise
      errormesg(sprintf('Data Type %d Unsupported. Aborting',hdr.datatype));
      return
      
   end

   
   fwrite(pFile, data, fmt);
   fclose(pFile);
   
   
 return

