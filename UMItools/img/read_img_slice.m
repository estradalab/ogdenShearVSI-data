function slice = read_img_slice(hdr, name, sl)
%
%	slice = read_img(hdr, name,sl)
%
% Luis hernandez
% last edit 6-30-98
%
% Loads a specific slice of data from an analyze format file 'name' 
% containing mutislice image data
%
% The function returns an array of slices
% The slices data structure is defined as 
%   slice = struct(...%
%     'xdim',hdr.xdim,...      % number of pixels in the x direction
%     'ydim',hdr.ydim,...      % number of pixels in the y direction
%     'slicenum', 1,...        % current slice number within the volume
%     'n_slices',hdr.zdim,...  % Number of slices in the volume
%     'data',zeros(1, 1));  % data proper


   
   slice = struct(...
     'xdim',hdr.xdim,...      % number of pixels in the x direction
     'ydim',hdr.ydim,...      % number of pixels in the y direction
     'slicenum', 1,...        % current slice number within the volume
     'n_slices',hdr.zdim,...  % Number of slices in the volume
     'data',zeros(1, 1));  % data proper

%% 

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

	    
   % Create a slice structure combining description and data
   
   slice.xdim = hdr.xdim;
   slice.ydim = hdr.ydim;
   slice.n_slices = hdr.zdim;	      
   
   
   [pFile,messg] = fopen(name, 'r');
   if pFile == -1
      disp(messg);   
      return;
   end
   
   skip = (sl-1)*(hdr.xdim*hdr.ydim)*(hdr.bits/8);
   fseek(pFile,skip,'bof');  
   slice.data = (fread(pFile,[hdr.xdim, hdr.ydim], fmt))' ;	
   
   fclose(pFile);
   
   %%%%
   
   return
   


