function slice = read_img(hdr, name)

%	slice = read_img(hdr, name)
%
% Luis hernandez
% last edit 1-7-98
%
% Loads the data from an analyze format file 'name' containing mutislice image data
%
% The function returns an array of slices
% The slices data structure is defined as 
%   slice = struct(...%
%     'xdim',hdr.xdim,...      % number of pixels in the x direction
%     'ydim',hdr.ydim,...      % number of pixels in the y direction
%     'slicenum', 1,...        % current slice number within the volume
%     'n_slices',hdr.zdim,...  % Number of slices in the volume
%     'data',zeros(1, 1));  % data proper


   [pFile,messg] = fopen(name, 'r');
   if pFile == -1
      msgbox(messg);   
      return;
   end
   
   slice = struct(...%
     'xdim',hdr.xdim,...      % number of pixels in the x direction
     'ydim',hdr.ydim,...      % number of pixels in the y direction
     'slicenum', 1,...        % current slice number within the volume
     'n_slices',hdr.zdim,...  % Number of slices in the volume
     'data',zeros(1, 1));  % data proper

%% vvvvvv

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
	for sl=1:hdr.zdim,
	      slice(sl).xdim = hdr.xdim;
	      slice(sl).ydim = hdr.ydim;
	      slice(sl).n_slices = hdr.zdim;	      
		slice(sl).data = (fread(pFile,[hdr.xdim, hdr.ydim], fmt))' ;	
	end
	fclose(pFile);

%%%%^^^^^   
   
   return



