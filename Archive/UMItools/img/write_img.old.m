function  write_img(slice, name)

% Luis hernandez
% last edit 1-7-98
%
% Writes the data to an analyze format file 'name' containing mutislice image data
%
% The function returns an array of slices
% The slices data structure is defined as 
%   slice = struct(...%
%     'xdim',hdr.xdim,...      % number of pixels in the x direction
%     'ydim',hdr.ydim,...      % number of pixels in the y direction
%     'slicenum', 1,...        % current slice number within the volume
%     'n_slices',hdr.zdim,...  % Number of slices in the volume
%     'data',zeros(1, 1));  % data proper



   [pFile,messg] = fopen(name, 'wb');
   if pFile == -1
      errormesg(messg);   
      return;
   end
   
  
   % Create a slice structure combining description and data
  	for sl=1:slice(1).n_slices,
        for i=1:slice(sl).xdim,
           for j=1:slice(sl).ydim,
              index = (sl-1)*slice(sl).xdim*slice(sl).ydim + (i-1)*slice(sl).xdim  + j;
              dummy(index) = slice(sl).data(i,j);
			end
		end
	end
   
   dummy = fread(pFile, 'uint16');
   fclose(pFile);
   
   return