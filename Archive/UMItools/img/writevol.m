function  writevol(slices)

% function writevol(slices)
% Luis hernandez
% last edit 12-29-97
%
% Saves the data in ANALYZE format.  (presently tha header structure 
% is unknown, but the *img file can be created from the structure data
% as a3-D array of shorts


   %Select file name and open for writing
   [name, path]= uiputfile('*.img','Save Analyze format Map As ...');
   pFile = fopen(name , 'w');
   
   % write data
   for sl=1:slices(1).n_slices,
      fwrite(pFile, slices(sl).data', 'short');
	end
   
   fclose(pFile);

   % This is where we will write the header from the analyze format file.
   % (YET TO BE DONE) 


return