function result = avg_fmri(inroot, outroot, first, last)
%
% function result = avg_fmri(inroot, outroot, first, last)
%
% Luis Hernandez
% last edit 4-6-98
%
% This function averages several time series into one.  
%
% inroot is the root name for the input image files.
% outroot is the root name for the outpur image files
% first is the first frame file of the time series
% last is the last frame of the first time series
% 
% This function will average all the following frames into the  
% first time series.  i.e - to be used when one time series is repeated
% several times.


hdr_name = strcat(inroot,addzeros(first));
hdr_name = strcat(hdr_name,'.hdr');
hdr = read_hdr(hdr_name);

files = dir(strcat(inroot,'*.img'));
sz = size(files);

timepoints = last - first + 1;
reps = (sz(1) - first+1)/timepoints;

for j = first:last
   
   data = zeros(1, hdr.xdim*hdr.ydim*hdr.zdim);
   for i=0:reps-1
  
      imgname = strcat(inroot,addzeros(i*timepoints + j )  );
      imgname = strcat(imgname,'.img');
      disp(strcat('Adding ... %s',imgname));
      data = data + read_img_data(hdr,imgname);
      
   end
   
   data = data / reps;
   outname = strcat( outroot, addzeros(j) );
   disp(outname);

   write_hdr( strcat(outname,'.hdr'), hdr );
   write_img_data(strcat(outname,'.img'), data );
   
end


return

%%%%%%%%%%%%%%%%%%

function string = addzeros(number)
%
% function string = addzeros(number)
%
% This function pads the frame number string with zeros to match the 
% format required for SPM, etc...
%

if (number < 100)
   zstring = '0';
   if (number < 10)
      zstring = '00';
   end
   
   string = strcat(zstring,num2str(number) );
else
   string = num2str(number);
end


return







