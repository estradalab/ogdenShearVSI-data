%
% A little macro to get the whole volume of a ginx series.
%
% Robert Welsh - Aug-17-00
% University of Michigan
% Department of Radiology
%  
% Syntax 
%
%    [ginxVol, errmsg] = getginxvol;
% 

% Modification to make : allow the person to pass a file name 
% and hence remove any interaction with graphics.
% If this happens, the path name will have to be parsed.

function [ginxVol, errmsg] = getginxvol

[fname pathname] = uigetfile('*i*','Pick an anatomy file.');

% Did they happen to point to an already existing img file?

if ( findstr(fname,'.img'))
  imgHDR = spm_vol([pathname fname]);
  [imgVOL imgXYZ] = spm_read_vols(imgHDR);
  clear imgXYZ;
  ginxVol = imgVOL;
  for iz = 1:size(imgVOL,3)
    ginxVol(:,:,iz) = rot90(imgVOL(:,:,iz),2);
  end
  clear imgVOL;
  errmsg = 0;
else
  if (fname == 0 )
    ginxVol = 0;
    errmsg  = 1;
  else
    fileToOpen = [pathname fname];
    
    [aimg nslices dimx dimy errmsg] = getginx(fileToOpen);
    
    if ( errmsg ) 
      fprintf('Something wrong with root file.\n');
      errmsg = 2;
    else
      ginxVol = zeros(dimx,dimy,nslices);
      for iimg = 1:nslices
        fileToOpen = newFileName(pathname,fname,iimg);
        [aimg ndum dumx dumy errmsgrtn] = getginx(fileToOpen); 
        if (errmsgrtn)
	  fprintf('Possibly missing or corrupted slices %d\n',iimg); 
        else
	  ginxVol(:,:,iimg) = aimg;
        end
      end
    end
  end
end

%
% All done;
%

function [newFileName] = newFileName(pathname,fname,iimg);

for ipos = length(fname):-1:1
  if (fname(ipos:ipos) == 'i')
     newFileName = [pathname '/' fname(1:ipos) sprintf('%d',iimg)];
     break
  end
end

%
% All done
%
