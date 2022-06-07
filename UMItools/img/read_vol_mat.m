function vol=read_vol_mat(anlzfile, dims, fmt)
% Usage ... vol=read_vol_mat(anlzfile,dimensions)
%
% anlzfile - is the analyze .img file
% dims - image dimansions [#slcs xdim ydim]
% fmt - file format (float, int16, ...etc)

anlzfid=fopen(anlzfile,'r');
if (anlzfid<3), error('Could not open file!'); end;

for m=1:dims(1),
  [im,cnt]=fread(anlzfid,dims(2:3),fmt);
  if (cnt~=prod(dims(2:3))), 
    error(['Could not read image ',int2str(m),'! (',int2str(cnt),')']); 
  end;
  vol(:,:,m)=im;
end;

fclose(anlzfid)

return
