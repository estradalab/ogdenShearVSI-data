function result = extractspm2()
%
% function result = extractspm2()
%
% Luis Hernandez
% lastedit: 3-24-98
% returns the xyz coordinates (in pixel units) and the
% intensities of all the non-zero elementsin an
% analyze format image set
%
% It also computes the center of mass of the spm maps
% and it counts the number of voxels on each side of the 
% center plane.

[name path] = uigetfile('*.img','Select SPM *.img file');

name = strcat(path,name)

sz = size(name);

imgname = strcat(name(1,1:sz(2)-4) , '.img');
hdrname = strcat(name(1,1:sz(2)-4) , '.hdr');

hdr = read_hdr(hdrname);
spm_volume = read_img2(hdr,imgname);
 
%read slice by slice to reduce size of data buffers
for i=1: hdr.zdim

	clear tmpx tmpy tmpx val
   
	[tmpx tmpy val]=find( spm_volume(:,:,i)  );
   tmpz=tmpx;
   tmpz(:) = i;
   
	data = [data; tmpx tmpy tmpz val];

end



whos data
disp('ending extractspm2');
result = data;

h = findobj('Tag','SPMFileName');
set(h,'String',imgname);

return


