function result = extractspm()

% function result = extractspm()
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
slices = read_img(hdr, imgname);

data = [0 0 0 0];

for z=1: slices(1).n_slices
	for y=1:hdr.ydim
		for x=1:hdr.xdim
			if slices(z).data(x,y) > 0
				data = [data; x y z slices(z).data(x,y)];				
			end
		end
	end
end


sz = size(data);
data = data(2:sz(1),:);
%  swap x and y for display purposes.
data = [data(:,2)  data(:,1) data(:,3)  data(:,4)];

sz = size(data);

result = data;

return


