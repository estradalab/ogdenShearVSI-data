function result=spm_info()

% This program allows the user to select an img file containing
% an SPM and performs different statistics on the number and location 
% of acivated voxels


	% Let user select filename ...
	[file path] = uigetfile('*.img','Select Analyze file');
	imgname = strcat(path,file);
	
	sz = size(imgname);

	hdrname = strcat(imgname(1,1:sz(2)-4) , '.hdr');
	hdr = read_hdr(hdrname);
	
	data=extractspm(imgname);

	plane=[	hdr.xdim/2  	hdr.xdim/2  	0;
		hdr.xdim/2  	0          	0;
		hdr.xdim/2	hdr.xdim/2	hdr.zdim]
	%tmswin
	rl = rightleft(data, plane)

return
