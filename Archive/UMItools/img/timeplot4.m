function tdata = timeplot4(path, file, xyz)

% function result = timeplot4(path,file, [ x,y,z] )
%
% Returns the intensity of the voxels in the points determined
% x y z are the coordinates of a set of voxels to be averaged 
% the data comes back as a row vector. eg -
%   [ x1 y1 z1 ;
%     x2 y2 z2 ;
%     x3 y3 z3 ]
%
% (c) 2005 Luis Hernandez-Garcia 
% University of Michigan
% report bugs to:  hernan@umich.edu
%
% 
oldpath=pwd;
cd (path);
x=xyz(:,1);
y=xyz(:,2);
z=xyz(:,3);

Nvoxels = max(size(x))

sz = size(file);
root = file(1,1:sz(2)-8);

files = dir(strcat(root,'*.img'));
if (size(files)==[0 1])
	tdata=0;
	fprintf('\n** ERROR:  %s ----- images not found **\n\n',file);
	return;
end

hfiles = dir(strcat(root,'*.hdr'));
sz = size(files);
hfiles(1).name;
hdr = read_hdr(hfiles(1).name);


switch hdr.datatype     
	case 0
		fmt = 'int8';
		bytes = 1;
		
	case 2
		fmt = 'uint8';
		bytes = 1;
	case 4
		fmt = 'short';
		bytes = 2;
	case 8
		fmt = 'int';
		bytes = 2;
	case 16
		fmt = 'float';
		bytes = 4;
	case 32
		fmt = 'float';
		xdim = hdr.xdim * 2;
		ydim = hdr.ydim * 2;
		bytes = 8;
		
	otherwise
		errormesg(sprintf('Data Type %d Unsupported. Aborting',hdr.bits));
		return
end

% extract data from files
tdata = zeros(1,sz(1));

for i=1:sz(1)
	fprintf('\rreading ... %s\n', files(i).name);
	[fp mesg]= fopen(files(i).name);
	if fp == -1
		disp(mesg);
		return
	end
	
	% calculate the position in the file corresponding to 
	% the voxel coordinates
	position = sub2ind([hdr.xdim hdr.ydim hdr.zdim],x,y,z);
		
	% average the ROI positions
	tmp = fread(fp,hdr.xdim*hdr.ydim*hdr.zdim,fmt);
	tmp = tmp(position);

	tdata(i) = mean(tmp);
	fclose(fp);
	
	
end
cd(oldpath);
%close

Nvoxels
%plot(tdata/max(tdata))  
%save timeseries tdata  

return





