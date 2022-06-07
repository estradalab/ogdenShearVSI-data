function tdata = timeplot5(path, root, index )

% function result = timeplot5(path,root, index)
%
% Returns the intensity of the voxel the point determined by the index
% 
% 
oldpath=pwd;
cd (path);

files = dir(strcat(root,'*.img'));
if (size(files)==[0 1])
	tdata=0;
	fprintf('\n** ERROR:  %s ----- images not found **\n\n',root);
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
	%fprintf('reading ... %s\n', files(i).name);
	[fp mesg]= fopen(files(i).name);
	if fp == -1
		disp(mesg);
		return
	end
	
    
    
    fseek(fp,index*bytes,'bof');     
    tdata(i) = fread(fp,1,fmt);
    fclose(fp);
    
    
	
end
cd(oldpath);
%close

return





