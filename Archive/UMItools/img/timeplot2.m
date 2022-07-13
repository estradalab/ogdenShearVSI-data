function tdata = timeplot2(path, file,x,y,z)
% function result = timeplot2(path,file,x,y,z)
%
% Returns the intensity of the voxels in the cube determined
% by the two points (x(1), y(1)), z(1)) and (y(1), y(2)) z(2)) in a time
% series of fmri???.img.  
% the data comes back as a row vector
%
    oldpath=pwd;
    cd (path);
    
    sz = size(file);
    root = file(1,1:sz(2)-8);
    
    files = dir(strcat(root,'*.img'));
    if (size(files)==[0 1])
        tdata=0;
        fprintf('%s-----images not found',file);
        return;
    end
    
    hfiles = dir(strcat(root,'*.hdr'));
    sz = size(files);
	hfiles(1).name;
	hdr = read_hdr(hfiles(1).name);
	
	

	% determine which positions are included into the ROI
	num = 1;

	startx = min(x);
	starty = min(y);
	startz = min(z);
    endx = max(x);
	endy = max(y);
    endz = max(z);
	
    
    for i=startx:endx
          for j=starty:endy
                for k=startz:endz
                      position(num) = (k-1)*(hdr.xdim*hdr.ydim) + (j-1)*hdr.xdim + i-1;
                      num = num +1;
                end
          end
    end

	num = num -1;

	% extract data from files
	tdata = zeros(1,sz(1));

	for i=1:sz(1)

		[fp mesg]= fopen(files(i).name);
      	%	disp(i)
		%disp  ( files(i).name)

		if fp == -1
			disp(mesg);
			return
		end
	   
      
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


		% average the ROI positions
		for n=1:num      
		 	fseek(fp,(position(n))*bytes,'bof');     
		 	data = fread(fp,1,fmt);
		  	tdata(i) = tdata(i) + data;
		end

		tdata(i) = tdata(i) / num;

	  	fclose(fp);
	end

    cd(oldpath);
	%close
	
	
	%plot(tdata/max(tdata))  
	%save timeseries tdata  
	
	return
			




