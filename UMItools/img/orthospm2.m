function orthospm2(roi, threshold, SMname,  ANAname, TSpath, TSname)
% function result = orthospm2( roi, 
%                           threshold, 
%                         [ SMname,  
%                           ANAname, 
%                           TSpath, 
%                           TSname ] )
%
% this function takes a file in analyze format and 
% displays orthogonal sections thru the planes interscting
% at x,y,z.  It then overlays a Statistical Parameter Map 
% thresholded at "threshold".
%
% additionally, this functiona allows you to extract a time series
% from a chose pixel, just by clicking on it.  The time series is saved
% in the file "tdata.dat" if you use the RIGHT MOUSE BUTTON.
% This file gets overwritten everytime you select a new pixel.
%
% assumption: time series and SPM map have the sampe dimensions
%
 
global SPM_scale_factor wscale myfig
figure

    %roi =0;  % this is the number of pixels around the pixel of interest
            % that we'll include in the time series avg.
    if (nargin<3)
        [name path] = uigetfile('*.img','Select SPM *.img file');
        name = strcat(path,name)
    else
        name = SMname;
    end
    sz = size(name);
    
    imgname = strcat(name(1,1:sz(2)-4) , '.img');
    hdrname = strcat(name(1,1:sz(2)-4) , '.hdr');
    
    spm_hdr = read_hdr(hdrname);
    spm_data = read_img2(spm_hdr,imgname);
    spm_data(find(spm_data==NaN))=0;
    spm_scale =  SPM_scale_factor;

    if nargin<3
        [name path] = uigetfile('*.img','Select anatomical *.img file');
        name = strcat(path,name)
    else
        name = ANAname;
    end
       
    sz = size(name);
    
    imgname = strcat(name(1,1:sz(2)-4) , '.img');
    hdrname = strcat(name(1,1:sz(2)-4) , '.hdr');
    
    hdr = read_hdr(hdrname);
    anat_data = read_img2(hdr,imgname);
    
    % interpolate the paramter map to fit the anatomical one
	% note the transpose....
    [x,y,z] = meshgrid(1:spm_hdr.ydim , 1:spm_hdr.xdim, 1:spm_hdr.zdim);
    [xi,yi, zi] = meshgrid(1:hdr.ydim , 1:hdr.xdim, 1:hdr.zdim);
    
    x = ( x - spm_hdr.origin(1) )* spm_hdr.xsize;
    y = ( y - spm_hdr.origin(2) )* spm_hdr.ysize;
    z = ( z - spm_hdr.origin(3) )* spm_hdr.zsize;
    
    xi = ( xi - hdr.origin(1) )* hdr.xsize;
    yi = ( yi - hdr.origin(2) )* hdr.ysize;
    zi = ( zi - hdr.origin(3) )* hdr.zsize;
% 
%     yi = yi * spm_hdr.ydim/hdr.ydim;
%     xi = xi * spm_hdr.xdim/hdr.xdim;
%     zi = zi * spm_hdr.zdim/hdr.zdim;
%    
    if hdr.zdim==1
	spm_data2 = interp2(x,y, spm_data, xi,yi,'nearest');
    else
	spm_data2 = interp3(x,y,z, spm_data, xi,yi,zi,'nearest');
    end 
    spm_data2(find(isnan(spm_data2)))=0;
    %spm_data = spm_data2;
    
   % whos
    
    % threshold the map at the 50%
    if nargin ==0
        threshold = mean(mean(mean(spm_data2))) * 3
    end
    spm_data2(find(spm_data2 < threshold )) = 0;
    
    % scale the maps to use the whole colormap
    anat_data = anat_data * 256 / max(max(max(anat_data)));
    spm_data2 = spm_data2 * 256 / max(max(max(spm_data2)));
    
    out_data = anat_data;
    % manual scale
    if ~isempty(wscale)
         out_data = (anat_data-wscale(1))*256 / wscale(2);
    end
    out_data(find(spm_data2)) =  256 + spm_data2(find(spm_data2))+1;
    
    %configure the colormap:
    set_func_colors

    
        
    d=out_data;

  % display the orthogonal sections
    x=hdr.xdim/2;
    y=hdr.ydim/2;
    z=hdr.zdim/2;

        xs = ceil(x*spm_hdr.xdim/hdr.xdim);
        ys = ceil(y*spm_hdr.ydim/hdr.ydim);
        zs = ceil(z*spm_hdr.zdim/hdr.zdim);
        
    xs=ceil(xs); ys=ceil(ys);z=ceil(zs);
	
    colordef black 	
    stretch = hdr.zsize/hdr.xsize;
    
     %[fig1, fig2, fig3] =  ov(hdr,d,x,y,z,roi);
    % determine which files contain the time series..:
	if nargin<3
        [file path] = uigetfile('*.img','Select one of the Analyze files in the series ');
    else
        file = TSname;
        path = TSpath;
    end
    %%%%%%
    fprintf('\nLoading the time series ...%s\n', file(1:end-8))
    wholedata = read_img_series(file(1:end-8));
    whos wholedata
    fprintf('\n...done\n')    
    if hdr.zdim==1
           tdata = xtract(spm_hdr, wholedata, [xs-roi x+roi],[ys-roi ys+roi],[1 1]);
    else
       tdata = xtract(spm_hdr, wholedata, [xs-roi xs+roi],[ys-roi ys+roi],[zs-roi zs+roi]);
    end
    %%%%%%%%%%%%%%%%%%%%%%%
    

    i=1;
    while i >= -10
        xs = ceil(x*spm_hdr.xdim/hdr.xdim);
        ys = ceil(y*spm_hdr.ydim/hdr.ydim);
        zs = ceil(z*spm_hdr.zdim/hdr.zdim);
        
        fprintf('\n picked x y z value :')
        disp([ xs ys zs spm_data(xs,ys,zs)*spm_scale ])
        
   
		[fig1, fig2, fig3] =  ov(hdr,d,x,y,z,roi);
		
		str=sprintf('(x,y,z)=  (%3.2f %3.2f %3.2f) mm,  \n (%3d %3d %3d)vx , \n val= %d',...
        	hdr.xsize*(xs-hdr.origin(1)), ...
        	hdr.ysize*(ys-hdr.origin(2)), ...
        	hdr.zsize*(zs-hdr.origin(3)), ...
			xs,ys,zs, ...
 			spm_data(xs,ys,zs)*spm_scale );
        
		fprintf('\n%s' , str);     
		        
        [i j button] = ginput(1);
		i=round(i);j=round(j);
		% exiting the program when you click outside bounds
		if (i<=-1)   
			colordef white
			return
		end
		
        fig = floor(gca);
        switch(fig)
        case floor(fig1)
            x=j;
            y=i;
        case floor(fig2)
            z=j;
            x=i;
        case floor(fig3)
            y=i;
            z=j;
        end
    	x=ceil(x); y=ceil(y);z=ceil(z);
	       
        tdata = xtract(spm_hdr, wholedata, [xs-roi xs+roi],[ys-roi ys+roi],[zs-roi zs+roi]);
  		% tdata = timeplot2(path, file,[xs-roi xs+roi],[ys-roi y+roi],[zs-roi z+roi]);
        subplot 224, plot(tdata);
        
        if button==3
            save tdata.dat tdata -ASCII
        end

    end 
    colordef white
    return













