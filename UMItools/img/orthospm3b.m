function orthospm3b(threshold, onsets, window ,spm_file, anat_file, tseries_path, func_root )
% function result = orthospm3b(  threshold , 
%                               onsets, 
%                               window 
%                               [,spm_file , 
%                                anat_file, 
%                                tseries_path, 
%                                func_root_name] )
%
% this function takes a file in analyze format and 
% displays orthogonal sections thru the planes interscting
% at x,y,z.  It then overlays a Statistical Parameter Map 
% thresholded at "threshold".
%
% additionally, this function allows you to extract a time series
% from a chose pixel, just by clicking on it.  The time series is saved
% in the file "tdata.dat" if you use the RIGHT MOUSE BUTTON.
% This file gets overwritten everytime you select a new pixel.
% 
% *** This version of the program calls event_avg.m 
% in order to average all the events into one.  No deconvolution!!
%
% onsets and window must be in scan units

global SPM_scale_factor 
close all
filter=1;
TR=1;

roi =0;  % this is the number of pixels around the pixel of interest
% that we'll include in the time series avg.
if nargin < 4
	[name path] = uigetfile('*.img','Select SPM *.img file');
	spm_file = strcat(path,name)
end

% figure out name for spm file
sz = size(spm_file);
imgname = strcat(spm_file(1,1:sz(2)-4) , '.img');
hdrname = strcat(spm_file(1,1:sz(2)-4) , '.hdr');
%label=sprintf('FUNCTIONAL:  ....%s',imgname(end-15:end));

spm_hdr = read_hdr(hdrname);
spm_data = read_img2(spm_hdr,imgname);
spm_data(find(spm_data==NaN))=0;
spm_scale =  SPM_scale_factor;

if nargin < 4
	[name path] = uigetfile('*.img','Select anatomical *.img file');
	anat_file= strcat(path,name)
end

% figure out name for anatomical file
sz = size(anat_file);

imgname = strcat(anat_file(1,1:sz(2)-4) , '.img');
hdrname = strcat(anat_file(1,1:sz(2)-4) , '.hdr');
%label=sprintf('%s\nANATOMY: ....%s', label, imgname(end-15:end));

hdr = read_hdr(hdrname);
anat_data = read_img2(hdr,imgname);

% interpolate the paramter map to fit the anatomical one
% note the transpose....
[x,y,z] = meshgrid(1:spm_hdr.ydim , 1:spm_hdr.xdim, 1:spm_hdr.zdim);
[xi,yi, zi] = meshgrid(1:hdr.ydim , 1:hdr.xdim, 1:hdr.zdim);

yi = yi * spm_hdr.ydim/hdr.ydim;
xi = xi * spm_hdr.xdim/hdr.xdim;
zi = zi * spm_hdr.zdim/hdr.zdim;

spm_data2 = interp3(x,y,z, spm_data, xi,yi,zi,'nearest');
spm_data2(find(isnan(spm_data2)))=0;
%spm_data = spm_data2;

%whos

% threshold the map at the 50%
if nargin ==0
	threshold = mean(mean(mean(spm_data2))) * 3
end
spm_data2(find(spm_data2 <= threshold )) = 0;

% scale the maps to use the whole colormap
anat_data = anat_data * 2^7 / max(max(max(anat_data)));
spm_data2 = spm_data2 * 2^7 / max(max(max(spm_data2)));

out_data = anat_data;
out_data(find(spm_data2)) =  2^7 + spm_data2(find(spm_data2));

%configure the colormap:
mygray = [0:1/127:1]' * [1 1 1]; 

myhot = [0:3:125]' * [1 0 0] ;
tmp =   [0:3:125]' * [0 1 0] ;
myhot = [myhot; tmp];
tmp =   [0:3:125]' * [0 0 1];
myhot =  [myhot;  tmp]/128;

myhot(round(128/3): 128, 1) = 1;
myhot(round(128*2/3):128,2) = 1;

mymap = [mygray; myhot];

colormap(mymap)
axis off


d=out_data;

fprintf('\nfirst display the orthogonal sections');
x=ceil(hdr.xdim/2);
y=ceil(hdr.ydim/2);
z=ceil(hdr.zdim/2);

%colordef black 	
stretch = hdr.zsize/hdr.xsize;

%colordef black
[fig1, fig2, fig3] =  ov(hdr,d,x,y,z,0);


if nargin < 4
	% determine which files contain the time series..:
	f='dummmyname';
	file=[];
	path = [];
	op=pwd;
	while f ~= 0
		[f p] = uigetfile('*.img','Select a file from the next run.  cancel when finished ');
		if f~=0
			file = [file ; f];
			path = [path ; p];
			cd(p);
		end
	end
	cd(op);
	numRuns = size(file,1);
end

[i j button] = ginput(1);
i=round(i);j=round(j);


%i=1;
while i >= -10
	
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
	
	
	xs = ceil(x*spm_hdr.xdim/hdr.xdim);
	ys = ceil(y*spm_hdr.ydim/hdr.ydim);
	zs = ceil(z*spm_hdr.zdim/hdr.zdim);
	
	str=sprintf('(x,y,z)=  (%3.2f %3.2f %3.2f) mm,  \n (%3d %3d %3d)vx , \n val= %d',...
		hdr.xsize*(xs-spm_hdr.origin(1)), ...
		hdr.ysize*(ys-spm_hdr.origin(2)), ...
		hdr.zsize*(zs-spm_hdr.origin(3)), ...
		xs,ys,zs, ...
		spm_data(xs,ys,zs)*spm_scale );
	
	%colordef black
	[ x y z] 
	[fig1, fig2, fig3] =  ov(hdr,d,x,y,z,0);
	subplot(221), title(str) %, xlabel(label)
	
	
	tmp=mean(mean(spm_data(xs-roi:xs+roi, ys-roi:ys+roi)));
	
	if nargin > 4
		numRuns=1;
		path=tseries_path;
		file=func_root;
	end
	
	tdata=[];
	
	g = make_gaussian(50,4,100);
	for n=1:numRuns
		tmp = timeplot2(path(n,:), file(n,:), xs,ys,zs);
%		tmp = smoothdata2(tmp, TR, 0.009, 0.3, 4); 
		tmp = detrend(tmp);
		tmp = conv(g,tmp);
		tmp = tmp(50:end-50);
	
		tdata = [tdata tmp];
	end
	
	
	[ev_avg ev_std] = event_avg(tdata,onsets,window,10);
	subplot 224, plot(ev_avg,'r');axis tight ; hold off; 
	%subplot 224, errorbar(ev_avg,ev_std);axis tight ; hold off; 
	set(gca, 'Xtick',[0:2:window])
	
	if button==3
		 tmp = [ev_avg ev_std];
                save avg.dat tmp -ASCII
                save tdata.dat tdata -ASCII
                save voxels.dat roi_xyz -ASCII

	end
	
	[i j button] = ginput(1);
	i=round(i);j=round(j);
	
end 
%colordef white
return
