function result = orthoft(roi, root, x, y, z, iscomplex)
% function result = orthoft( roi, root, x, y, z [,iscomplex])
% ...or    result=orthoft()
%
% this function thakes a file in analyze format and 
% displays orthogonal sections thru the planes interscting
% at x,y,z
%
% additionally, this functiona allows you to extract a time series
% from a chose pixel, just by clicking on it.  The time series is saved
% in the file "tdata.dat" if you use the RIGHT MOUSE BUTTON.
% This file gets overwritten everytime you select a new pixel.
%

%roi =0;  % this is the number of pixels around the pixel of interest
% that we'll include in the time series avg.
global ACTIVE mainfig ftfig
global wscale

ACTIVE = 1;

if nargin==0       
	roi=0;
end

if nargin < 2
	[name path] = uigetfile('*.img','Select Anatomical *.img file');
	name = strcat(path,name)
	
	sz = size(name);
	
	imgname = strcat(name(1,1:sz(2)-4) , '.img');
	hdrname = strcat(name(1,1:sz(2)-4) , '.hdr');
	h = read_hdr(hdrname);
	d = read_img2(h, imgname);
	
	[file path] = uigetfile('*.img','Select one of the Analyze files in the series ');
else
	
	%figure out the whole name of the header and image files
	% and read them.  Note the use of strcat
	
	h = read_hdr(strcat(root,'.hdr'));
	d = read_img2(h, (strcat(root,'.img') ));   
	
	file = sprintf('%s.img',root);
	path = pwd;
	
end


if nargin < 3
	x=ceil(h.xdim/2);
	y=ceil(h.ydim/2);
	z=ceil(h.zdim/2);
end

if nargin <6
	iscomplex=0
end

%%%%%%
fprintf('\nLoading the time series ...%s\n', file(1:end-8))
wholedata = read_img_series(file(1:end-8));
whos wholedata
fprintf('\n...done\n')    
tdata = xtract(h, wholedata, [x-roi x+roi],[y-roi y+roi],[z-roi z+roi]);

if (iscomplex)
    % read phase images and form complex time series.

    wholePdata = read_img_series( sprintf('p_%s', file(1:end-8)) );   
    ptmp = xtract(h, wholePdata, [x-roi x+roi],[y-roi y+roi],[z-roi z+roi]);
    whos ptmp mtmp
    
    tdata = tdata .* exp(-i*ptmp/1000);
end
%%%%%%%%%%%%%%%%%%%%%%%


stretch = h.zsize/h.xsize;

% display the orthogonal sections
mygray = [0:1/255:1]' * [1 1 1]; 
figure;
mainfig=gcf;
colordef black 	
colormap(mygray)
%%%%%%
% d2 = d* 2^8 / max(max(max(d)));
% if ~isempty(wscale)
%          d2 = (d-wscale(1))*256 / wscale(2);
% end

% scale image to fit colormap
range= max(max(max(d))) - min(min(min(d)));
d2 = (d-min(min(min(d))))*256/range;
if ~isempty(wscale)
    d2 = (d-wscale(1))*256 / wscale(2);
end

if ~isempty(wscale)
    d2 = (d-wscale(1))*256 / wscale(2);
end

nx=0;
[fig1, fig2, fig3] =  ov(h,d2,x,y,z,0);

figure;
ftfig=gcf;
set(gcf,'Position',[1 1 560,420]);
figure(mainfig);
k=0; j=0;
while (ACTIVE==1)
    [k j button] = ginput(1);
    k=round(k);j=round(j);
    fig = floor(gca);

    switch(fig)
        case floor(fig1)
            x=j;
            y=k;
        case floor(fig2)
            z=j;
            x=k;
        case floor(fig3)
            y=k;
            z=j;
    end 
   

    if k < 0 , ACTIVE=0, end;
    
	[fig1, fig2, fig3] =  ov(h,d2,x,y,z,0);
   
    title(fig3,sprintf('coords: %d %d %d',x,y,z));
    colordef white
    
    tdata = xtract(h, wholedata, [x-roi x+roi],[y-roi y+roi],[z-roi z+roi]);    
    if (iscomplex)
        % read phase images and form complex time series.        
        ptmp = xtract(h, wholePdata, [x-roi x+roi],[y-roi y+roi],[z-roi z+roi]);    
        tdata = tdata .* exp(-i*ptmp/1000);
    end
    
    subplot 224, plot(tdata),  title('Time series Magnitude');
	

	%show the spectrum ...
	figure(ftfig),hold off,  
    subplot(211), plot(linspace(-1,1,length(tdata)), abs(fftshift(fft(tdata - mean(tdata)))),'r')
    title('FFT magnitude (mean removed)')
    subplot(212), plot(linspace(-1,1,length(tdata)), angle(fftshift(fft(tdata - mean(tdata)))),'r')
    title('FFT Phase(mean removed)')
	figure(mainfig) 
	
	if button==3
		save tdata.dat tdata -ASCII
	end
	
end 


return

