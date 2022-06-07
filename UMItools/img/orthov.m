function result = orthov( roi, root, x, y, z)
% function result = orthov(roi_size, rootname, x, y, z)
%
% this function thakes a file in analyze format and
% displays orthogonal sections thru the planes interscting
% at x,y,z
global wscale
global ACTIVE mainfig
ACTIVE=1;

mainfig=figure;
colormap(gray)

if nargin==0
    roi=0;
end

if nargin < 2
    [name path] = uigetfile('*.img','Select  *.img file');
    name = strcat(path,name)

    sz = size(name);

    imgname = strcat(name(1,1:end-4) , '.img');
    hdrname = strcat(name(1,1:end-4) , '.hdr');
    h = read_hdr(hdrname);
    d = read_img2(h, imgname);
else
    %figure out the whole name of the header and image files
    % and read them.  Note the use of strcat
    if root(end-3:end)=='.hdr'
        root=root(1:end-4);
    elseif root(end-3:end)=='.img'
        root=root(1:end-4);
    end
    h = read_hdr(strcat(root,'.hdr'));
    d = read_img2(h, (strcat(root,'.img') ));
end


if nargin < 4
    x=h.xdim/2;
    y=h.ydim/2;
    z=h.zdim/2;
end

fprintf('\nheader info: %d x %d x %d', h.xdim, h.ydim, h.zdim);
fprintf('\nvox size: %d x %d x %d', h.xsize, h.ysize, h.zsize);

x=ceil(x); y=ceil(y); z=ceil(z);

i=0;
% configure colormap
my_map=(0:255)';
my_map=[my_map my_map my_map]/256;
colormap(my_map);

% scale image to fit colormap
range= max(max(max(d))) - min(min(min(d)));
dd = (d-min(min(min(d))))*256/range;
if ~isempty(wscale)
    dd = (d-wscale(1))*256 / wscale(2);
end



while (ACTIVE==1)
    [fig1, fig2, fig3] = ov(h,dd,x,y,z, roi);
    subplot(224), hist(reshape(d,h.xdim*h.ydim*h.zdim,1),100)
    if roi>0
        tmp=mean(mean(d(x-roi:x+roi, y-roi:y+roi,z-roi:z+roi)));
    else
        tmp=d(x,y,z);
    end

    str = sprintf('\n(x,y,z)=  (%d %d %d) , val= %6.2f  \n', x, y, z, tmp);
   fprintf('\n%s' , str);

    [i j] = ginput(1);
    i=round(i);j=round(j);
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
    if i <=-1 , ACTIVE=0, end;




end


return
