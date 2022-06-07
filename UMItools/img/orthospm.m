function orthospm(threshold,spm_file, anat_file )
% function orthospm(threshold,spm_file, anat_file )

global SPM_scale_factor wscale w2scale

if nargin<3
    [name path] = uigetfile('*.img','Select SPM *.img file');
    spm_file = strcat(path,name)

end
sz = size(spm_file);

imgname = strcat(spm_file(1,1:sz(2)-4) , '.img');
hdrname = strcat(spm_file(1,1:sz(2)-4) , '.hdr');


spm_hdr = read_hdr(hdrname);
spm_data = read_img2(spm_hdr,imgname);
spm_data(find(spm_data==NaN))=0;
spm_scale =  SPM_scale_factor

if nargin<3
    [name path] = uigetfile('*.img','Select anatomical *.img file');
    anat_file = strcat(path,name)
end
sz = size(anat_file);

imgname = strcat(anat_file(1,1:sz(2)-4) , '.img');
hdrname = strcat(anat_file(1,1:sz(2)-4) , '.hdr');

hdr = read_hdr(hdrname);
anat_data = read_img2(hdr,imgname);

% interpolate the paramter map to fit the anatomical one
% note the transpose....A weird quirk of meshgrid
[x,y,z] = meshgrid(1:spm_hdr.ydim , 1:spm_hdr.xdim, 1:spm_hdr.zdim);
[xi,yi, zi] = meshgrid(1:hdr.ydim , 1:hdr.xdim, 1:hdr.zdim);

% optional in case the origin is not set in the image headers...
%
% spm_hdr.origin = [0 0 0];
% hdr.origin = [0 0 0];
%

x = ( x - spm_hdr.origin(1) )* spm_hdr.xsize;
y = ( y - spm_hdr.origin(2) )* spm_hdr.ysize;
z = ( z - spm_hdr.origin(3) )* spm_hdr.zsize;

xi = ( xi - hdr.origin(1) )* hdr.xsize;
yi = ( yi - hdr.origin(2) )* hdr.ysize;
zi = ( zi - hdr.origin(3) )* hdr.zsize;


whos anat_data spm_data
hdr.origin'
spm_hdr.origin'

if hdr.zdim==1
    spm_data_sc = interp2(x,y, spm_data, xi,yi,'nearest');
else
    spm_data_sc = interp3(x,y,z, spm_data, xi,yi,zi,'nearest');
end
spm_data_sc(find(isnan(spm_data_sc)))=0;
spm_data_vals = spm_data_sc;


% threshold the map at the 50%
if nargin ==0
    threshold = mean(mean(mean(spm_data_sc))) * 3
end
spm_data_sc(find(spm_data_sc < threshold )) = 0;


% scale the maps to use the whole colormap
anat_data = anat_data * 256 / max(max(max(anat_data)));

if ~isempty(w2scale)
    spm_data_sc = (spm_data_sc)*256 / w2scale(2) ;
else
    spm_data_sc = (spm_data_sc) * 256 / max(max(max(spm_data_sc)));
end

out_data = anat_data;
% manual scale
if ~isempty(wscale)
    out_data = (anat_data - wscale(1))*256 / wscale(2);
end
out_data(find(spm_data_sc)) =  256 + spm_data_sc(find(spm_data_sc))+1;
% 
% %configure the colormap:
figure
set_func_colors

d=out_data;

% display the orthogonal sections
x = round(hdr.xdim/2);
y = round(hdr.ydim/2);
z = round (hdr.zdim/2);


stretch = hdr.zsize/hdr.xsize;

xs=round(x);ys=round(y);zs=round(z);

i=1;
while i >= -10
    % calculate voxel coords. for spm data point
    xs = round((x - hdr.origin(1))*hdr.xsize/spm_hdr.xsize + spm_hdr.origin(1));
    ys = round((y - hdr.origin(2))*hdr.ysize/spm_hdr.ysize + spm_hdr.origin(2));
    zs = round((z - hdr.origin(3))*hdr.zsize/spm_hdr.zsize + spm_hdr.origin(3));

    % display the data:
    [fig1, fig2, fig3] =  ov(hdr,d,x,y,z,0);

    str=sprintf('(x,y,z)=  (%3.2f %3.2f %3.2f) mm, \n (x,y,z)=(%3d %3d %3d)vx , val= %6.2f',...
        hdr.xsize*( x - hdr.origin(1)), ...
        hdr.ysize*( y - hdr.origin(2)), ...
        hdr.zsize*( z - hdr.origin(3)), ...
        xs,ys,zs, ...
        spm_data_vals(x,y,z)*spm_scale );
   fprintf('\n%s' , str);

    [i j] = ginput(1);
    i=round(i);j=round(j);
    % exiting the program when you click outside bounds
    if (i<=-1)
        %colordef white
        return
    end

    fig = round(gca);
    switch(fig)
        case round(fig1)
            x=j;
            y=i;
        case round(fig2)
            z=j;
            x=i;
        case round(fig3)
            y=i;
            z=j;
    end



end
%colordef white
return












