function mshowspm()
%

[name path] = uigetfile('*.img','Select SPM *.img file');
name = strcat(path,name)

sz = size(name);

imgname = strcat(name(1,1:sz(2)-4) , '.img');
hdrname = strcat(name(1,1:sz(2)-4) , '.hdr');

spm_hdr = read_hdr(hdrname);
spm_data = read_img2(spm_hdr,imgname);

    
[name path] = uigetfile('*.img','Select anatomical *.img file');
name = strcat(path,name)

sz = size(name);

imgname = strcat(name(1,1:sz(2)-4) , '.img');
hdrname = strcat(name(1,1:sz(2)-4) , '.hdr');

hdr = read_hdr(hdrname);
anat_data = read_img2(hdr,imgname);

[x,y] = meshgrid(1:spm_hdr.xdim , 1:spm_hdr.ydim);
[xi,yi] = meshgrid(1:hdr.xdim , 1:hdr.ydim);

yi = yi * spm_hdr.ydim/hdr.ydim;
xi = xi * spm_hdr.xdim/hdr.xdim;

spm_data2 = [];
for i=1:hdr.zdim
    %spm_data2(:,:,i) =linterp(spm_data(:,:,i), 4);
    spm_data2(:,:,i) = interp2(x,y,spm_data(:,:,i), xi,yi);
    fprintf('\rinterpolating ...%d', i);
end


spm_data = spm_data2;
spm_data(find(spm_data < 500 )) = 0;

anat_data = anat_data * 2^7 / max(max(max(anat_data)));
spm_data = spm_data * 2^7 / max(max(max(spm_data)));

out_data = anat_data;
out_data(find(spm_data)) =  2^7 + spm_data(find(spm_data))-1;
	


%xx=1;
%while(xx > 0)
%	orthov2(out_data, xx, yy, zz);

	mygray = [0:1/127:1]' * [1 1 1]; 
	myhot = [1:3:127]' * [1 0 0] ;
	tmp =   [1:3:127]' * [0 1 0] ;
	myhot = [myhot; tmp];
	tmp =   [1:3:127]' * [0 0 1];
	myhot =  [myhot;  tmp]/128;

	myhot(round(128/3): 128, 1) = 1;
	myhot(round(128*2/3):128,2) = 1;

 	mymap = [mygray; myhot];
	colormap(mymap)
	axis off
    
    
    
  % display the orthogonal sections
    colordef black 	
    colormap(gray)
    subplot 221, imagesc(squeeze(d(:,:,z))), axis image,axis xy 
    subplot 222, imagesc(squeeze(d(:,y,:))'), axis image,axis xy 
    subplot 223, imagesc(squeeze(d(x,:,:))'), axis image,axis xy
   

    i=0;
    while i >= -10
        [i j] = ginput(1);
        i=round(i);j=round(j);
        fig = floor(gca)
        switch(fig)
        case 99
            x=j;
            y=i;
        case 105
            z=j;
            x=i;
        case 111
            y=i;
            z=j;
        end
        [ x y z d(x,y,z)]
        
        colordef black
        
        subplot 221, imagesc(squeeze(d(:,:,z))), axis image,axis xy 
        subplot 222, imagesc(squeeze(d(:,y,:))'), axis image,axis xy 
        subplot 223, imagesc(squeeze(d(x,:,:))'), axis image,axis xy
        colordef white
    end 
    
return













