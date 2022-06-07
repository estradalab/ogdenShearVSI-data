function results = animate_series(root,x,y,z, scale)
%function results = animate_series(root,x,y,z, scale)


    
    files = dir(strcat(root,'*.img'));
	if (size(files)==[0 1])
		tdata=0;
        fprintf('%s-----images not found',file);
		return;
	end

	hfiles = dir(strcat(root,'*.hdr'));
	sz = size(files)
	hfiles(1).name;
	h = read_hdr(hfiles(1).name);
    
	if nargin <4
        x=ceil(h.xdim/2);
        y=ceil(h.ydim/2);
        z=ceil(h.zdim/2);
    end
        
    if nargin < 5
        scale  = floor(255/100);
    end
    
    mygray = [1:255]';
    mygray = [mygray mygray mygray]/255;
    colormap(mygray);
    
    cor=moviein(sz);
    sag=moviein(sz);
    tra=moviein(sz);
    
    for count = 1:sz
        a = read_img2(h, files(count).name);
        a = abs(a);
        
%         cor(count) = im2frame(a(:,y,:)*scale +1,mygray);
%         sag(count) = im2frame(a(x,:,:)*scale +1, mygray);
%         tran (count)= im2frame(a(:,:,z)*scale +1, mygray);
        image(a(:,y,:)*scale +1);
        cor(count) = getframe(gcf);
        %pause
        image(a(x,:,:)*scale +1);
        sag(count) = getframe(gcf);
        %pause
        image(a(:,:,z)*scale +1);
        tran (count)= getframe(gcf);
        %pause
        fprintf('\nadding image ...%s',files(count).name);
    end
    
    movie(tran,10,5,[1 1 100 100])
    movie2avi(tra,'transverse.avi')
  

    
   return 
