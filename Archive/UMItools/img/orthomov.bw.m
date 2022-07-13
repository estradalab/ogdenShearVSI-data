function result = orthomov( roi, root, x, y, z)
% function result = orthomov(roi_size, rootname, x, y, z)
%
% this function thakes a file in analyze format and 
% displays orthogonal sections thru the planes interscting
% at x,y,z
global wscale

origFig = figure;
movieFig = figure;
figure(origFig);
colormap(gray)

nargin = 0;
if nargin==0
  roi=0;
end
       
if nargin < 2
  [fname path] = uigetfile('*.img','Select  HIGHEST NUMBERED *.img file in time series');
  name = strcat(path,fname)
 
  sz = size(name);
   
  imgname = strcat(name(1,1:sz(2)-4) , '.img');
  hdrname = strcat(name(1,1:sz(2)-4) , '.hdr');
  h = read_hdr(hdrname);
  d = read_img2(h, imgname);
else 
  %figure out the whole name of the header and image files
  % and read them.  Note the use of strcat
  h = read_hdr(strcat(root,'.hdr'));
  d = read_img2(h, (strcat(root,'.img') ));
end

if nargin < 4
  x=h.xdim/2;
  y=h.ydim/2;
  z=h.zdim/2;
end

x=ceil(x); 
y=ceil(y); 
z=ceil(z);

i=0;
% configure colormap
my_map=(0:255)';
my_map=[my_map my_map my_map]/256;
colormap(my_map);

%gather all of the images in a time-series
numImgs = str2num(fname(end-7:end-4));  %img names start at 0001
imgSize = size(d);
imgSeries = zeros([imgSize numImgs]);

for i = 1:numImgs
  j = i - 1;    %img names start at 0000
  j = i;        %img names start at 0001
  imgname = strcat(name(1,1:sz(2)-8),sprintf('%0.4d',j),'.img');
  d = read_img2(h,imgname);
  
  % scale image to fit colormap
  range= max(max(max(d))) - min(min(min(d)));
  dd = (d-min(min(min(d))))*256/range;
  if ~isempty(wscale)
    dd = (d-wscale(1))*256 / wscale(2);
  end

  imgSeries(:,:,:,i) = dd;
end

subplot(224), hist(d(:),100)
while i >= -10
  [fig1, fig2, fig3] = ov(h,dd,x,y,z, roi);
  if roi>0
    tmp=mean(mean(d(x-roi:x+roi, y-roi:y+roi,z-roi:z+roi)));
  else
    tmp=d(x,y,z);
  end
    
  str = sprintf('\n(x,y,z)=  (%d %d %d) , val= %6.2f  \n', x, y, z, tmp);
  subplot(221), title(str) 

  [i j button] = ginput(1);
  
  i=round(i);
  j=round(j);
  fig = floor(gca);
        
  if button == 1
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
	
      %exiting the program when you click outside bounds
      if (i<=-1)   
        colordef white
        return
      end
  elseif button == 3
      %play movie of selected view
      figure(movieFig)
      aviF = avifile('ortho.avi')
      clear M;
      for idx = 1:numImgs
        dd = imgSeries(:,:,:,idx);
        [fig1, fig2, fig3] = ovmov(h,dd,x,y,z, roi, my_map);
        M(idx) = getframe(gcf);
	aviF=addframe(aviF, M(idx));
	fprintf('Frame : %d  ',idx);
      end
      aviF=close(aviF);
      %fps = input('Enter Frame Rate (in FPS): ');
      fps = 4;
      movie(movieFig,M,1,fps,[0 0 0 0]);
      figure(origFig);
  end
end

return
