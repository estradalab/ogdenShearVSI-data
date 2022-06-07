%
% (c) 2005 Luis Hernandez-Garcia 
% University of Michigan
% report bugs to:  hernan@umich.edu
%

h=read_hdr('avol_0001.hdr');
im=read_img2(h,'avol_0001.img');
show(im(:,:,3));


m=zeros(100,1);
%m(5:8)=.8;
%m(15:18) = 0.8;
%m(30:33) = 0.8
m(1:33)= spm_hrf(1);
m(31:63) =spm_hrf(1);
m(61:93) =spm_hrf(1);
m = (m+10)/10;

plot(m)

for i=1:size(m,1)
	
      im2 = im;
      im2(h.xdim*0.25:h.xdim*0.75, h.ydim*0.25:h.ydim*0.75, h.zdim/2:end) = ...
          im2(h.xdim*0.25:h.xdim*0.75, h.ydim*0.25:h.ydim*0.75, h.zdim/2:end) * m(i);
      str=sprintf('vol_%04d.img',i);
      write_img_data(str,im2,h);
      
      str=sprintf('vol_%04d.hdr',i);
      write_hdr(str,h);
      
      
end

return

datasize=15
x=30; y=30; z=3;

index= h.xdim*h.ydim* (z-1) + h.xdim*(y-1) + x-1 ;
for i=1:datasize
   str=sprintf('fmri%04d.img',i)
   pFile = fopen(str,'rb');
   fseek(pFile,index*2,'bof');
   series(i)=fread(pFile,1,'int16');
   fclose(pFile);
end

% testing the Tmapper
h=read_hdr('tsamp.hdr');
im=read_img2(h,'tsamp.img');
show(im(:,:,3));


m=rand(10,1);
plot(m)

for i=1:size(m,1)
	
      im2 = im;
      im2(h.xdim*0.25:h.xdim*0.75, h.ydim*0.25:h.ydim*0.75, 6:10) = im2(h.xdim*0.25:h.xdim*0.75, h.ydim*0.25:h.ydim*0.75, 6:10) * m(i) + 10;
      str=sprintf('asamp_%04d.img',i);
      write_img_data(str,im2,h);
      
      str=sprintf('asamp_%04d.hdr',i);
      write_hdr(str,h);
      
      im2(h.xdim*0.25:h.xdim*0.75, h.ydim*0.25:h.ydim*0.75, 6:10) = im2(h.xdim*0.25:h.xdim*0.75, h.ydim*0.25:h.ydim*0.75, 6:10) * m(i) + 15;
      str=sprintf('bsamp_%04d.img',i);
      write_img_data(str,im2,h);
      
      str=sprintf('bsamp_%04d.hdr',i);
      write_hdr(str,h);
      
end
