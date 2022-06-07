!rm movie.*
!rm base.*


h = read_hdr('fmri002.hdr');
baseline = read_img_slice(h,'fmri049.img',7);
b = baseline.data;

hm = read_hdr('spm.hdr');
mask = read_img_slice(hm,'spm.img',7);
m =mask.data;


for i=1:hm.xdim
   for j=1:hm.ydim
      if (m(i,j)~=0 )
         m(i,j) = 1;
      end
   end
end

   


outf1 = fopen('movie.img','a');
outf2 = fopen('base.img','a');


for i=2:61
   
   if i<10
      name = strcat('sfmri00',num2str(i));
   
   else
      name = strcat('sfmri0',num2str(i));
   end
   
   name = strcat(name,'.img');
   
   disp(name)
   
   signal = read_img_slice(h, name,7);
   s = signal.data;
   
   r = (s - b + 700 ).*m *10;
  
   
   max(max(r))
   min(min(r))
   
   fwrite(outf1,r', 'int16');
   fwrite(outf2,b', 'int16');

end



fclose(outf1);
fclose(outf2);


outh = h;
outh.zdim = 60;
write_hdr('movie.hdr',outh);
!cp movie.hdr base.hdr
disp('done.')





