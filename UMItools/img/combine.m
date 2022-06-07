!rm combo.*
   
outf1 = fopen('combo.img','a');


for i=1:60

h = read_hdr('T1stack0.hdr');
base = read_img_slice(h,'T1stack0.img',i);
b = base.data;

hm = read_hdr('movie0.hdr');
spm = read_img_slice(hm,'movie0.img',i);
m =spm.data;


for i=1: hm.xdim
   
   for j=1: hm.ydim
      
      if (m(i,j) ~= 128 )
         
         b(i,j) = m(i,j);
         
      end
      
   end
end


   
   
 fwrite(outf1, b', 'short');

end



fclose(outf1);


outh = hm;
outh.glmin = 0;
write_hdr( 'combo.hdr', outh)


disp('done.')






