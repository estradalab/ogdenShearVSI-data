function XYZ_print(fn,xyz_loc)
fp = fopen(fn,'w');
fprintf(fp,'%5d %5d %5d\n',xyz_loc'); 
fclose(fp)
