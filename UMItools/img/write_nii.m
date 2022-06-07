function  write_nii(name, data, h, doCompress)

% Luis hernandez
% last edit 9-26-2006
%
% function  write_nii(name, data, h, doCompress)
%
% Writes the data to a NIFTI format file  name  containing image data 
% this also handles a timeseries in one file (ie -each image is a row )
%
%
% (c) 2006 Luis Hernandez-Garcia 
% University of Michigan
% report bugs to:  hernan@umich.edu
%



   [pFile,messg] = fopen(name,  'wb','ieee-le' );
   if pFile == -1
      errormesg(messg);   
      return;
   end
   
      
fwrite(pFile,  h.sizeof_hdr, 'int32' );      % should be 348!
fwrite(pFile,  h. data_type      ,  'ubit8'  );
fwrite(pFile,  h. db_name        ,  'ubit8'  );
fwrite(pFile,  h. extents        , 'int32' ) ;
fwrite(pFile,  h. session_error  , 'int16' ) ;
fwrite(pFile,  h. regular        , 'ubit8' ) ;
fwrite(pFile,  h. dim_info       , 'ubit8' ) ;
fwrite(pFile,  h.dim         ,'int16' ) ;
fwrite(pFile,  h.intent_p1   ,'float32' ) ;
fwrite(pFile,  h.intent_p2   ,'float32' ) ;
fwrite(pFile,  h.intent_p3   , 'float32' ) ;
fwrite(pFile,  h.intent_code  ,'int16' ) ;
fwrite(pFile,  h.datatype    ,'int16' ) ;
fwrite(pFile,  h.bitpix      ,'int16' ) ;
fwrite(pFile,  h.slice_start ,'int16' ) ;
fwrite(pFile,  h.pixdim      ,'float32' ) ;
fwrite(pFile,  h.vox_offset  ,'float32' ) ;
fwrite(pFile,  h.scl_slope   ,'float32' ) ;
fwrite(pFile,  h.scl_inter   ,'float32' ) ;
fwrite(pFile,  h.slice_end   ,'int16' ) ;
fwrite(pFile,  h.slice_code  ,'ubit8' ) ;
fwrite(pFile,  h.xyzt_units  ,'ubit8' ) ;
fwrite(pFile,  h.cal_max     ,'float32' ) ;
fwrite(pFile,  h.cal_min     ,'float32' ) ; 
fwrite(pFile,  h.slice_duration  ,'float32' ) ;
fwrite(pFile,  h.toffset     ,'float32' ) ;
fwrite(pFile,  h.glmax       ,'int32' ) ;
fwrite(pFile,  h.glmin       ,'int32' ) ;
fwrite(pFile,  h.descrip      ,'ubit8'  );
fwrite(pFile,  h.aux_file     , 'ubit8');
fwrite(pFile,  h.qform_code   ,'int16' ) ;
fwrite(pFile,  h.sform_code   ,'int16' ) ;
fwrite(pFile,  h.quatern_b    ,'float32' ) ;
fwrite(pFile,  h.quatern_c    ,'float32' ) ;
fwrite(pFile,  h.quatern_d    ,'float32' ) ;
fwrite(pFile,  h.qoffset_x    ,'float32' ) ;
fwrite(pFile,  h.qoffset_y    ,'float32' ) ;
fwrite(pFile,  h.qoffset_z    ,'float32' ) ;
fwrite(pFile,  h.srow_x       ,'float32' ) ;
fwrite(pFile,  h.srow_y       ,'float32' ) ;
fwrite(pFile,  h.srow_z       ,'float32' ) ;
fwrite(pFile,  h.intent_name  ,'ubit8'  );
fwrite(pFile,  h.magic        ,'ubit8' );
if length(h.magic)==3
	fwrite(pFile, 0         ,'ubit8' );
end

%fwrite(pFile,  h.originator   ,'int16' );
%fwrite(pFile,  h.esize        ,'ubit8');
%fwrite(pFile,  h.ecode        , 'ubit8');
%fwrite(pFile,  h.edata        , 'ubit8');   

fwrite(pFile,  0.0        , 'float32');   

switch h.datatype     
case 2
      fmt =  'ubit8' ;
case 4
      fmt =  'short' ;
case 8
      fmt =  'int' ;
case 16
      fmt =  'float' ;
case 32
      fmt =  'float' ;
           
otherwise
      errormesg(sprintf('Data Type %d Unsupported. Aborting' ,  h.datatype));
      return
end

% a kludge to make sure we can go to the right offset in the file
fwrite(pFile, '             ','ubit8');
fseek(pFile,h.vox_offset,'bof');
if h.dim(5) >1
     fwrite(pFile, data', fmt);
     else
     fwrite(pFile, data, fmt);
end

fclose(pFile);
if (doCompress)
     str = sprintf('!gzip %s',name)
     eval(str)
end  

 return

