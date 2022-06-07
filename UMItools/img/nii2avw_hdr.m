function hdr = nii2avw_hdr(niih)
%funtion avwh = nii2avw_hdr(niih)

hdr = define_avw_hdr;

hdr.sizeof_hdr = niih.sizeof_hdr;
%hdr.pad1 = niih.?  
hdr.extents= niih.extents;
%hdr.pad2  = niih.?
%hdr.regular= niih.?
%hdr.pad3  = niih.?
hdr.dims = niih.dim(1);
hdr.xdim = niih.dim(2);
hdr.ydim = niih.dim(3);
hdr.zdim = niih.dim(4);
hdr.tdim = niih.dim(5);
%hdr.pad4 = niih.?
hdr.datatype = niih.datatype;
hdr.bits  = niih.bitpix;
%hdr.pad5 = niih.?
hdr.xsize = niih.pixdim(2);
hdr.ysize = niih.pixdim(3);
hdr.zsize = niih.pixdim(4);
%hdr.pad6 = niih.?
hdr.glmax = niih.cal_max;
hdr.glmin = niih.cal_min;
hdr.descrip = niih.descrip;
hdr.aux_file = niih.aux_file;
%hdr.orient = niih.?
hdr.origin(1) = niih.qoffset_x / hdr.xsize ; 
hdr.origin(2) = -niih.qoffset_y / hdr.ysize; 
hdr.origin(3) = -niih.qoffset_z / hdr.zsize; 
%hdr.generated = niih.?
%hdr.scannum = niih.?
%hdr.patient_id = niih.?
%hdr.exp_date = niih.?
%hdr.exp_time = niih.?
%hdr.hist_un0 = niih.?
%hdr.views = niih.?
%hdr.vols_added = niih.?
%hdr.start_field = niih.?
%hdr.field_skip = niih.?
%hdr.omax = niih.?
%hdr.omin = niih.?
%hdr.smax = niih.?
%hdr.smin = niih.?

return
