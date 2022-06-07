function niih = avw2nii_hdr(hdr)
%funtion avwh = avw2nii_hdr(hdr)

niih = define_nii_hdr;

niih.sizeof_hdr = hdr.sizeof_hdr ;
%hdr.pad1 = niih.?  
niih.extents = hdr.extents;
%hdr.pad2  = niih.?
%hdr.regular= niih.?
%hdr.pad3  = niih.?
if hdr.tdim>1
    niih.dim(1)= 4 ;
else
    niih.dim(1)= 3 ;
end
niih.dim(2) = hdr.xdim ;
niih.dim(3) = hdr.ydim ;
niih.dim(4) = hdr.zdim ;
niih.dim(5) = hdr.tdim ;
%hdr.pad4 = niih.?
niih.datatype = hdr.datatype ;
niih.bitpix = hdr.bits ;
%hdr.pad5 = niih.?
niih.pixdim(2) = hdr.xsize ;
niih.pixdim(3) = hdr.ysize ;
niih.pixdim(4)= hdr.zsize ;
%hdr.pad6 = niih.?
niih.cal_max = hdr.glmax ;
niih.cal_min = hdr.glmin ;
niih.descrip = hdr.descrip ;
niih.aux_file = hdr.aux_file ;
%hdr.orient = niih.?
%hdr.origin = ?
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
