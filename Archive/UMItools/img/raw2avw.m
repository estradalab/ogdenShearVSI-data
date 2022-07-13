function raw2avw(raw,name)
h = define_avw_hdr;
h.xdim=size(raw,1);
h.ydim=size(raw,2);
h.zdim=size(raw,3);
h.tdim=size(raw,4);

h.xsize=1;
h.ysize=1;
h.zsize=1;
h.tsize=1;

h.datatype=4;
h.bits=16;

raw = raw(:);
raw=reshape(raw, h.xdim*h.ydim*h.zdim, h.tdim);

write_img(name, raw', h);

return
