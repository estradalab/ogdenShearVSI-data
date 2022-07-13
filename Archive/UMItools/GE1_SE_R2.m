%survey

timestamp='1229';

data=varianms([timestamp '_dwgesege_us01.fid'],'zf','shphase',0.1,'flipud');
data=hardwarephase(data);
data=softshim(data);
data=blur3d(data,'thk',1); %1mm virtual slice thickness

imagesize=[size(data.image) 1 1];
nechos=imagesize(4);
nslices=imagesize(3);
%Generate R2star, and a spin echo amplitude;
[out.SEamplitude, out.R2star]=matrix_expfit(data.image(:,:,:,2:nechos), data.pars.te00 * [0: (nechos-2)]);
out.T1amplitude=exp(+min(out.R2star,2/data.pars.te00)*data.pars.te0).* squeeze(data.image(:,:,:,1));



