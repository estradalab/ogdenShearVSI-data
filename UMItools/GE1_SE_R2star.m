function out=GE1_SE_R2star(timestamp)

data=varianms([timestamp '_dwgesege_us01.fid'],'zf','shphase',0.1,'flipud');
data=hardwarephase(data);
data=softshim(data);
data=blur3d(data,'thk',1); %1mm virtual slice thickness

imagesize=[size(data.image) 1 1];
nechos=imagesize(4);
nslices=imagesize(3);
%Generate R2star, and a spin echo amplitude;
[out.SEamplitude, out.T2star, out.R2star]=matrix_expfit(data.image(:,:,:,2:nechos), data.pars.te00 * (0: (nechos-2)));

%generate a weighted sum of the last three echoes. for averaging and nopise
%reduction. The effective TE for this weighted sum is TE_eff = TE + TE00,
%where TE00 = TE2 is the spacing betwenn the last three echoes
%1. generate the echo weighting factor
echofactor=ones(size(data.image(:,:,:,2:nechos)));
for jj=1:nechos-1;
    echofactor(:,:,:,jj)=exp(-abs(out.R2star)*(jj-1)*data.pars.te00);
end
out.SEbar = sum(abs(data.image(:,:,:,2:nechos)).*echofactor,4)./sum(echofactor,4);

%out.GE_Tzero=exp(+min(out.R2star,1/data.pars.te00)*data.pars.te0).* squeeze(data.image(:,:,:,1));  %T2star correction of the first echo, to be checked
%2. put the first echo into a field
out.GE_Tzero= squeeze(data.image(:,:,:,1));
out.image=data.image;


out.pars=data.pars;



