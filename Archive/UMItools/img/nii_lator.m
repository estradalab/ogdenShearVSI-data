% this script converts all the files in this directory 
% from AVW to NII format.

fnames = dir('*.img');

for n=1:length(fnames)
	[d h] = read_img(fnames(n).name);
	niih = avw2nii_hdr(h);
	
	outname = fnames(n).name;
        outname(end-3:end)= '.nii';

        fprintf('\rWriting out %s', outname);
	write_nii(outname, d, niih,0);

	str = ['!rm ' fnames(n).name];
	eval(str);
	hname = fnames(n).name;
        hname(end-3:end)= '.hdr';
	str = ['!rm ' hname];
	eval(str)
end

