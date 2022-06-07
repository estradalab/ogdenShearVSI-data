% example :

files = [ ...
{'/export/subjects/060405lh/results/con_0001.hdr'}
{'/export/subjects/060405lh/anatomy/T1/et1spgr.hdr'}
{'/export/subjects/060406bh/func/run_01/ra_img/ravol_0001.img'}
% ... put your file names here....

];

for f=1:length(files)
	swap_hdrx(files(f))
end


