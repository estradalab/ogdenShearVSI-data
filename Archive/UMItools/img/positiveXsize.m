subjs = [
'060518ak'
'060518mw'
'060524aw'
'060524mm'
'060524ms'
'060606jm'
'060613cg'
'060613cj'
'060620ag'
'060620cs'
];

rootdir = '/export/subjects/'
for s=1:size(subjs,1)
	cd ([rootdir subjs(s,:) '/anatomy/'])
	files = dir('t1*.hdr');
	for f=1:length(files)
		fprintf('\n file ... %s', files(f).name)
		str = sprintf('!cp %s e%s', char(files(f).name), char(files(f).name) )
		eval(str)
		str = sprintf('!cp %s he%s', char(files(f).name), char(files(f).name) )
		eval(str)
		str = '!rm *_sn.mat whet*.img whet*.hdr t1*.mat et1*.mat het1*.mat '
		eval(str)
	end
	% now let's do the funcs
	cd ([rootdir subjs(s,:) '/func/run_01/img'])
	files = dir('*.hdr')
	for f=1:length(files)
		h = read_hdr(files(f).name);
		h.xsize = abs(h.xsize);
		write_hdr(files(f).name , h);
	end

	cd ([rootdir subjs(s,:) '/func/run_01/ra_img'])
	files = dir('*.hdr')
	for f=1:length(files)
		h = read_hdr(files(f).name);
		h.xsize = abs(h.xsize);
		write_hdr(files(f).name , h);
	end
	pwd
end


