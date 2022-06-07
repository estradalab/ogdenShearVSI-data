function normImg(imgname, option, sfactor)
% function normImg(imgname, option)
%
% normalize to:
% 1 - the total energy of the image
% 2 - the maximum of the image
% 3 - an arbitrary value
% 4 - rank based on histogram with 100 bins

[img, h] = read_img(imgname);

switch(option)
	case 1
	outimg = img/sum(img(find(img)));

	case 2
	outimg = img /max(img(find(img(:))));

	case 3
	outimg = img * sfactor;

	case 4
	cleanimg = zeros(size(img));
	cleanimg(find(img)) = img(find(img));
	range = (max(cleanimg) - min(cleanimg)) ;
	binsize = range / 100;
	[N, X] = hist(cleanimg,100);
	outimg = zeros(size(img));
	for p = 1: length(outimg)
		val = cleanimg(p);
		if cleanimg(p)~=0
			outimg(p) = find( abs(X-val) < binsize/2 );
		end
	end
end

write_img(sprintf('N_%s', imgname) , outimg,  h);

return


