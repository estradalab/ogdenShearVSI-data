function [series_data, h] = read_img_series(root)
%
% 	[series_data , h] = read_img_series(root)
% 
% Luis hernandez
% last edit 6/26/2004
%
% hdr		header information
% root	root name of the image files
% 
% ouput: 	a matrix with pixels in each column and
% 		 rows means position in time
%
% similar to read_img, but used when there are multiple files
% for the time series
%
% (c) 2005 Luis Hernandez-Garcia 
% University of Michigan
% report bugs to:  hernan@umich.edu
%

% Figure out the direcory and go there
curDir = pwd;
slash = find(root=='/' | root=='\');
if ~isempty(slash)
    dirStr = root(1:slash(end));
    cd (dirStr);
end

names = dir(sprintf('%s*.img', root));
hnames = dir(sprintf('%s*.hdr', root))

if isempty(names)
	fprintf('\n Did not find ANALYZE series - trying NIFTI');	
	names = dir(sprintf('%s*.nii*', root));
	if ~isempty(names)
		h = read_nii_hdr(names(1).name); 
		h = nii2avw_hdr(h);
	else
		fprintf('\n No NIFTI either');
		series_data=0;
		h=0;
		return
	end
else
	h = read_hdr(hnames(1).name);
end

if h.tdim >1
    fprintf('\n   WARNING.  THE FIRST FILE CONTAINS A TIME SERIES');
    fprintf('\n   READING ONLY ONE FILE \n');
    [series_data h] = read_img(names(1).name);
    whos series_data
else

    series_data = zeros(length(names), h.xdim*h.ydim*h.zdim);

    for count=1:length(names)
        fprintf('\r Reading .... %s', names(count).name);
        tmp = read_img(names(count).name);
        series_data(count,:) = tmp(:);
    end
end

% come back from the direcotry
cd(curDir)

return
