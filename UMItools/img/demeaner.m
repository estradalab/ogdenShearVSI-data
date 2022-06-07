function demeaner(root)
% function demeaner(root)
%
% output:  an image that has been demeaned by its slices
%
% slices can be cruel
%
[data hdr] = read_img(root);
Npix = hdr.xdim *hdr.ydim;
outdata = zeros(size(data));

for sl = 1:hdr.zdim
    slMean = mean ( data(:, (sl-1)* Npix+1 : sl*Npix) , 2) ;
    outdata( : , (sl-1)* Npix +1 : sl*Npix )  = ...
        data( : , (sl-1)* Npix+1 : sl*Npix ) ...
        - repmat(slMean, 1, Npix);
end

outdata = outdata + 100;
write_img('demeaned.img', outdata, hdr);
return
