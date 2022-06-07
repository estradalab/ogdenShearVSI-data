function sliceTimer(root, TR, Ttag)
% interpolates the data to the center of the TR
% Assumes The first slice is acquired at TR+Tag
%

[data h] = read_nii_img(root);
Nslices = h.dim(4);
Nframes = size(data,1);
Npix = h.dim(2) * h.dimi(3);
Tslice = (TR-Ttag)/Nslices;
TRtimes = TR * ([1:Nframes] -0.5 );
outdata = zeros(size(data));

for sl = 1:Nslices
	fprintf('\rslice timing : ... %d', sl);
	AQtimes= TR*(0:Nframes-1) + Ttag + (sl-1)*Tslice;
	for p=1:Npix
		ts = data(:, Npix*(sl-1) + p);	
		outTS = interp1(AQtimes, ts, TRtimes,'sinc');
		outdata(:, Npix*(sl-1) + p) = outTS;	
	end

end
warning off
h.dim(5) = size(outdata,1);
write_nii(['a' root ], outdata, h)
%keyboard
warning on
return
