function Tmapper(rootname)
% function Tmapper(rootname)
%
%      Luis Hernandez-Garcia at U Michigan
%      12-17-2003
%
% Tmap calculator for SPM images avoiding the mask issues in SPM99
%
% this function takes a set of images with a common rootname
% and calculates:
%
% df.img:
%    the number of available measurements at each pixel
%    (not NaN's or Infs): this is the  degrees of freedom 
%    at that pixel
% mean.img:
%     the mean of the pixels excluding NaNs and Infs
% var.img:
%     the variance of the pixels exluding NaNs and Infs
% tscore.img:
%     the t score wherever possible, taking into account the 
%     available number of measurements.  
%     If there is only one measurement, we call it a NaN.  
%     (not much of a  t-test if there is only one measurement.
% pval.img
%     the corresponding p values to the tscores using the
%     degrees of freedom calculated earlier
%
% 10lnP:  
%     this one is just -10*ln(p)  for bettern display purposes.
%
% USAGE example:  Tmapper('rfx_')
%
% (when doing RFX, it's useful to put links to the individual
% contrast images in the same directroty)
%
% warning:  it's slow! 
%

warning off

names = dir(sprintf('%s*.hdr',rootname));
Nfiles = size(names,1);

raw = [];
for count=1:Nfiles
	
	h = read_hdr(names(count).name);
	str = names(count).name;
	str(length(str)-3:end) = '.img';
	fprintf('\rreading ...%s', str);
	raw = [raw ; read_img(h, str)];

end

whos raw
Npix = size(raw,2);

df=zeros(Npix,1);
ave=zeros(Npix,1);
vr=zeros(Npix,1);
tscor=zeros(Npix,1);
pval=zeros(Npix,1);
%tic

for count = 1:Npix
	tmp = raw(:,count);
	tmp = tmp(find(tmp));
	tmp = tmp(find(~isnan(tmp)));	
	tmp = tmp(find(~isinf(tmp)));

	if ~isempty(tmp)
		n = length(tmp);
		a = mean(tmp);
		v = var(tmp);
		if n==1
			t=NaN;
			p=NaN;
		else
			t = a/sqrt(v/n);
			p = 1-tcdf(t,n-1);
		end
	else
		n = NaN;
		a = NaN;
		v = NaN;
		t = NaN;
		p = NaN;
	end
	if (rem(count,1000)==0)
		%toc
		%tic
		fprintf('\r*-%d', count);
	end
	df(count) =  n-1;
	ave(count) = a;
	vr(count) =  v;
	tscor(count) = t;
	pval(count)  = p;
%keyboard
end

h.tdim = 1;

write_hdr('tscore.hdr',h);
write_img_data('tscore.img',tscor,h);

write_hdr('pval.hdr',h);
write_img_data('pval.img',pval,h);

write_hdr('10lnP.hdr',h);
write_img_data('10lnP.img',-10*log(pval),h);

write_hdr('var.hdr',h);
write_img_data('var.img',vr,h);

write_hdr('mean.hdr',h);
write_img_data('mean.img',ave,h);

write_hdr('df.hdr',h);
write_img_data('df.img',df,h);


return

