function um7t_stab01(pathname)
% function um7t_stab01(pathname)
% e.g. um7_stab01('~vnmr1/vnmrsys/QA_data/s_20111110401/epip03.img/')
% 

spike_thresh = 1 ;


% first read in the file and translate it into NII format
[data h ] = fdf2nii(pathname, 0);

data = (reshape(data,  h.dim(2)*h.dim(3)*h.dim(4), h.dim(5)))';

[x,y,z] = ndgrid([31:33],[31,33],[4:6]) ;

inds = sub2ind( [h.dim(2), h.dim(3), h.dim(4)], x(:), y(:), z(:));

% test # 1 - look for spikes in the time series:
vdata = var(data);
mdata = mean(data,2);
snr_map = mdata ./ sqrt(vdata);
ddata = zeros(size(data));

for t=2:size(data,1)
	ddata(t,:) = data(t, :) - data(t-1,:);
	
	tmp1 = reshape(ddata(t,:), h.dim(2), h.dim(3), h.dim(4));
	tmp2 = reshape(data(t,:), h.dim(2), h.dim(3), h.dim(4));

	figure(33)
	subplot(322)
	lightbox(tmp1,[],3);
    title(['Difference image ' num2str(t)])
	
	subplot(321)
	lightbox(tmp2, [], 3);
    title(['Image ' num2str(t)])
	drawnow
    pause(0.1)
end

stab_plot = data(:, inds);

x=[0:length(stab_plot)-1]';
y = mean(stab_plot,2);
stdDev = std(y);

[p s] = polyfit( x,y , 3); 
y2 = polyval(p,x);

subplot(323)
plot(stab_plot)
hold on
plot(x,y2,'k') 
title('ROI and fitted mean of the ROI')


subplot(326)
stab_plot = stab_plot - repmat(mean(stab_plot), size(stab_plot,1),1);

%lightbox(reshape(vdata, h.dim(2), h.dim(3), h.dim(4)), [],3);
%title('The Variance Map')
lightbox(reshape(snr_map, h.dim(2), h.dim(3), h.dim(4)), [],3);
title('SNR Map')

mstab_plot = mean(stab_plot,2);

fstab_plot = fftshift(fft(fftshift(mstab_plot)));
subplot(325)
w = linspace(-1,1,length(fstab_plot));
plot(w, abs(fstab_plot))
title('FFT of the ROI')



subplot(324)
cla
axis off
text(0, 0.8,['standard deviation = ' num2str(stdDev,3)])
text(0, 0.6,['poly coeffs = ' num2str(p,3)])
text(0, 0.4, date)

title('Stats')

str=['print -depsc ' date '_stability']
 eval(str)


