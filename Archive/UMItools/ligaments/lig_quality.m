function out=lig_quality(DESTE) 
%quality script to assess meaning of angles
% DESTE is a structure produced in lid_deste_3d
% at this point it contains smoothed complex data for the three
% encoding directions, as well as axes, and a mask


m=mean(mean(abs(alldata),5),4);
d=diff(abs(alldata),1,4);

for nd=1:ND
	t=squeeze(d1048.image(:,:,sli,:));	%the image on which to compare stuff
	
	r=abs(t(:,:,1))./abs(t(:,:,2));		%the ratio of intensities
	d=abs(t(:,:,1))-abs(t(:,:,2));		%the difference of intensities
	m=mean(abs(t),3);					%the mean intensity
	p=angle(t(:,:,1)./t(:,:,2));
	
	figure('colormap',jet(64));
	subplot(1,4,1);
	imagesc(r,[0.6 1.4]);
	title('r');
	colorbar;
	
	subplot(1,4,2)
	imagesc(log10(m),[2 4]);
	title('log10(m)');
	colorbar;
	
	
	subplot(1,4,3)
	imagesc(d,[-500 500]);
	colorbar;
	title('d');
	
	subplot(1,4,4)
	imagesc(p,[-pi pi]);
	colorbar;
	title('p');
end

%%
sems=lig_sems(1237,'gradlag',40);
plotms(angle(sems.orig_data.image(:,:,28:2:50)),'col');
%plotms(log10(sems.orig_data.image(:,:,28:2:50)/200),'clim',[0.5 2.5]);
