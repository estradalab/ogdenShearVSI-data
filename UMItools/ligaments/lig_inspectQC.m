function fh = lig_inspectQC(in)
% out=lig_inspectQC(in)
% is is either a structure containing a DESTE workspace, 
% or the filename of a DESTE workspace with element QC, or a QC structure; 
% anything goes, almost

if ischar(in)
	load(in,'QC');
elseif isstruct(in)
	if isfield(in,'QC')
		QC=in.QC;
	elseif isfield(in,'data')
		QC=in;
	end
end

[sig,noise]=estimate_snr(QC.data);

hc=0;

polstr='+-';
for encdir=1:3
	for polarity=1:2
		plotms(QC.data(:,:,:,polarity,encdir),'layout',[2 32],'clim',[noise 2*sig],'col');
		
		titlestring=['Magnitude: Encode direction= ' num2str(encdir) ' ,  encode polarity= ' polstr(polarity)];
		set(gcf,'Name',titlestring);
		hc=hc+1;
		fh(hc)=gcf;
	end
end

for encdir=1:3
	for polarity=1:2
		plotms(angle(QC.data(:,:,:,polarity,encdir)),'layout',[2 32],'col');
		
		titlestring=['Phase: Encode direction= ' num2str(encdir) ' ,  encode polarity= ' polstr(polarity)];
		set(gcf,'Name',titlestring);
		hc=hc+1;
		fh(hc)=gcf;
	end
end


end

