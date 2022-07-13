function brain_synth_flair(infit,taxis,TR)
%calculate synthetic flair signal

numti=numel(taxis);

figure;
colormap gray;
for jj=1:numti;
    
    synthimage=infit.amplitude.*(1-2.*exp(-taxis(jj)./infit.T1)+exp(-TR./infit.T1));
    imagesc(synthimage);
    set(gca,'clim',[-5000 5000]);
    pause(0.1);
end