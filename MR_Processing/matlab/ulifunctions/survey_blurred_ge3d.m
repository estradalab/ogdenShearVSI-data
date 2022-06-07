%survey the effect of different blurring on ge3d image quality
close all
clear


%% read and blur;
data=varianms(1037,'apod');
pic=data.image;

voxvec={[0 0 0]; [1 2 2]; [2 3 3]; [3 3 3]};

si=size(pic);
sb=numel(voxvec);
blurredata=zeros([si sb]);

for nt=1:sb; 
    blurreddata(:,:,:,nt)=blur3d(pic,'vox',voxvec{nt});
    [signallevel(nt), noiselevel(nt)] = estimate_noiselevel(blurreddata(:,:,:,nt));
    
    blurreddata(:,:,:,nt)=log10(abs(blurreddata(:,:,:,nt)./noiselevel(nt)));
end


%% plot blurred
for ii=50:5:150;
    h=plotms(squeeze(blurreddata(ii,:,:,:)),'layout',[1 sb],'col','clim',[0 3]);
    
    panelaxes=get(h,'children');
    panelaxes=panelaxes(end:-1:1);  %reversethe order so axes are first generated to alst generated
    
    for jj=1:numel(panelaxes);
        axes(panelaxes(jj));
        title(num2str(voxvec{jj},'%d'));
    end
    
    set(gcf,'position', [50 100 1900 500]);
end