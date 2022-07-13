%survey the effect of different blurring on ge3d image quality
close all
clear

%% read and blur;
wksp=load('DESTE_strains_1059_1115_1131_1037.mat','origdata','origpars');
pic=wksp.origdata(:,:,:,3);  %take the RO encode data
pars=wksp.origpars;

voxvec={[0 0 0]; [1 1 1]; [2 2 2]; [3 3 3]};

si=size(pic);
sb=numel(voxvec);
blurredata=zeros([si sb]);

for nt=1:sb; 
    blurreddata(:,:,:,nt)=blur3d(pic,'vox',voxvec{nt},'grid',[128 64 64]);
    [signallevel(nt), noiselevel(nt)] = estimate_noiselevel(blurreddata(:,:,:,nt));
end


%% plot blurred
for ii=5:5:125;
    figure('position',[50 100 1900 500]);
    for nt=1:sb;
        subplot(2,sb,nt);
        imagesc(squeeze(abs(blurreddata(ii,:,:,nt)))',[0 6000]);
        title(num2str(voxvec{nt},'%d'));
        axis image;
        
        subplot(2,sb,nt+sb);
        imagesc(squeeze(angle(blurreddata(ii,:,:,nt))+pi)');
        title(num2str(voxvec{nt},'%d'));
        axis image;
    end
end