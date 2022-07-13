%test navigators;
clear

nonav=varianms('1345'); 
nav=varianms('1345','nav');
yvec=25:55;
xvec=60:120;

%%
for sl=1:12;
    for ndiff=1:4;
        patch=abs(nonav.image(yvec,xvec,sl,ndiff));
        nonavnoise(sl,ndiff) = std(abs(patch(:)))/mean(abs(patch(:)));
        
        patch=abs(nav.image(yvec,xvec,sl,ndiff));
        navnoise(sl,ndiff) = std(abs(patch(:)))/mean(abs(patch(:)));
    end
end

figure;
subplot(1,2,1);
mesh(nonavnoise);
set(gca,'zlim',[0 0.3]);

subplot(1,2,2);
mesh(navnoise);
set(gca,'zlim',[0 0.3]);

colormap gray
%%
figure;
clim=[2100 700 700 700];
for nd=1:4;
    subplot(2,2,nd);
    imagesc(abs(nav.navigators(:,:,nd)));
    set(gca,'clim',[0 clim(nd)]);
end

figure;
for nd=1:4;
    subplot(2,2,nd);
    imagesc(angle(nav.navigators(:,:,nd)));
end