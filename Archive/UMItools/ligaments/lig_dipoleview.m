function out=lig_dipoleview(D,slicenumber)
% side by side plot of various aspect of the data, in a slice

clim_log=[min(log(abs(D.data(:)))) max(log(abs(D.data(:))))];

for jj=1:numel(slicenumber)
    sl=slicenumber(jj);
    
    axis1=D.axis1;
    axis2=D.axis2;
    
    figure('Position', [300 300 1800 700]);
    
    subplot(1,8,1);
    imagesc(axis2,axis1,log(abs(D.data(:,:,sl,1))),clim_log);
    axis image;
    grid on;
    set(gca,'YDir', 'normal');
    title('log(|STE|)_1');
    
    subplot(1,8,2);
    imagesc(axis2,axis1,log(abs(D.data(:,:,sl,2))),clim_log);
    axis image;
    grid on;
    set(gca,'YDir', 'normal');
    title('log(|STE|)_2');
    
    subplot(1,8,3);
    imagesc(axis2,axis1,log(abs(D.data(:,:,sl,3))),clim_log);
    axis image;
    grid on;
    set(gca,'YDir', 'normal');
    title('log(|STE|)_3');
    
    
    subplot(1,8,4);
    imagesc(axis2,axis1,log(abs(D.HIRES.magnitude(:,:,sl))));
    axis image;
    grid on;
    set(gca,'YDir', 'normal');
    title('log(|HIRES|)');
    
    for kk=1:3
        subplot(1,8,4+kk);
        imagesc(axis2,axis1,D.HIRES.magnitude(:,:,sl));
        axis image;
        grid on;
        set(gca,'YDir', 'normal');
        lig_heatoverlay(D.Lagrange(:,:,sl,kk,kk));
        axis image;
        grid on;
        set(gca,'YDir', 'normal');
        title(['E_{' num2str(kk*11) '}']);
    end
    
    subplot(1,8,8);
    imagesc(axis2,axis1,std(log(abs(D.data(:,:,sl,:))),0,4));
    set(gca,'colormap', gray);
    axis image;
    grid on;
    set(gca,'YDir', 'normal');
    title('std(log(|STE|))');
    
end
    