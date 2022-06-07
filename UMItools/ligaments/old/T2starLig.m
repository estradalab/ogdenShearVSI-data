function out=T2starLig(timestamp)
%takes a mge3d data set of a ligament, data id bu time stamp

%1. read the data from FID
data=varianms(timestamp,'noiso','zf2','zf3');  %noiso suppresses the automatic isotropy on the first two dimensions
                                              %flags zf and zf3 fill by 1.
                                              %in the three dimensions
                                              
%2. apply a 3 d blur to get the voxel size right
vox=0.1*[1 1 1];
data=blur3d(data,'vox',vox);
pic=data.image;

%3. calculate the time axis
taxis=data.pars.te+[0:data.pars.ne-1]*data.pars.te2;

%4. pick the region of interest where the relaxation times and amplitudes
%are to be calculated use the central fifth of the readout direction
si=size(pic);
xaxis=(1:si(3))/si(3)*data.pars.lpe2*10; %mm
yaxis=(1:si(2))/si(2)*data.pars.lpe*10;  %mm
zaxis=(1:si(1))/si(1)*data.pars.lro*10;  %mm

centvec=round((2*si(1)/5) : (3*si(1)/5));
test=squeeze(sum(squeeze(sum(abs(pic(centvec,:,:,:)),4)),1));
maxamp=max(max(test))/numel(centvec);

%5. now there may be wrap around in the two phase encode directions
%generate a larger test panel to pick the ROI
largetest=zeros(3*size(test));
sit=size(test);
for xx=0:2;
    for yy=0:2;
    xvec=xx*sit(2)+(1:sit(2));
    yvec=yy*sit(1)+(1:sit(1));
    largetest(yvec,xvec)=test;
    end
end
                                              
                                              
                                              
[digiout, xlim, ylim,ROI]=pickroi(largetest);
xlimvec=mod(xlim(1):xlim(2),sit(2))+1;
ylimvec=mod(ylim(1):ylim(2),sit(1))+1;

ROIvolume=pic(:,ylimvec,xlimvec,:);
roixaxis=xaxis(1:numel(xlimvec));
roiyaxis=yaxis(1:numel(ylimvec));
roizaxis=zaxis;


%6. calculate x y and zaxes for plotting



si=size(ROIvolume);
mag=zeros([si(1:3)]);
t2star=mag;

out=matrix_expfit(ROIvolume,taxis);
mag=out.amplitude;
t2star=out.Tau;



%%
threshold=min(min(test))/numel(centvec)*3;
if isnumeric(timestamp);
filename=[num2str(timestamp) '_movie'];
else
    filename=[timestamp '_movie'];
end


figure('Position',[ 700         500         1250         520]);
for jj=1:si(1); 
    
    MM=squeeze(mag(jj,:,:));
    mask=MM<threshold;
    TT=squeeze(t2star(jj,:,:));
    %MM(mask)=0;
    %TT(mask)=0;
    
    subplot(1,2,1);
    imagesc(roixaxis,roiyaxis,MM,[0 maxamp]); 
    colorbar;
    title(['magnitude, z=' num2str(roizaxis(jj),'%0.3f') 'mm']);
    axis square;
    axis xy;
    xlabel('x (mm)');
    ylabel('y (mm)');
    
    subplot(1,2,2);
    imagesc(roixaxis,roiyaxis,TT*1000,[0 20]); 
    colorbar;
    title('T2* (ms)');
    axis square;
    axis xy;
    xlabel('x (mm)');
    ylabel('y (mm)');
    
    
    %colormap gray;
    pause(0.2); 
    
    %make the gif;
    frame=getframe(gcf);
    im=frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
     % Write to the GIF File
    if jj == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append');
    end
    
end;
