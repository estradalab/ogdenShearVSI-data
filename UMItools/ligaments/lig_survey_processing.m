function h=lig_fig2b_sagittal(varargin)
%coronal and sagittal cuts 

if ischar(varargin{1})
    P=load(varargin{1});
end

if isstruct(varargin{1})
    P=varargin{1};
end

if nargin>1
    sli=varargin{2};
else
    sli=32;
end

D=squeeze(P.QC.data(:,sli,:,:,1));
Dcombined=squeeze(P.data(:,sli,:,1));

%%
hf=0;
if isfield(P,'mask_hr');
    mask=P.mask_hr.*P.mask; %and the masks of the GE3D data and the DESTE data
else
    mask=P.mask;
end
%identify a default center of the sample, using the mask;

D=D.*squeeze(mask(:,sli,:));
Dcombined=Dcombined.*squeeze(mask(:,sli,:));

if isfield(P.HIRES,'orig_data')
	[signal_estimate,noise_estimate]=estimate_noiselevel(P.HIRES.orig_data.image);
else
	[signal_estimate,noise_estimate]=estimate_noiselevel(P.HIRES.magnitude);
end

D=permute(D,[2 1 3]);
Dcombined=permute(Dcombined,[2 1 3]);
xaxis=P.QC.axis1{1};
yaxis=P.QC.axis3{1};

labxo=-34;

figure('Position',[440   144   560   410]);
colormap gray;
gap=20;
%image plus
tsubplot(3,2,1,gap);
imagesc(xaxis,yaxis,abs(D(:,:,1)));
axis image;
lig_angleoverlay(angle(D(:,:,1)));
set(gca,'Xticklabel','','YTick',[-5 0 5],'Ydir','normal');
hyl=ylabel('z (mm)'); set(hyl,'Position', [labxo-2 0 1]); 
hold on;
plot(xaxis, zeros(size(xaxis)),'w');
title('\theta^+_y for \lambda^+_y=+1.0mm');
hcb(1)=cbar;

%line plot of phase through center plus
tsubplot(3,2,3,gap);
plot(xaxis(:),angle(D(32,:,1)),'k.-');
set(gca,'xlim',[xaxis(end) xaxis(1)]);
a=0.26; pos=get(gca,'position'); set(gca,'Position',[pos(1) pos(2)+a pos(3) 0.1]);
set(gca,'xlim',[xaxis(end) xaxis(1)], 'ylim',[-pi pi], 'Ytick', [-pi 0 pi],...
    'Yticklabel',{'-\pi','0','+\pi'});
hyl=ylabel('\theta^+_y (rad)'); set(hyl,'Position', [labxo 0 1]); 
xlabel('y (mm)');
%text(-27,2.2,'z=0');


%image phase with negative polarity
tsubplot(3,2,2,gap);
imagesc(xaxis,yaxis,abs(D(:,:,2)));
axis image
lig_angleoverlay(angle(D(:,:,2)));
set(gca,'Xticklabel','','Ytick',[-5 0 5],'Ydir','normal');
a=0; %y offset
pos=get(gca,'position'); set(gca,'Position',[pos(1) pos(2)+a pos(3) pos(4)]);
hyl=ylabel('z (mm)'); set(hyl,'Position', [labxo-2 0 1]); 
hold on;
plot(xaxis, zeros(size(xaxis)),'w');
title('\theta^-_y for \lambda^-_y=-1.0mm');
hcb(2)=cbar;
%phase minus
tsubplot(3,2,4,gap);
plot(xaxis(:),angle(D(32,:,2)),'k.-');
set(gca,'xlim',[xaxis(end) xaxis(1)], 'ylim',[-pi pi], 'Ytick', [-pi 0 pi],...
    'Yticklabel',{'-\pi','0','+\pi'});
a=0.26; 
pos=get(gca,'position'); set(gca,'Position',[pos(1) pos(2)+a pos(3) 0.1]);
hyl=ylabel('\theta^-_y (rad)'); set(hyl,'Position', [labxo 0 1]); 
xlabel('y (mm)');
%text(-27,2.2,'z=0');


% plot the combined phase image
tsubplot(3,2,3,gap);
imagesc(xaxis,yaxis,abs(Dcombined(:,:,1)));
axis image;
lig_angleoverlay(angle(Dcombined(:,:,1)));
set(gca,'Xticklabel','','YTick',[-5 0 5],'Ydir','normal');
hyl= ylabel('z (mm)'); set(hyl,'Position', [labxo-2 0 1]); 
hold on;
plot(xaxis, zeros(size(xaxis)),'w');
title('\theta_y for \lambda_y=+0.5mm');
%shift y downwards
a=-0.12; pos=get(gca,'position'); set(gca,'Position',[pos(1) pos(2)+a pos(3) pos(4)]);
hcb(3)=cbar;
%phase plus
tsubplot(3,2,5,gap);
plot(xaxis(:),angle(Dcombined(32,:,1)),'k.-');
set(gca,'xlim',[xaxis(end) xaxis(1)]);
a=-0.2; set(gca,'Position',[pos(1) pos(2)+a pos(3) 0.1]);
set(gca,'xlim',[xaxis(end) xaxis(1)], 'ylim',[-pi pi], 'Ytick', [-pi 0 pi],...
    'Yticklabel',{'-\pi','0','+\pi'});
hyl= ylabel('\theta_y (rad)'); set(hyl,'Position', [labxo 0 1]); 
xlabel('y (mm)');
%text(-27,2.2,'z=0');

%derivative
tsubplot(3,2,4,gap);
imagesc(xaxis,yaxis,sum(abs(D),3));
axis image
lig_angleoverlay(angle(Dcombined(:,1:end-1)./(Dcombined(:,2:end)+eps)),'clim',[0 pi/8]);
set(gca,'Xticklabel','','Ytick',[-5 0 5],'Ydir','normal');
a=-0.12; 
pos=get(gca,'position'); set(gca,'Position',[pos(1) pos(2)+a pos(3) pos(4)]);
hyl= ylabel('z (mm)'); set(hyl,'Position', [labxo-2 0 0]); 
hold on;
plot(xaxis, zeros(size(xaxis)),'w');
title('\delta\theta_y/pixel, in the y-direction');
hcb(4)=cbar([0 pi/8]);

%derivative line plot
tsubplot(3,2,6,gap);
ydata=angle(Dcombined(32,1:end-1)./(Dcombined(32,2:end)+eps));
ydata(abs(ydata)>pi/6)=NaN;
plot(xaxis(1:end-1),ydata,'k.-');
set(gca,'xlim',[xaxis(end) xaxis(1)],'Ytick', [0 0.5 1],'ylim',[-0.1 0.6]);
a=-0.2; 
set(gca,'Position',[pos(1) pos(2)+a pos(3) 0.1],'Ytick',[0 pi/16 pi/8],...
    'yticklabel',{'0'; '\pi/16'; '\pi/8'});
hyl=ylabel('\delta\theta_y/pix (rad)'); set(hyl,'Position', [labxo 0.15 -1]); 
xlabel('y (mm)');
grid on;
%text(-27,pi/8,'z=0');


annotation(gcf,'rectangle',...
    [0.0242 0.5585 0.95 0.3975]);

annotation(gcf,'rectangle',...
    [0.0242 0.11 0.45 0.3975]);

annotation(gcf,'rectangle',...
    [0.5222 0.11 0.45 0.3975]);



for jj=1:4;
    set(hcb(jj),'colormap',jet(200));
end


function h=cbar(varargin)
p=get(gca,'position');
h=axes;

if nargin~=0
    clim=varargin{1};
    dummy=(clim(1)+(diff(clim)/200)*(0:200)')*[1 1];
    yticklabel={'0'; '\pi/16';'\pi/8'};
    ytick=[0 pi/16 pi/8];
else
    dummy=(pi*(-1:0.01:1)' * [1 1]);
    yticklabel={'-\pi','0','+\pi'};
    ytick=[-pi 0 pi];
end
imagesc([0 1],dummy(:,1),dummy);
set(gca,'position',[p(1)+p(3)+0.01 p(2)+0.04 0.02 p(4)-0.08]);
axis on;

set(gca,'xticklabel','','colormap',jet(256),'YTick',ytick,...
    'YaxisLocation','right','Yticklabel',yticklabel,'YDir','normal');


