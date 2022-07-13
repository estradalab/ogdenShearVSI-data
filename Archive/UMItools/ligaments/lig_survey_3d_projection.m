function out=lig_surf(varargin)
% playing around with surface and volume and 3d plots, to see if a
% reasonable way can be found to plot strains in 3d

P=varargin{1};



xlim=[-12 12];
ylim=[-25 30];
zlim=[-12 12];
clim=0.05*[-1 1];
title_str='E11';


[X,Y,Z]=meshgrid(P.axis2,P.axis1,P.axis3);

mask=P.mask.*P.mask_hr;
nanmask=mask;
nanmask(nanmask(:)==0)=NaN;

%Lagrange=P.strains(:,:,:,2,2).*nanmask;
%Lagrange=angle(P.data(:,:,:,1)).*nanmask;
Lagrange=P.Lagrange(:,:,:,1,1).*nanmask;

wcal=P.wiggle(1)/2;

if nargin>1;
    switch varargin{2};
        case '11'
            Lagrange=P.Lagrange(:,:,:,1,1).*nanmask;
            clim=0.08*[0 1]*wcal;
            title_str='E_{yy}';
        case '12'
            Lagrange=P.Lagrange(:,:,:,1,2).*nanmask;
            clim=0.02*[-1 1]*wcal;
            title_str='E_{yx}';
        case '22'
             Lagrange=P.Lagrange(:,:,:,2,2).*nanmask;
            clim=0.04*[-1 0]*wcal;
            title_str='E_{xx}';
        case '33';
             Lagrange=P.Lagrange(:,:,:,3,3).*nanmask;
            clim=0.04*[-1 0]*wcal;
            title_str='E_{zz}';
        case '13'
             Lagrange=P.Lagrange(:,:,:,1,3).*nanmask;
            clim=0.03*[0 1]*wcal;
            title_str='E_{yz}';
        case '23'
             Lagrange=P.Lagrange(:,:,:,2,3).*nanmask;
            clim=0.01*[-1 1]*wcal;
            title_str='E_{xz}';
    end
end
            
title_str=[title_str ', stretch= '  num2str(P.wiggle(1)) 'mm'];

figure('Position',[1780 850 560 420]);
%Ds=smooth3(mask);
Ds=mask;

hiso = patch(isosurface(X,Y,Z,Ds,0.5),...
   'FaceColor',[0.5,0.5,0.5],...
   'FaceAlpha',0.4,...
   'EdgeColor','none');
   
isonormals(Ds,hiso);

%Define the view
view(43,32) 
%axis tight 
%daspect([1,1,.4])
axis image;

%Add Lighting
lightangle(45,30);
lighting flat;
hiso.SpecularColorReflectance = 0;
hiso.SpecularExponent = 50;


%now plot the central projections of planes through the origin onto the
%respective planes
hold on;
set(gca,'colormap',jet(256));
caxis([-0.1 0.1]);
hx = slice(X,Y,Z,Lagrange,0,[],[]);
hx.FaceColor = 'interp';
hx.EdgeColor = 'none';
%move it to the edge
hxc=hx;
hx.XData=xlim(1)*ones(size(hx.XData));



hy = slice(X,Y,Z,Lagrange,[],0,[]);
hy.FaceColor = 'interp';
hy.EdgeColor = 'none';
hy.YData=ylim(2)*ones(size(hy.YData));


hz = slice(X,Y,Z,Lagrange,[],[],0);
hz.FaceColor = 'interp';
hz.EdgeColor = 'none';
hz.ZData=zlim(1)*ones(size(hz.ZData));


%now plot the lines around the gray outline
perimeterflag=true;
if perimeterflag;
    dummy=zeros(size(Ds));
    bw=dummy;
    bw(:,:,32)=bwperim(Ds(:,:,32),8);
    hold on;
    peri = patch(isosurface(X,Y,Z,bw,0.5),...
        'FaceColor',[0.5,0.5,0.5]/3,...
        'FaceAlpha',0.2,...
        'EdgeColor','none');
    
    dummy=zeros(size(Ds));
    bw=dummy;
    bw(:,32,:)=bwperim(squeeze(Ds(:,32,:)),8);
    hold on;
    peri = patch(isosurface(X,Y,Z,bw,0.5),...
        'FaceColor',[0.5,0.5,0.5]/3,...
        'FaceAlpha',0.2,...
        'EdgeColor','none');
    
    dummy=zeros(size(Ds));
    bw=dummy;
    bw(96,:,:)=bwperim(squeeze(Ds(96,:,:)),8);
    hold on;
    peri = patch(isosurface(X,Y,Z,bw,0.5),...
        'FaceColor',[0.5,0.5,0.5]/3,...
        'FaceAlpha',0.2,...
        'EdgeColor','none');
end








set(gca, 'xlim',xlim,'ylim',ylim,'zlim',zlim);
xlabel('x (mm)');
ylabel('y (mm)')
zlabel('z (mm)');
caxis(clim);
title(title_str);

%
grid on;
colorbar;



