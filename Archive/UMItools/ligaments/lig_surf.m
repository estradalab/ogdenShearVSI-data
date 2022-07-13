function out=lig_surf(varargin)
% playing around with surface and volume and 3d plots, to see if a
% reasonable way can be found to plot strains in 3d
% options '11' '12' '22' '13' '23' '33' 'angle1' 'angle2' 'angle3' 'sub'

P=varargin{1};


xlim=1*[-12 12];
ylim=1*[-25 30];
zlim=1*[-12 12];
clim=0.05*[0 1];

title_str='E11';


[X,Y,Z]=meshgrid(P.axis2,P.axis1,P.axis3);

mask=P.mask_hr;
nanmask=mask;
nanmask(nanmask(:)==0)=NaN;

%Lagrange=P.strains(:,:,:,2,2).*nanmask;
%Lagrange=angle(P.data(:,:,:,1)).*nanmask;
Lagrange=P.Lagrange(:,:,:,1,1).*nanmask;

wcal=P.wiggle(1)/2;

load spectral

%default color appearance
cmap=spectral(110:255,:);
shade='interp';

if nargin>1;
    switch varargin{2};
        case '11'
            Lagrange=P.Lagrange(:,:,:,1,1).*nanmask;
            clim=0.08*[0 1]*wcal;
            title_str='\epsilon_{yy}';
        case '12'
            Lagrange=P.Lagrange(:,:,:,1,2).*nanmask;
            clim=0.02*[-1 1]*wcal;
            title_str='\epsilon_{yx}';
        case '22'
             Lagrange=P.Lagrange(:,:,:,2,2).*nanmask;
            clim=0.04*[-1 0]*wcal;
            title_str='\epsilon_{xx}';
            cmap=spectral(1:128,:);
        case '33';
             Lagrange=P.Lagrange(:,:,:,3,3).*nanmask;
            clim=0.04*[-1 0]*wcal;
            title_str='\epsilon_{zz}';
            cmap=spectral(1:128,:);
        case '13'
             Lagrange=P.Lagrange(:,:,:,1,3).*nanmask;
            clim=0.03*[0 1]*wcal;
            title_str='\epsilon_{yz}';
        case '23'
             Lagrange=P.Lagrange(:,:,:,2,3).*nanmask;
            clim=0.01*[-1 1]*wcal;
            title_str='\epsilon_{xz}';
        case 'angle1'
            Lagrange=angle(P.data(:,:,:,1)).*nanmask;
            clim=pi*[-1 1];
            title_str='\theta_y';
            cmap=hot(255);
            shade='interp';
            ytickvalue=(-4:1:4)*pi/4;
            yticklabellist={'-\pi'; '';'-\pi/2';'';'0';'';'+\pi/2';'';'+\pi'};
        case 'angle2'
            Lagrange=angle(P.data(:,:,:,2)).*nanmask;
            clim=pi*[-1 1];
            title_str='\theta_x';
            cmap=hot(255);
            shade='interp';
            ytickvalue=(-4:1:4)*pi/4;
            yticklabellist={'-\pi'; '';'-\pi/2';'';'0';'';'+\pi/2';'';'+\pi'};
        case 'angle3'
            Lagrange=angle(P.data(:,:,:,3)).*nanmask;
            clim=pi*[-1 1];
            title_str='\theta_z';
            cmap=hot(255);
            shade='interp';
            ytickvalue=(-4:1:4)*pi/4;
            yticklabellist={'-\pi'; '';'-\pi/2';'';'0';'';'+\pi/2';'';'+\pi'};
            
    end
end

%title_str=[title_str ', \deltay = '  num2str(P.wiggle(1)) 'mm'];

% make a new figure; unless 'sub' is specified
if ~any(strcmp(varargin,'sub'))
    figure('Position',[1780 850 560 420],'Name', [ '\deltay = '  num2str(P.wiggle(1)) 'mm']);
else
    set(gcf,'Name', [ '\deltay = '  num2str(P.wiggle(1)) 'mm']);
end



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
lightangle(40,40);
lighting flat;
hiso.SpecularColorReflectance = 0;
hiso.SpecularExponent = 50;



%now plot the central projections of planes through the origin onto the
%respective planes
hold on;
%set(gca,'colormap',jet(256));
set(gca,'colormap',cmap);
caxis(clim);
hx = slice(X,Y,Z,Lagrange,0,[],[]);
hx.FaceColor = shade;
hx.EdgeColor = 'none';
%move it to the edge
hxc=hx;
hx.XData=xlim(1)*ones(size(hx.XData));



hy = slice(X,Y,Z,Lagrange,[],0,[]);
hy.FaceColor = shade;
hy.EdgeColor = 'none';
hy.YData=ylim(2)*ones(size(hy.YData));


hz = slice(X,Y,Z,Lagrange,[],[],0);
hz.FaceColor = shade;
hz.EdgeColor = 'none';
hz.ZData=zlim(1)*ones(size(hz.ZData));


%now plot the lines around the gray outline
perimeterflag=false;
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
    
else
%     hy = slice(X,Y,Z,Lagrange,[],0,[]);
%     hy.FaceColor = 'none';
%     hy.EdgeColor = 'k';
end

vertexflag=true;
if vertexflag          %draw an empty polygon to markt the planes we are looking at
    graylevl=0.7;
    v=[0 -20 -5.5; 0 23 -5.5; 0 23 6; 0 -20 6]; 
    f= [1 2 3 4];
    patch('Faces',f,'Vertices',v,'Edgecolor','k','Facecolor','none','Linewidth',1,'Linestyle','-.');
    patch('Faces',f,'Vertices',v+repmat([xlim(1) 0 0],4,1) ,'Edgecolor',graylevl*[1 1 1],'Facecolor','none','Linewidth',1,'Linestyle','-');
    
    v=[4 -20 0; 4 23 0; -4 23 0; -4 -20 0]; 
    f= [1 2 3 4];
    patch('Faces',f,'Vertices',v,'Edgecolor','k','Facecolor','none','Linewidth',1,'Linestyle','-.');
    patch('Faces',f,'Vertices',v+repmat([0 0 zlim(1)],4,1) ,'Edgecolor',graylevl*[1 1 1],'Facecolor','none','Linewidth',1,'Linestyle','-');
    
    v=[4 0 -5.5; 4 0 6; -4 0 6; -4 0 -5.5]; 
    f= [1 2 3 4];
    patch('Faces',f,'Vertices',v,'Edgecolor','k','Facecolor','none','Linewidth',1,'Linestyle','-.');
    patch('Faces',f,'Vertices',v+repmat([0 ylim(2) 0],4,1) ,'Edgecolor',graylevl*[1 1 1],'Facecolor','none','Linewidth',1,'Linestyle','-');
else
    text(-12,-20,4.5,'x=0','rotation',22);
    text(-1,30,5.5,'y=0');
    text(5,-19,-12,'z=0','rotation',30);
end



set(gca, 'xlim',xlim,'ylim',ylim,'zlim',zlim);
hxl=xlabel('x (mm)');
pos=get(hxl, 'Position');
set(hxl,'Position',pos+[-2 5 0]);
hyl=ylabel('y (mm)');
pos=get(hyl, 'Position');
set(hyl,'Position',pos+[-5 5 0]);
zlabel('z (mm)');
caxis(clim);
%title(title_str);

%
grid on;

ti=text(xlim(1),ylim(1)+5,zlim(2)+14,title_str);
ti.FontSize=16;
% make and resize colorbar
colorbar_flag=true;
if colorbar_flag;
    posax=get(gca,'Position');
    cboff=[posax(1)+0.32 posax(2)+0.08 0.01 0.5*posax(4)];
    hcb=colorbar('Location','manual', 'Position',cboff);
    pos=get(hcb,'Position');
    %vertshrink=0.1;
    %set(hcb,'Position', [pos(1) pos(2)+vertshrink 0.01 pos(4)-2*vertshrink]);

    if exist('ytickvalue','var');
        set(hcb,'Ytick',ytickvalue,'YTicklabel',yticklabellist);
    end
end
%set(hcb,'Position', [pos(1) pos(2)+vertshrink pos(3) pos(4)-vertshrink]);





