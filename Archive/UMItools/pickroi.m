function [digiout, xlim, ylim,ROI]=pickroi(varargin)
inmatrix=varargin{1};
inmatrix=abs(inmatrix);
if ndims(inmatrix)==3;
	inmatrix=mean(inmatrix,3);
end

%default: take the center of the roi, and the smaller range x and y, to
%create a square ROI;
centerflag=false;

fh=figure('Toolbar','figure','Position',[350    60   700   560]);
h=imagesc(inmatrix);
axis image;
set(gca,'YDir','normal');
if any(strcmp(varargin, 'gray'));
    colormap gray;
end
title('Zoom in on ROI and click the ROI button','Fontsize', 14);
zoom on
uih=uicontrol('Style','pushbutton','String', 'ROI', 'callback','uiresume');
uiwait(gcf);

uicontrol(uih,'String', 'done');
title('ROI accepted, processing now','Fontsize', 14);


xlim=get(gca,'xlim'); 
ylim=get(gca,'ylim'); 

if centerflag;
    xcenter=    mean(xlim);
    ycenter=    mean(ylim);
    xrange=     diff(xlim)/2;
    yrange=     diff(ylim)/2;
    range=round(min(xrange,yrange));
    
    xlim=[ceil(xcenter-range) floor(xcenter+range)];
    ylim=[ceil(ycenter-range) floor(ycenter+range)];
else
    xlim=[ceil(xlim(1)) floor(xlim(2))];
    ylim=[ceil(ylim(1)) floor(ylim(2))];
end
    


si=size(inmatrix);
ivec=round(max(1,ylim(1))):round(min(ylim(2),si(1)));
jvec=round(max(1,xlim(1))):round(min(xlim(2),si(2)));

digiout=false(si);
digiout(ivec,jvec)=true;

ROI={'ROI limits'; xlim; ylim;jvec;ivec};

close(fh);