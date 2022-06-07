function update_monitorplot
%callback routine to update the monitoring plot produced by
%matrix_invrecfit, with option 'monitorplot'
XYdummy=get(gca,'CurrentPoint');
callingaxes=gca;
callingfigure=gcf;
xi=round(XYdummy(1,1));
yi=round(XYdummy(1,2));

global chimon;
global devmon;
global taumon;
persistent crosshandle1 crosshandle2

if ishandle(crosshandle1);
    delete(crosshandle1);
end

if ishandle(crosshandle2);
    delete(crosshandle2);
end


si=size(chimon);




dummy=findobj('Tag','taumonitor');
htau=dummy(1);
axes(htau);
plot(squeeze(taumon(yi,xi,:)),'-o');
set(htau,'Tag','taumonitor','YLim',[0 3],'Xlim',[1 si(3)]);
title('\tau');

dummy=findobj('Tag','chimonitor');
hchi=dummy(1);
axes(hchi);
plot(squeeze(chimon(yi,xi,:)),'-o');
set(hchi,'Tag','chimonitor','YLim',[0 2],'Xlim',[1 si(3)]);
title('\chi^2');

dummy=findobj('Tag','devmonitor');
hdev=dummy(1);
axes(hdev);
plot(squeeze(devmon(yi,xi,:)),'-o');
set(hdev,'Tag','devmonitor','YLim',[-0.1 1],'Xlim',[1 si(3)])
title('deviation');

figure(callingfigure);
maghandle=findobj('Tag', 'Magimage');
axes(maghandle);
title(['(x,y)=(' num2str(xi) ',' num2str(yi) ')']);
crosshandle1=plot(xi,yi,'+k');
xlim=get(maghandle,'XLim');
ylim=get(maghandle,'YLim');
tauhandle=findobj('Tag', 'Tauimage');
set(tauhandle,'Xlim',xlim,'ylim',ylim);
axes(tauhandle);
crosshandle2=plot(xi,yi,'+k');








