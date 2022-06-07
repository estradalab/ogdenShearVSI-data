function [fig1, fig2, fig3] =  oq(h,d, dx,dy,dz,x,y,z,roi)
%function [fig1, fig2, fig3] =  oq(h,dx,dy,dz,,x,y,z,roi)
%
% (c) 2005 Luis Hernandez-Garcia
% University of Michigan
% report bugs to:  hernan@umich.edu
%


% if isempty(h)


fig1=subplot (223);axis xy ,axis square
image(squeeze(d(:,:,z))')
hold on; plot(y,x,'go');
quiver(squeeze(dx(:,:,z))', squeeze(dy(:,:,z))')
hold off;
set(gca,'Position',[ 0.05    0.01    0.45    0.45]);
axis xy ,axis square

fig2=subplot (222);
image(squeeze(d(:,y,:))'),
hold on; plot(x,z,'go');
quiver(squeeze(dx(:,y,:))', squeeze(dz(:,y,:))')
set(gca,'Position',[ 0.5    0.5    0.45    0.45]);
hold off;
axis xy ,axis square

fig3=subplot (221);
image(squeeze(d(x,:,:))')
hold on; plot(y,z,'go');
quiver(squeeze(dy(x,:,:))', squeeze(dz(x,:,:))')
set(gca,'Position',[ 0.05    0.5    0.45    0.45]);
hold off;
axis xy ,axis square

if roi>0
    value=mean(mean(d(x-roi:x+roi, y-roi:y+roi,z-roi:z+roi)));
else
    value=d(x,y,z);
end
fprintf('\n(x,y,z)=  (%d %d %d) , val= %6.2f  \n', x, y, z, value);
% else
%     stretch = h.zsize/h.xsize;
%     %colordef black
%     subplot(2,2,4),
%     set(gca,'Position',[ 0.5    0.05    0.4   0.3])
%
%     fig1=subplot (223);
%     image(squeeze(d(:,:,z))), axis ([1 h.ydim 1 h.xdim]) ,axis xy ,axis square
%     hold on; plot(y,x,'go');hold off;set(gca, 'XTick', []); set(gca, 'YTick', []);
%     set(gca,'Position',[ 0.05    0.01    0.45    0.45]);
%
%     fig2=subplot (222);
%     image(squeeze(d(:,y,:))'), axis ([1 h.xdim 0 h.xdim/stretch]),axis xy , axis square
%     hold on; plot(x,z,'go');hold off;set(gca, 'XTick', []); set(gca, 'YTick', []);
%     set(gca,'Position',[ 0.5    0.5    0.45    0.45]);
%
%     fig3=subplot (221);
%     image(squeeze(d(x,:,:))'), axis ([1 h.ydim 0 h.xdim/stretch]), axis xy, axis square
%     hold on; plot(y,z,'go');hold off;set(gca, 'XTick', []); set(gca, 'YTick', []);
%     set(gca,'Position',[ 0.05    0.5    0.45    0.45]);
% end

%

%str = sprintf('\n(x,y,z)=  (%d %d %d) , val= %6.2f  \n', x, y, z, tmp);
%subplot(221), title(str)

return
