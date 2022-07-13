function [fig1, fig2, fig3] =  ov(h,d,x,y,z,roi,cmap)
%function [fig1, fig2, fig3] =  ov(h,d,x,y,z,roi)

stretch = h.zsize/h.xsize;

fig1=subplot (221); 
  image(squeeze(d(:,:,z)));
  colormap(cmap);
  axis ([1 h.ydim 1 h.xdim]);
  axis xy;
fig2=subplot (222); 
  image(squeeze(d(:,y,:))');
  colormap(cmap);
  axis ([1 h.xdim 0 h.xdim/stretch]);
  axis xy;  
fig3=subplot (223); 
  image(squeeze(d(x,:,:))');
  colormap(cmap);
  axis ([1 h.ydim 0 h.xdim/stretch]);
  axis xy;
 return
