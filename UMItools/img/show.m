function out = show(a,windfact)
% usage .. show(a,f);
% displays matrix "a" as a greyscale image in the matlab window
% and "f" is the optional window factors [min,max]

doColorBar=0;

if exist('windfact') == 0, 
  amin = min(min(a));
  amax = max(max(a));
  minmax = [amin,amax];
  a = (a  - amin);
else
  amin = windfact(1);
  amax = windfact(2);
  minmax = [amin,amax];
  a = (a  - amin);
  %a = a .* (a > 0);
end

colormap(gray(256));
imageHandle = image((a)./(amax-amin).*256);
axis('image');
axis('on');
grid on;

if doColorBar
    ncbar_labels=5;
    c1 = colorbar;
    set(c1,'YTick',linspace(0,255,ncbar_labels),...
        'YTickLabel',linspace(amin,amax,ncbar_labels),...
        'FontSize',12);
end

if nargout==0,
  disp(['min/max= ',num2str(minmax(1)),' / ',num2str(minmax(2))]);
end;
