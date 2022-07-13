function h=textn(varargin)
% textn(varargin)
% takes the same arguments as text, but the x and y coordinates are
% normalized, so you can place the text in the smae corner

xlim=get(gca,'xlim');
ylim=get(gca,'ylim');

varargin{1}=xlim(1)+diff(xlim)*varargin{1}; %replace normalized x with coordinate x
varargin{2}=ylim(1)+diff(ylim)*varargin{2}; %replace normalized y with coordinate y

h=text(varargin{:});
