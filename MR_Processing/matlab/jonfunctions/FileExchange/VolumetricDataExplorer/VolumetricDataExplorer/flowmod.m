function [x,y,z,v] = flowmod(x,y,z,t)
%FLOW A simple function of 3 variables.
%   FLOW, a function of three variables, generates fluid flow data 
%   that is useful for demonstrating SLICE, INTERP3, and other 
%   functions that visualize scalar volume data.
%   
%   There are several variants of the calling sequence:
%       V = FLOW; produces a 25-by-50-by-25 array.
%       V = FLOW(N) produces a N-by-2N-by-N array.
%       V = FLOW(X,Y,Z) evaluates the speed profile at the points (X,Y,Z).
%       [X,Y,Z,V] = FLOW(...) returns the coordinates as well.

%   Reference: Fluid Mechanics, L. D. Landau and E. M. Lifshitz.
%   Copyright 1984-2013 The MathWorks, Inc.
%   $Revision: 1.8.4.6 $ $Date: 2011/11/17 21:06:44 $
v(length(y),length(x),length(z),length(t))=0;
if nargin==0,
  [x,y,z] = meshgrid(.1:.2:10,-3:.25:3,-3:.25:3);
elseif nargin==1, % flow(N)
  n = x;
  [x,y,z] = meshgrid(linspace(.1,10,2*n),linspace(-3,3,n),linspace(-3,3,n));
elseif nargin==2, 
  error('MATLAB:flow:InputCountIncorrect','%s',...
      getString(message('MATLAB:demos:shared:InputCountIncorrect')));
elseif nargin==4
    [x,y,z] = meshgrid(x,y,z);
end

% Convert to spherical coordinates (with x as the axis).
% A = 2; nu = 1;

[~,phi,r] = cart2sph(y,z,x);
for ii = 1:length(t)
    A = 2*sin(2*t(ii))+4;
    nu = cos(2*(t(ii)))+1.5;
    vr = 2*nu/r.*((A^2-1)/(A-cos(phi)).^2 - 1);
    vphi = -2*sin(t(ii))*nu*sin(phi)./(A-cos(phi))./r;
    vth = zeros(size(r));

    [vx,vy,vz] = sph2cart(vth,vphi,vr);
    v(:,:,:,ii) = -log(sqrt(vx.^2 + vy.^2 + vz.^2));
end
if nargout < 2, x = v; end
