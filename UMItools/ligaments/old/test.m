% testing the shape fit

%define a test shape, deformed egg with bumps

aspectratio=5;      %aspectratio
phi=-pi/7;          %rotation of the shape
nharmonics=8;       %degrees of freedom in Fourier fit: 2*nharmonics+1

noiselevel=0.1;

xoffset=10;         
yoffset=5;


theta=[0:1000]/1000*2*pi;
%define R(theta) for the shape, principal axes along the xy coordinates
R=0.1*cos(theta.^3)+2*sin(theta)+ 0.3*cos(3*theta).^2 +5;
x=aspectratio*R.*cos(theta)+noiselevel*randn(size(theta));
y=R.*sin(theta)+noiselevel*randn(size(theta));
x=x(:);
y=y(:);

%rotate the test shape for arbitrary principal axis direction
xy=[x y]*[cos(phi) sin(phi); -sin(phi) cos(phi)]; 
xr=xy(:,1)+xoffset;
yr=xy(:,2)+yoffset;

%call the fitting function, with optional 'plot' to overlay 
%data and fit, geostruct contains fit and data about the shape
[xfit,yfit,outstruct]=ligafit(xr,yr,nharmonics,'plot');
display(outstruct)

