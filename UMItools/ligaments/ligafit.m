function [xfit,yfit,outstruct]=ligafit(varargin)
% routine to fit harmonic surcumference to a roundish object
% geostruct=ligafit(x,y,nharmonics);
% x and y are 1-d vectors containing the x,y coordinates of the measured
% boundary. nharmonics is how many harmonics are being fitted.

%input checking
x=varargin{1};
y=varargin{2};
x = x(:);							%x
y = y(:);	
if numel(x)~=numel(y);
    display('The x and y vectors you handed down do not have the same length!');
end
if nargin>2;
    nharmonics=varargin{3};
else
    nharmonics=4;
    display(['No number of harmonics specified, defaulting to ' num2str(nharmonics) ' harmonics.']);
end

%find the provisional center;
xm=mean(x);     %mean x coordinate
ym=mean(y);     %mean y coordinate

xs=(x-xm);      %shifted and scaled x coordinate
ys=(y-ym);      %shifted and scaled y-coordinate

%find the optimal rotation angle and rotate
for jj=1:36;
    phi(jj)=jj*5/180*pi;
    xy=[xs ys]*[cos(phi(jj)) sin(phi(jj)); -sin(phi(jj)) cos(phi(jj))]; 
    eli(jj,:)=std(xy);
end
[maxx,maxxind]=max(eli(:,1));
[miny,minyind]=min(eli(:,2));
maxxind=maxxind(1);
minyind=minyind(1);

if maxxind==minyind;
    display(['Consistent maxind minind for rotation and scaling']);
else
    display(['maxind=' num2str(maxind) ', minind=' num2str(minind)]);
end

phiindex=maxxind;
xy=[xs ys]*[cos(phi(phiindex)) sin(phi(phiindex)); -sin(phi(phiindex)) cos(phi(phiindex))];
sx=std(xy(:,1));
sy=std(xy(:,2));
xs=xy(:,1)/sx;
ys=xy(:,2)/sy;


%define R and theta with respect to the provisional center
R=sqrt(xs.^2+ys.^2);
theta=angle(xs+1i*ys);
[theta,sortvec]=sort(theta);
R=R(sortvec);

%fit the polar harmonics;
[phiout,Rout]=polarfit(theta,R,nharmonics);
residual=std(R-interp1(phiout,Rout,theta))*sqrt(sx*sy)/sqrt(2);

%calculate x & y
xfit=Rout.*cos(phiout);
yfit=Rout.*sin(phiout);
%scale it
xfit=xfit*sx;
yfit=yfit*sy;
%rotate it back
xy=[xfit yfit]*[cos(-phi(phiindex)) sin(-phi(phiindex)); -sin(-phi(phiindex)) cos(-phi(phiindex))];

%calculate shape properties xcenter, ycenter, area
Rout=sqrt(xy(:,1).^2+xy(:,2).^2);
phiout=angle(xy(:,1)+1i*xy(:,2));
[phiout,phiind]=sort(phiout);
Rout=Rout(phiind);
dphi=diff(phiout);

Rave=(Rout(1:(numel(Rout)-1))+Rout(2:numel(Rout)))/2;
phiave=(phiout(1:numel(phiout)-1)+phiout(2:numel(phiout)))/2;
Area=sum((dphi.*Rave.^2/2));
xcenter=2/3*(Rave.*cos(phiave))' *(dphi.*Rave.^2/2)/Area +xm;
ycenter=2/3*(Rave.*sin(phiave))' *(dphi.*Rave.^2/2)/Area +ym;


outstruct.center_x=xcenter;
outstruct.center_y=ycenter;
outstruct.Area=Area;
outstruct.residual=residual;

%shift the origin back
xfit=xy(:,1)+xm;
yfit=xy(:,2)+ym;


if any(strcmp('plot',varargin));
    figure;
    plot(x,y);
    hold on;
    plot(xfit,yfit,'r');
    title(['Radial deviation ' num2str(residual)])
    grid on;
    
end





function [phiout,Rout,outstructure]=polarfit(varargin)
%[phiout,Rout,outstructure]=polarfit(theta,R,harmonicorder)
%phiout and rout are the fitted circumference
%outstrucure contains the coefficients of the harmonic fit, in the order
%[DC sin(theta) cos(theta) sin(2*theta) cos(2*theta) sin(n*theta) .....]


%[powercoefficients,errors, out]=powersfit(x,y,powers);
%[powercoefficients,errors, out]=powersfit(x,y,ey,powers);
%out.C is the curvaturematrix
%out.P can be used with polyval to generate the fit


if nargin<3;
	display('Input must have at least three arguments: theta, R, n');
	display('where n is the degree of the highest harmonic');
	return
end

x=varargin{1}(:);
y=varargin{2}(:);
%ey=ones(size(y));
harmonics=varargin{3};

%check if there are sufficient data points to constarin the degrees of freedom
if numel(x)<(2*harmonics+1)
	%the data is not constrained sufficiently;
	phiout=[1:360]';
	Rout=eps*ones(360,1);
	outstructure='insufficient points to constrain radial fit';
	return
end

if ~isequal(size(x),size(y))
    error('theta and R vectors must be the same size.')
end

x = x(:);							%x is theta
y = y(:);							%y is R
if nargin<4;
	ey=ones(size(y));			  %this could be upgraded to improve fitting
else
	ey=varargin{4}(:);
	if numel(ey)~=numel(y);
		display('The optional error vector you handed down does not have the same length as the coordinate vectors!');
		return
	end
end

%rescale x for accuracy
x = mod(x,2*pi);

%%%%%%%%%%%
% calculate alpha, beta and x scale matrix
for j = 1:2:(2*harmonics+1);
	%A(:,j)=x.^powers(j)./ey;
	%V(:,j) = x.^powers(j)./ey.^2; 
	%Vm(:,j)= x.^powers(j);
	if j==1;
		A(:,j)=ones(size(x))./ey;
		V(:,j)=ones(size(x))./ey.^2;
		Vm(:,j)=ones(size(x));
	else
		n=(j-1)/2;							%harmonic index				
		%sine
		A(:,j-1)=sin(n*x)./ey;
		A(:,j)   =cos(n*x)./ey;
		
		V(:,j-1)=sin(n*x)./ey.^2;
		V(:,j)   =cos(n*x)./ey.^2;
		
		Vm(:,j-1)=sin(n*x);
		Vm(:,j)   =cos(n*x);
	end
end

beta=y'*V;
alpha=A'*A;

%inversion
p=beta/alpha;

%calculate residuals
r = y - Vm*p';
		
%re-scale the powercoefficients to get back to real units
harmoniccoefficients = p;

%calculates chi squared from the definition
chi2=(r./ey)'*(r./ey);
		
%evaluates the quality of the fit, using chi2 and the number of degrees of freedom
Q=gammainc((length(x)-length(2*harmoniccoefficients+1))/2, chi2/2);

%calculate the fitted curve, using one degree increments
phiout=(-180:180)'/360*2*pi;
for j = 1:2:(2*harmonics+1);
	if j==1;
		Vmo(:,j)=ones(size(phiout));
	else
		n=(j-1)/2;							%harmonic index				
		Vmo(:,j-1)=sin(n*phiout);
		Vmo(:,j)   =cos(n*phiout);
	end
end
Rout=Vmo*p';

outstructure.Q=Q;
outstructure.coeff=p;