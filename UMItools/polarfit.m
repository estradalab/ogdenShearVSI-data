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
phiout=(1:360)'/360*2*pi;
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