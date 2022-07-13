function result = mydetrend( data, saveCoeffs)
%function result = mydetrend( data [,saveCoeffs])
% this function fits and subtracts a third order polynomial
% from the data.
%
% (c) 2005 Luis Hernandez-Garcia 
% University of Michigan
% report bugs to:  hernan@umich.edu
%
	num = max(size(data));
	%t=[1:max(size(data))]';
	t = linspace(0,num, num);
    t = reshape(t,size(data));
	[coeffs, error] = polyfit(t,data,4);
	result = data-polyval(coeffs,t);
    
    if nargin<2
        save coeffs
    end
return
