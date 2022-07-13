function hrf = make_hrf(delta,tau, npoints)
% function hrf = make_hrf(delta,tau, npoints)
%
% Generates a gamma variate function for the hemodynamic response function.
% points is the number of points to be included in the function.
% (ie - the length of the response to be calculated)
% Normalized so that it peaks at 1
%
% delta 	- delay in the function
% tau		- decay/uptake constant
% npoints	- number of points to generate
%
% NOTE:  must be consistent in units!
%
% (c) 2005 Luis Hernandez-Garcia 
% University of Michigan
% report bugs to:  hernan@umich.edu
%
	   
	t = linspace(0,npoints-1,npoints);
	hrf = ((t-delta)./tau).^2  .*  exp(-(t-delta)./tau)  .* ((t-delta) > 0);
   	hrf = hrf/max(hrf);
   
return
