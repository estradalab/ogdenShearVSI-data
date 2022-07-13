function [t, beta_est, v, RSS]= my_glm(x,y,c)
%function [t, [beta_est, v, RSS]] = my_glm(x,y,c)
% t is the resulting T score of the comparison
%
% x is the design matrix
% c is a vector of contrasts
% y is the vector of data
%
% (c) 2005 Luis Hernandez-Garcia 
% University of Michigan
% report bugs to:  hernan@umich.edu
%
% unlike its predecessor, this version will not add an intercept to
% the design matrix automatically.  Add your own damn intercept.

    y = reshape(y,length(y),1);
	x = reshape(x,length(y),length(c));
	c = reshape(c,length(c),1);

	%x = [ones(size(x,1) , 1)  x];
	%c = [0; c];
	
	p = size(x,2);
	n = size(x,1);

	xtx_inv = pinv(x'*x);
	beta_est = xtx_inv*x'*y;

	RSS = y - x*beta_est;
	RSS = sum(RSS.^2);

	var_est = RSS/(n-p);
	cov_beta = var_est * xtx_inv;

	% ---- %

	v = var_est*(c' * xtx_inv * c);
	t = (beta_est' * c) ./ sqrt(v);

    Nargs = nargout;
    if Nargs <= 2
        save glm.mat beta_est var_est RSS
    end
	% For testing only.  Change this when you are done !!!
	%t = beta_est' * c;
	%t = v;
	%t = RSS;

return
