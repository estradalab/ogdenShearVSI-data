function [func_hat, var_hat] = hrf_deconv(data, onsets, length)
% function [func_hat, var_hat] = hrf_deconv(data, onsets, length)
%
% deconvolves the hrf from a time series, given:
%
% data - (of course)this needs to be a column
% onsets -  the onsets of the events 
% length - of the HRF to estimate desired
% 
% the function uses a linear regression approach to 
% deconvolve the function.
%
% note - it also adds another regressor to the matrix in order to 
% estimate the DC offset.  This is not reported in the final estimate
%
% it returns:
% fun - the resulting HRF function
% v - the variance estimate for the function
%
%data = data + 100;


x = zeros(size(data,1)+length, length);

for count=0:length-1
	x(onsets+count, count+1) = 1;
end
x = x(1:size(data,1) , 1:end);
% add a regressor to estimate the baseline 
x = [x ones(size(data,1),1)];

[func_hat, var_hat]= linreg(x,data);

baseline = func_hat(end);
func_hat = func_hat(1:end-1);
% add the baseline back in.
func_hat = func_hat + baseline;
return

function [beta_est, var_est] = linreg(x,y)
%function [beta_est, var_est] = linreg(x,y)
%
% gives the least squares estimate for beta in the 
% problem :
%           y = x*beta + error
%
%
% x is the design matrix
% y is the vector of data
% beta_est is the parameter estimates
% var_est is the estimate of the variance on those betas
%
	
	p = size(x,2);
	n = size(x,1);

	%beta_est = (pinv(x'*x))*x'*y;
	beta_est = pinv(x)*y;

	RSS = y - x*beta_est;
	RSS = sum(RSS.^2);

	var_est = RSS/(n-p);
	%cov_beta = var_est * xtx_inv;

	% ---- %

	%v = var_est*(c' * xtx_inv * c);
	%t = (beta_est' * c) ./ sqrt(v);

	%save glm.mat beta_est var_est

	% For testing only.  Change this when you are done !!!
	%t = beta_est' * c;
	%t = v;
	%t = RSS;

return
