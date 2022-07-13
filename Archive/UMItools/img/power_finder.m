function Result = Spower_finder(x, alpha, bf, x_effect, Spower_threshold)
% function [Spower N_min]= Spower_finder(x, alpha, bf, effect, Spower_threshold)
%
% computes the Number of measurements
% to achieve a given Spower for a t-test accross
% several images of every voxel.
%
% Luis Hernandez
% University of Michigan
% Last Edit 12 - 3 -2009
%
% The null hypothesis is that each images has a value equal
% to the mean value of all the images in the group.
% The alternative hypothesis is that they are not (double sided test)
% a p-value less than alpha means that they are different.
%
% This function computes the probability of being right if we reject the
% null hypothesis of being equal, based on a given effect size.
%
% Arguments and defaults:
%   x:       2D data array. rows = pixels, cols = subjects
%   alpha: 	significance level (0.05)
%   bf :   	Bonferroni correction factor.
%   effect:   Effect size (%) to be assumed, when calculating the Spower (5%)
%   Spower_threshold:       Desired Spower level
%
%

if nargin == 0
    x=[ 12 12 11 14 15 16 7 3 4 ]';
    bf = 1 ;
    alpha = 0.05;
    x_effect = mean(x);
    Spower_threshold = 0.8;
    
end


% remove zeros and NaNs.  ignore them...
x = x(find(~isnan(x)));
x = x( find(abs(x) > 0) );

if isempty(x)
    result = [nan nan nan];
    return
end


% if the effect size is entered as zero, we choose 
% the effect size automatically to be the mean of the data
if x_effect==0
    x_effect=mean(x);
end

sigma_x = var(x);
N_subjects  = size(x,1);
df = N_subjects;
alpha = alpha/(bf*2);  % 2-tailed Significance level after Bonferroni correction
Spower = 0;

% find critical score for this significance level(alpha)
% the variable 'probabilities' is a table of values corresponding to
% p-values, and the number of standard deviations between the
% critical point (x_alpha) and the mean.

tcrit =  spm_invTcdf(1-alpha, df)
xcrit = tcrit * sigma_x;
q = spm_Ncdf(xcrit, x_effect, var(x));
Spower = 1-q;
current_Spower = Spower;

% Loop for an increasing number of subjects (up to three times as many subjects as we have...)
% set the upper bound to 50 data points required.
N = N_subjects;
while ( ( N < 50)   & (Spower < Spower_threshold)  )
    N = N+1;
    df = N;
    
    sigma_x = var(x)/(sqrt(N / N_subjects));
    
    % what is the critical T score
    tcrit =  spm_invTcdf(1-alpha, df);
    % translate that into a parameter
    xcrit = tcrit * sigma_x;
    
    q = spm_Ncdf(xcrit, x_effect, sigma_x);
    Spower = 1-q;
    fprintf('\n Effect = %2.2f xcrit = %2.2f   sigma_x = %2.2f   N=  %d Power= %2.2f \%', ...
        x_effect, xcrit, sigma_x, N, Spower*100);

end


Result=[current_Spower N Spower];



return

