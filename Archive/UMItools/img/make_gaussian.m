function result=make_gaussian(m,sd,pts)
%function result=make_gaussian(m,sd,pts)
%
% x = linspace(0,pts,pts);
% result = exp( - (x-m).^2 / sd^2);
% result = result/sum(result);
%
% (c) 2005 Luis Hernandez-Garcia 
% University of Michigan
% report bugs to:  hernan@umich.edu
%
x = linspace(0,pts-1,pts);
% result = exp( - (x-m).^2 / sd^2);
% result = result/sum(result);
result = normpdf(x,m,sd);
return
