% function power_img = spm_power(alpha, df, error)
% function power_img = spm_power(alpha, df, error)
%
% computes the power of a t-test accross
% several images of every voxel.
%
% Luis Hernandez
% University of Michigan
% Last Edit 1-27-2000
%
% The null hypothesis is that each images has a value equal
% to the mean value of all the images in the group.
% The alternative hypothesis is that they are not (double sided test)
% a p-value less than alpha means that they are different. 
%
% Arguments and defaults:
% 	alpha: 	significance level (0.05)
% 	df :   	degrees of freedom (number of images * average number of Ressels)
% 	error: 	Numbber of standard deviations of error to assume, when
%			calculating the power (1)
% 	files.txt: a text file containing the filenames (with path) of
%			of the images to test
%
% Output:
% 	power.img 
%	p-value.img
%


hold off

% generate dummy data
 data = [100:5:130];

% Read in the data
% pFile = fopen('./files.txt','r');
% if pFile == -1
%   msgbox('I need files.txt, which must contain the names of the files to analyze');
%   return
% end

%old_path=pwd;
%done=0;
%
%while ~done  
%   [filename newpath] = uigetfile('*.img','SELECT IMAGES FOR TEST, hit CANCEL when finished',200,200);
%   if newpath==0,
%      done=1;
%   else%
%	   cd (newpath)
%	   hdrname = [filename([1:size(filename,2)-4]) '.hdr']
%	   h = read_hdr(hdrname);   
%	   img_data = read_img_data(h,filename);
%	   if(~exist('data'))
%	      data = img_data';
%	   else 
%	      data = [data img_data'];
%      end
%      clear img_data;
%   end   
%end



% compute NUmber Average Number of Ressels
avgRessels = size(data,2)*100;
% avgRessels = 150

%if nargin == 0,
   alpha = 0.05/2;  % Significance level
   df = avgRessels-1;
	error = 1;
%end

% compute Z-scores for all images
x_mean  = mean(data,2);
x_var = (var(data'))';
x_sigma = sqrt(x_var);
%clear x_var;

% change the size of the matrices to avoid for loops
x_mean = x_mean * ones(1,size(data,2));
x_sigma = x_sigma * ones(1,size(data,2));

disp('computing t_score ...')
t_score = (data - x_mean) ./ (x_sigma/sqrt(avgRessels) )
disp('done')

disp('computing p_value ...')
p_value = (1 - tcdf(t_score,df))
mean_p_value =  mean(p_value,2)
%clear p_value;
disp('done')

% make a plot of the t-ditribution
t=[-100:100];
f = tpdf(t,df);
plot( t , f,'b')
hold on
% pause

% compute power assuming that we are wrong and the truth is one 
% some standard deviations higher. The t_score in this ideal distribution
% is then given by
% t_ideal = t_score - error * sqrt(avgRessels)

% make a plot of the ideal ditribution
t=[-100:100];
f = tpdf(t - error * sqrt(avgRessels),df);
plot(t,f,'g')
%
	
	
% find critical t-score for this significance level(alpha)
% and plot it
t = 0;
p=alpha*2;  %(some initial value)
while (p >= alpha),
  		p = (1-tcdf(t,df))*2;
  		t = t + 0.01;
end
t_alpha = t;

plot(t_alpha, tpdf(t_alpha, df), 'r*')
plot(t_alpha, tpdf(t_alpha - (error) * sqrt(avgRessels), df), 'r*')


disp('computing beta ...')
beta = tcdf(t_alpha - error * sqrt(avgRessels), df)
power = 1 - beta
disp('done')

mean_power =  mean(power,2);
   

%cd(old_path);
%
%write_hdr('power.hdr',h);
%write_img_data('power.img',mean_power,h);
%write_hdr('p-value.hdr',h);
%write_img_data('p-value.img',mean_p_value,h);


%return
