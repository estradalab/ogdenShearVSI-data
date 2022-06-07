function N = power_ranger(alpha, bf, effect_size, power_threshold)
% function N = power_ranger(alpha, bf, effect_size, power_threshold)
%
% computes the power of a t-test accross
% several images of every voxel.
%
% Luis Hernandez
% University of Michigan
% Last Edit 1 - 25 -2005
%
% The null hypothesis is that each images has a value equal
% to the mean value of all the images in the group.
% The alternative hypothesis is that they are not (double sided test)
% a p-value less than alpha means that they are different. 
%
% USAGE:  the program will put up a file selection box.  
%     select the images of the activation maps you have and 
%     hit the CANCEL button once you have selected all of them
%     The program will compute the power of your test with your
%     current number of samples and the number of samples needed
%     to achive the desired power.  This is done pixel by pixel,
%     so what you get is two images.
%
% Arguments and defaults:
%   alpha: 	significance level (default ... 0.05)
%   bf :   	Bonferroni correction factor (default ... 1). 
%   effect_size: 	what we guess that the Effect's value is.
%   power_threshold:  the desired power level (from 0 to 1) for which the program
%            will calculate the minimum number of subjects (default ...0.8)
% Output:
%   power.img 
%   N.img
%

%hold off



old_path=pwd;
done=0;
data = [];

while ~done  
   [filename newpath] = uigetfile('*.img','SELECT IMAGES, hit CANCEL when finished',200,200);
   if newpath==0,
      done=1;
   else
	   cd (newpath)
       data = [data; read_img(filename)];
   end   
end
root = filename(1:end-4);
h =read_hdr(sprintf('%s.hdr', root);
cd(old_path);


% Set up the defaults
if nargin == 0,
   bf = 1;     
   alpha = 0.05;
   effect_size = 0;
   power_threshold = 0.8;
end

% if bf ~= 1
% compute NUmber Average Number of Ressels
% (Not yet implemented, so for now, it's just one at a time)
% bf = 1;  % approximation:  Bonferroni correction factor could be the Number of Ressels
% end

N = size(data,1);

current_power = zeros(1, h.xdim * h.ydim * h.zdim);
req_N = current_power;
end_power = current_power;

for v=1:size(data,2)
    
    buf = power_finder(data(:,v), alpha, bf, effect_size, power_threshold);
    
    current_power(v) = buf(1);
    req_N(v) = buf(2);
    end_power(v) = buf(3);
    
end

h.glmax = 1;
h.glmin = 0;
write_hdr('power.hdr',h);
write_img_data('power.img', current_power , h);

write_hdr('end_power.hdr',h);
write_img_data('end_power.img',end_power, h);


h.glmax = N*4;
h.glmin = 0;
h.datatype = 2;
h.bits = 8;
write_hdr('N.hdr',h);
write_img_data('N.img', req_N(:,2),h);



return 

