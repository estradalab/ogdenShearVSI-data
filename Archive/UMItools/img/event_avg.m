function [m,s]=event_avg(indata, onsets, window, sampRatio)
% function [m,s]=event_avg(data,onsets, duration, sampRatio)
%
% data:     is a single ROW vector with the time series
% onsets:   is a vector of the onset times (in scan units)
% window:   is the number of scans to average in each event
% sampRatio:   is the interpolation factor if the window, or the 
%               onset times are non-integer.
% result    : is a matrix containing the average single event, and their
%			standard deviation.  The units are % change realtive to the
%			mean signal before detrending
%
% the first scan is scan #1, not #0
%

%datamean=mean(data);
%data = detrend(data);


if sampRatio > 1
    fprintf('Resampling the data by ... %02f \n', sampRatio);
    % I'm finding some interpolators do better than others
    % depending on the ocasion ....
    %data = resample(data,sampRatio,1);
    x=[0:length(indata)-1];
    xi = [0:sampRatio*length(indata)-1]/ sampRatio;
    data = interp1(x, indata, xi);
else
    data=indata;
end

% make sure that the data are in a ROW:
data=reshape(data,length(data),1);

onsets = onsets*sampRatio;
window = window*sampRatio;

result=[];

for n=1:length(onsets)
    %     onsets(n)+window 
    %     max(size(data))
    if (onsets(n)+window <= length(data) )   
        r = data(onsets(n) : onsets(n) + window-1 );
        %r=r-r(1);
        %plot(r); pause
        result = [result ; r'];
    end
end

if ~exist('r') 
    fprintf('No data in buffer.  Check onsets times');
end

m = mean(result,1)';
s = std(result,0,1)';

allevents = result;
save allevents.mat allevents

% m=m*100/datamean;
% s=s*100/datamean;
if sampRatio > 1
    %m = resample(m,1,sampRatio);
    %s = resample(s,1, sampRatio);
    %m = decimate(m,sampRatio);
    %s = decimate(s, sampRatio);
    %m = m(1:sampRatio:end);
    %s = s(1:sampRatio:end);
    xi = [0:window-1];
    x = [0:window/sampRatio-1]*sampRatio;
    m = interp1(xi,m,x);
    s = interp1(xi,s,x);
    
end

return
