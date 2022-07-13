function mkMIP(img, stdT, mskImg)
% function mkMIP(img, stdT, mskImg)
% 
% Display crude maximum intensity projections of THRESHOLDED
% 3D brain maps
%
% img:  the  data
% stdT:  How many standard deviations to use for the threshold

DEBUG = 0;
    
% discard the zeros, and take the abs()
img = abs(img);
buffer = img(:);
% remove the zeros from the data
buffer = buffer(find(buffer));

% basic stats on the image
m = mean(buffer(:));
sd = std(buffer(:));

threshold = m + stdT *sd;

if DEBUG
    fprintf('\nInfo about MIP: -- \nMean: %6.2f  | StdDev: %6.2f | Threshold: %6.2f | Max:  %6.2f | Min: %6.2f\n' , ...
        m, sd, threshold, max(img(:)), min(img(:)) );
end

% we use 10% of the mean to create a background for overlaying the data
% buffer = zeros(size(img));
% buffer(find(img)) = m/10;

img(find(img < threshold)) = 0;

% do the MIP here
f1 = squeeze(sum(img,1));   f1 =  f1 / max(f1(:));
f1m = squeeze(sum(mskImg,1));  f1m = f1m / max(f1m(:));

f2 = squeeze(sum(img,2));    f2 = f2 / max(f2(:));
f2m = squeeze(sum(mskImg,2));   f2m = f2m / max(f2m(:));

f3 = squeeze(sum(img,3));    f3 = f3 / max(f3(:));
f3m = squeeze(sum(mskImg,3));     f3m = f3m / max(f3m(:));

f1m (find(f1)) = 1 + f1(find(f1));
f2m (find(f2)) = 1 + f2(find(f2) );
f3m (find(f3)) = 1 + f3(find(f3));

colormap([gray ; jet])

set(gcf, 'Position',[100 100 800 400])
subplot (231), imagesc(f1m'), axis xy, axis equal, axis tight
subplot (232), imagesc(f2m'), axis xy, axis equal, axis tight
subplot (233), imagesc(f3m'), axis xy, axis equal, axis tight


return