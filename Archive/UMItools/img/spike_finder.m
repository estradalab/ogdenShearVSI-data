function bad_imgs = spike_finder(file_name, doCorrect)
% function bad_imgs = spike_finder(file_name, doCorrect)
%
% this function returns the index number of those images in the time series
% whose mean value is more than three standard deviations different from the
% grand mean
%
% file_name   : the name of any of the images in the time series.
% doCorrect   : replace the suspicious image by the average of its naighbors

data = read_img_series(file_name(1:end-8));

m = mean(data,2);
gm = mean(m);
sd = std(m);


bad_imgs = find( abs(m-gm) > 3*sd );
fprintf('\n Grand Mean ... %f', gm);
fprintf('\n standard deviation ... %f', sd);

if ~isempty(bad_imgs)
    
    fprintf('\n Look out for the following images.');
    fprintf('\n Their mean value deviates more than three times the std. dev');
    bad_imgs'
    m(bad_imgs)'
    
    if(doCorrect)
        
        files = dir(sprintf('%s*.img', file_name(1:end-8)));
        h=read_hdr(sprintf('%s.hdr', file_name(1:end-4)));
        
        for count = 1:length(bad_imgs);
            
            b = bad_imgs(count);
            outname = files(b).name;
            fprintf('\nReplacing ... %s', outname);
            outdata = ( data(b-1,:) + data(b+1,:) )/2;
            write_img(outname, outdata, h); 
            
        end
    end        
end


return