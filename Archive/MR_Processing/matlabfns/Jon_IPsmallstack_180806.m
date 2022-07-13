clear all;
close all;
try
delete(findall(0));
end

    %load('1804_ref_vol_mask.mat');
    load('vol_mask.mat');
    load('DESTE_strains_1059_1115_1131_1037.mat');    
    %load('DESTE_strains_1205_1221_1237_1244.mat');
    %load('DESTE_strains_1319_1336_1352_1359.mat');
        %load('DESTE_strains_1501_1536_1548_1436_quadrupled_lambda.mat');
    z=69;
    
    %temp overwrite of data
    origdata;% = data;
    
    try
        mag = sqrt(real(origdata(:,:,:,3)).^2+imag(origdata(:,:,:,3)).^2);
        mag = mag/max(mag(:));
        %mask = ones(size(mag));
    catch
        origdata = data(:,1:2:end,1:2:end,:);
        mag = sqrt(real(origdata(:,:,:,3)).^2+imag(origdata(:,:,:,3)).^2);
        mag = mag/max(mag(:));
    end
    
    
    [tfx, tfy, tfz] = meshgrid(HIRES.axis2,HIRES.axis1,HIRES.axis3);
%    [dmx, dmy, dmz] = meshgrid(axis2(1:2:end-1), axis1, axis3(1:2:end-1));
    [dmx, dmy, dmz] = meshgrid(origpars(1).axis2,origpars(1).axis1,origpars(1).axis3);
    newmaskmaybe = interp3(tfx,tfy,tfz,double(tf),dmx,dmy,dmz);
    tr = 0.5;
    newmaskmaybe(newmaskmaybe<tr)=nan; newmaskmaybe(newmaskmaybe>=tr) = 1;
    
    mask_ = zeros(size(tf)); mask_(tf) = 1;
    
    %%
    % mask = nan(size(mag));
    % mask(mag>.2) = 1;
    for i = 1:3
        %wu{i} = ndnanfilter((origdata(:,:,:,i)),'hamming',[5 3 3]);
        %wu{i} = ndnanfilter((origdata(:,:,:,i)),'hamming',[5 5 5]);
        wu{i} = origdata(:,:,:,i);
        nu{i} = data(:,:,:,i);
    end
    
    %u = WeightedUnwrap3(angle(wu{1}).*flip(newmaskmaybe,1),ones(size(newmaskmaybe)),flip(newmaskmaybe,1));
    u = WeightedUnwrap3(angle(wu{1}).*newmaskmaybe,ones(size(newmaskmaybe)),newmaskmaybe)*lambda(1)/(2*pi);
    v = WeightedUnwrap3(angle(wu{2}).*newmaskmaybe,ones(size(newmaskmaybe)),newmaskmaybe)*lambda(2)/(2*pi);
    w = WeightedUnwrap3(angle(wu{3}).*newmaskmaybe,ones(size(newmaskmaybe)),newmaskmaybe)*lambda(3)/(2*pi);
    
    Nu = WeightedUnwrap3(angle(nu{1}).*mask,ones(size(mask)),mask)*lambda(1)/(2*pi);
    Nv = WeightedUnwrap3(angle(nu{2}).*mask,ones(size(mask)),mask)*lambda(2)/(2*pi);
    Nw = WeightedUnwrap3(angle(nu{3}).*mask,ones(size(mask)),mask)*lambda(3)/(2*pi);
    
    
    
    filt = 'optimalX';
    [Qij{1,1}(:,:,:), Qij{1,2}(:,:,:), Qij{1,3}(:,:,:)] = complexDivideN(wu{1}.*newmaskmaybe,filt);
    [Qij{2,1}(:,:,:), Qij{2,2}(:,:,:), Qij{2,3}(:,:,:)] = complexDivideN(wu{2}.*newmaskmaybe,filt);
    [Qij{3,1}(:,:,:), Qij{3,2}(:,:,:), Qij{3,3}(:,:,:)] = complexDivideN(wu{3}.*newmaskmaybe,filt);
       
    GUij = cell(3,3);
    delta = [abs(origpars(1).axis2(1)-origpars(1).axis2(2)) abs(origpars(1).axis1(1)-origpars(1).axis1(2)) abs(origpars(1).axis3(1)-origpars(1).axis3(2))];

    for e =1:3
        for g =1:3
            try
            GUij{e,g}(:,:,:) = -Qij{e,g}(:,:,:)/(2*pi)*lambda(e)/(2*delta(g));
            catch
            GUij{e,g}(:,:,:) = -Qij{e,g}(:,:,:)/(2*pi)*lambda/(2*delta(g));
            end
        end
    end
    
    epsij = calculateEpsij(GUij);
    Eij = calculateEij(GUij);
    Em = calculateEm(Eij);

    view3dgui(epsij{1,1}.*newmaskmaybe);
    
    %%
    view3dgui(u.*newmaskmaybe);
    view3dgui(v.*newmaskmaybe);
    view3dgui(w.*newmaskmaybe);
    

    
    fu = ndnanfilter(u,'hamming',[3 3 3]);
    fv = ndnanfilter(v,'hamming',[3 3 3]);
    fw = ndnanfilter(w,'hamming',[3 3 3]);
    
    %Give a histogram of the noise in the mask region
    figure(1);
    histogram(log10(mag),-3:0.01:0); hold on;
    histogram(log10(mag(~isnan(newmaskmaybe))),-3:0.01:0); hold on;

    
    figure(2);
    subplot(1,3,1);
    histogram(u,-1:0.05:1); hold on;
    subplot(1,3,2);
    histogram(v,-1:0.05:1); hold on;
    subplot(1,3,3);
    histogram(w,-1:0.05:1); hold on;
    
    %figure(3);
    %mesh(angle(wu{1}(:,:,2*z).*flip(newmaskmaybe(:,:,2*z),1))); 

%%
view3dgui(fu.*newmaskmaybe);
view3dgui(fv.*newmaskmaybe);
view3dgui(fw.*newmaskmaybe);
view3dgui(log10(mag).*newmaskmaybe);

view3dgui(((fu.^2+fv.^2+fw.^2).^0.5).*newmaskmaybe);

% %Histogram of the total sample and the cropped region
% figure; histogram(log10(mag),-3:0.01:0); hold on;
% histogram(log10(mag.*newmaskmaybe),-3:0.01:0);

%view3dgui(flip(newmaskmaybe,1));
%view3dgui(angle(wu{3}).*flip(newmaskmaybe,1));
%view3dgui(HIRES.magnitude);
%view3dgui(mag);