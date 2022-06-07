    FWIN = 'hamming';
    F = [5 5];
    N = 100;
    Pnoise = 0.30;
    PNaNs  = 0.20;
    X = peaks(N);                                     % original
    Y = X + ((rand(size(X))-0.5)*2)*max(X(:))*Pnoise; % add noise
    Y(round(1 + (N^2-1).*rand(N^2*PNaNs,1))) = NaN;   % add NaNs
    [Z0,W] = ndnanfilter(Y,FWIN,F);                   % filters
    Z1 = Z0; Z2 = Y; inan = isnan(Y);
    Z1(inan) = NaN;
    Z2(inan) = Z0(inan);  
    subplot(231), imagesc(X), clim = caxis; axis equal tight
                  title('Original data')
    subplot(232), imagesc(Y),  caxis(clim), axis equal tight 
                  title('Data + NOISE + NaNs')
    subplot(234), imagesc(Z0), caxis(clim), axis equal tight 
                  title('FILTERS + NaNs interpolation')
    subplot(235), imagesc(Z1), caxis(clim), axis equal tight 
                  title('FILTERS ignoring NaNs')
    subplot(236), imagesc(Z2), caxis(clim), axis equal tight 
                  title('GAP-filling with interpolated NaNs')
    subplot(233), imagesc(-F(1):F(1),-F(2):F(2),W), axis equal tight, 
                   title([upper(FWIN) ' 2D window']), view(2) 