%The goal of this script is to determine what role Gaussian noise in the
%independent lambda+ and lambda- experiments plays in total error in
%strain

clear; close all;

saving = false;
plots = true;
rbVals = flipud(brewermap(256,'RdYlBu'));

%% Synthetic : linear displacements

allSNRs = 4:2:20;
magmaVals = colormap(magma); 

filtnums = 5;%3:2:19;
for q = 1:length(filtnums)
    for p = 1:length(allSNRs)
        
        %total number of wraps
        numwraps = 8;%[1 2 4 8 16];
        %-6.5 is a max width of ~300 -> -26.5 is ~3000, 13.5 is ~30.
        SNRs = allSNRs(p);%-16.5; %-26.5:2.5:13.5;
        %indices for optimal filters to test
        filtnum = filtnums(q);
        
        Isize = [256 2048];
        
        phi = [rand() rand()];
        ampmax = 2;
        rng = 10.^(linspace(log10(sqrt(ampmax)),log10(sqrt(1/ampmax)),Isize(1))); invrng = 1./rng;
        sc = numwraps*2;
        amp = 4000;
        x1min_ = amp*sin(linspace(0,-sc*pi/2+phi(1),Isize(1))).*rng;
        y1min_ = amp*cos(linspace(0,-sc*pi/2+phi(1),Isize(1))).*rng;
        x1plu_ = amp*sin(linspace(0,sc*pi/2+phi(2),Isize(1))).*invrng;
        y1plu_ = amp*cos(linspace(0,sc*pi/2+phi(2),Isize(1))).*invrng;
        
        %2D
        x1mo = repmat(x1min_,[Isize(2) 1])';
        y1mo = repmat(y1min_,[Isize(2) 1])';
        CIm0 = y1mo+1i*x1mo;
        
        x1po = repmat(x1plu_,[Isize(2) 1])';
        y1po = repmat(y1plu_,[Isize(2) 1])';
        CIp0 = y1po+1i*x1po;
        
        
        %Check of the non-noisy data
        % filt = 'optimal3';
        % [GUij(:,:,1), GUij(:,:,2)] = complexDivideN(CIm0,filt);
        % U1NN = cumsum(GUij(:,:,1),1);
        % u1c1NN = GUij(:,:,1);
        
        x1m = gNoise(x1mo, SNRs, amp);
        y1m = gNoise(y1mo, SNRs, amp);
        CIm1 = y1m+1i*x1m;
        
        x1p = gNoise(x1po, SNRs, amp);
        y1p = gNoise(y1po, SNRs, amp);
        CIp1 = y1p+1i*x1p;
        
        if plots
        figure(1); subplot(1,2,1); plot(x1m(:),y1m(:),'.'); hold on; axis image;
        title('\lambda_-'); xlabel('Real(Z_-)'); ylabel('Imag(Z_-)');
        subplot(1,2,2); plot(x1p(:),y1p(:),'.'); axis image;
        title('\lambda_+'); xlabel('Real(Z_+)'); ylabel('Imag(Z_+)');
        
        figure(201); subplot(1,2,1); imagesc(angle(CIm0)); axis image;
        title('CI Phase, \lambda_-');
        subplot(1,2,2); imagesc(angle(CIp0)); axis image;
        title('CI Phase, \lambda_+');  colormap(rbVals); truesize;
        figure(202);
        subplot(1,2,1); imagesc(abs(CIm0)); axis image;
        title('CI Magnitude, \lambda_-');
        subplot(1,2,2); imagesc(abs(CIp0)); axis image;
        title('CI Magnitude, \lambda_+'); colormap gray; truesize;
        
        figure(211); subplot(1,2,1); imagesc(angle(CIm1)); axis image;
        title('CI Phase, \lambda_-');
        subplot(1,2,2); imagesc(angle(CIp1)); axis image;
        title('CI Phase, \lambda_+');  colormap(rbVals); truesize;
        figure(212);
        subplot(1,2,1); imagesc(abs(CIm1)); axis image;
        title('CI Magnitude, \lambda_-');
        subplot(1,2,2); imagesc(abs(CIp1)); axis image;
        title('CI Magnitude, \lambda_+'); colormap gray; truesize;
        end
        
        % CI0 = (CIp0./CIm0)./abs(CIp0./CIm0).*(abs(CIp0)+abs(CIm0))/2;
        % CI1 = (CIp1./CIm1)./abs(CIp1./CIm1).*(abs(CIp1)+abs(CIm1))/2;
        CI0 = (CIp0./CIm0)./abs(CIp0./CIm0).*sqrt((abs(CIp0).*abs(CIm0)));
        CI1 = (CIp1./CIm1)./abs(CIp1./CIm1).*sqrt((abs(CIp1).*abs(CIm1)));
        
        if plots
        figure(2); imagesc(angle(CIp0./CIm0)); axis image; truesize; 
        colormap(rbVals);
            
        figure(30); imagesc(sqrt(abs(CIp1).*abs(CIm1))); axis image;
        title('Weighted Amplitude'); colorbar; colormap(rbVals); truesize;
        figure(31); imagesc(angle(CI1)); axis image;
        title('QC Phase'); colorbar; colormap(rbVals); truesize;
        figure(32); 
        imagesc(abs(CIp1./CIm1)); axis image; colorbar; colormap gray; caxis([1/ampmax ampmax]);
        title('QC Amplitude'); truesize;
        end
        % subplot(1,2,1); imagesc(sqrt((abs(CIp0).*abs(CIm0)))); axis image; colorbar;
        % title('Mean Amplitude');
        
        filt = ['optimal' num2str(filtnum)];
            [GUij0(:,:,1), GUij0(:,:,2)] = complexDivideN(CI0,filt);
            [GUij1(:,:,1), GUij1(:,:,2)] = complexDivideN(CI1,filt);
            if plots
            figure(30)
            imagesc(GUij1(:,:,1)); colorbar; axis image;
            title('u_{1,1}');
            end
            
            GU11 = GUij0(:,:,1); mean_u1c1 = nanmean(GU11(:));
            GU22 = GUij0(:,:,2); mean_u2c2 = nanmean(GU22(:));

            
            err_u1c1 = GUij1(:,:,1)-GU11;
            if plots
            figure(9); subplot(1,2,1); imagesc(err_u1c1); colorbar; colormap(rbVals); axis image;
            title('error in u_{1,1}'); caxis([-0.4 0.4]); hold on;
            subplot(1,2,2); ksdensity(err_u1c1(:)); xlim([-0.2 0.2]); hold on;
            end
            
            for i=1:Isize(1)
            u1c1_std(i) = std(err_u1c1(i,:));
            end
            crop = 5;
            if plots
            figure(52); plot(log10(rng((1+crop):end-crop).^2), u1c1_std((1+crop):end-crop),'.','Color',magmaVals(round(p/length(allSNRs)*63+1),:));
            end
            
            log10rng = log10(rng((1+crop):end-crop)); u1c1_std_ = (u1c1_std((1+crop):end-crop));
            hold on;
            cf = fit(log10rng', u1c1_std_', 'a*x^2+c');
            
            n_amp = (amp / SNRs);
            SNRpts = min([amp*rng; amp*invrng]/n_amp);
            
            figure(ampmax*10000+filtnum*100+numwraps);
            plot((SNRpts((1+crop):end-crop)), u1c1_std_,'.','Color',magmaVals(round(p/length(allSNRs)*63+1),:)); hold on;
            xlabel('Minimum SNR'); ylabel('Standard deviation strain error');
            
            cropGU11 = GUij1(crop:end-crop, :,1);
            stdu1c1(p) = std(cropGU11(:));
    end
end

figure; plot(allSNRs, stdu1c1, '.');