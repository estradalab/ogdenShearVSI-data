clear; close all;

saving = false; 

%% Case 1: linear displacements

%total number of wraps
numwraps = 8;%[1 2 4 8 16];
%-6.5 is a max width of ~300 -> -26.5 is ~3000, 13.5 is ~30.
SNRs = -16.5; %-26.5:2.5:13.5;
%indices for optimal filters to test
filtnum = 3:2:19;

Isize = [256 96];

%Initializing vars
%vars of size Image x numwraps x SNRs x num of filters
u1c1 = zeros([Isize length(numwraps) length(SNRs) 1+length(filtnum)]);
U1 = u1c1;

%vars of size numwraps x SNRs x num of filters
u1c1std = zeros([length(numwraps) length(SNRs) 1+length(filtnum)]);
u1c1mean = u1c1std;
%u0std = u1c1mean;
u1c1error = u1c1std;

u1c1NN = zeros([Isize length(numwraps)]);
U1NN = u1c1NN;

%vars of size numwraps x 1
u1c1meanNN = zeros([length(numwraps),1]);
u1c1stdNN = zeros([length(numwraps),1]);

tic
for i=1:length(numwraps)
    sc = numwraps(i)*2;
    x1_ = 4000*sin(linspace(0,sc*pi,Isize(1)));
    y1_ = 4000*cos(linspace(0,sc*pi,Isize(1)));
    %2D
    x1o = repmat(x1_,[Isize(2) 1])';
    y1o = repmat(y1_,[Isize(2) 1])';
    CI0 = y1o+1i*x1o;
    
    %even shot noise
    % x1 = x1_+2*(2*rand(1,256)-1); y1 = y1_+2*(2*rand(1,256)-1);
    % for i=2:96
    % x1 = cat(1,x1,x1_+2*(2*rand(1,256)-1));
    % y1 = cat(1,y1,y1_+2*(2*rand(1,256)-1));
    % end
    
    %Check of the non-noisy data
    filt = 'optimal3';
    % [GPhij(:,:,1), GPhij(:,:,2)] = complexDivideN(CI0,filt);
    % GUij = -angle(GPhij);
    [GUij(:,:,1), GUij(:,:,2)] = complexDivideN(CI0,filt);
    % figure(110); subplot(1,2,1);
    % imagesc(GUij(:,:,1)); colormap cividis; axis image; caxis([0.10 0.30]);
    % title(['complexDivideN - e11, ' filt]);
    U1NN(:,:,i) = cumsum(GUij(:,:,1),1);
    % subplot(1,2,2)
    % imagesc(U1(:,:,1));  colormap cividis; axis image;
    % title(['complexDivideN - u1, ' filt]);
    u1c1NN(:,:,i) = GUij(:,:,1);
    % save(['nonnoisy_' filt '.mat'], 'CI1', 'GUij', 'U1', 'u1c1');
    %u1c1NN_ = u1c1NN(u1c1NN>sc/8*0.05 & u1c1NN<sc/8*0.35);
    u1c1NN_ = u1c1NN(2:end,2:end,i);
    u1c1meanNN(i) = nanmean(u1c1NN_(:));
    u1c1stdNN(i) = nanstd(u1c1NN_(:));%/u0mean
   
    for j=1:length(SNRs)
        SNR = SNRs(j);
        
        x1 = gNoise(x1o, SNR);
        y1 = gNoise(y1o, SNR);
        CI1 = y1+1i*x1;
        
        % figure (10); hold on;
        % subplot(1,2,1);
        % imagesc(x1); colormap cividis; axis image;
        % subplot(1,2,2)
        % imagesc(y1); colormap cividis; axis image;
        % figure(11);
        % plot(y1(:),x1(:), '.'); xlim([-5000 5000]); ylim([-5000 5000]); axis image;
        
        %% Uli's code as a benchmark
        %CI = cat(4,CI1,CI1,CI1);
        for k=1
            si=size(CI1);
            complex_strainimage=zeros([si 2]);
            CSI = complex_strainimage;
            ind1=2:si(1)-1;
            ind2=2:si(2)-1;
            %ind3=[2:si(3)-1];
            
            phaseimage=CI1./abs(CI1);
            
            % RO encode  lambda  subtraction
            % complex_strainimage(ind1,:,:,1,1)=(phaseimage(ind1+1,:,:,1)./phaseimage(ind1-1,:,:,1));
            % complex_strainimage(:,ind2,:,1,2)=(phaseimage(:,ind2+1,:,1)./phaseimage(:,ind2-1,:,1));
            % complex_strainimage(:,:,ind3,1,3)=(phaseimage(:,:,ind3+1,1)./phaseimage(:,:,ind3-1,1));
            complex_strainimage(ind1,:,1)=(phaseimage(ind1+1,:)./phaseimage(ind1-1,:));
            complex_strainimage(:,ind2,2)=(phaseimage(:,ind2+1)./phaseimage(:,ind2-1));
            SI = angle(complex_strainimage);
            U1(:,:,i,j,k) = cumsum(SI(:,:,1),1);
            
            % figure(1)
            % subplot(2,5,1)
            % imagesc(SI(:,:,1));  colormap cividis; axis image; caxis([0.10 0.30]);
            % title('complexDivideManual - e11');
            % subplot(2,5,6)
            % imagesc(U1(:,:,1));  colormap cividis; axis image;
            % title('complexDivideManual - u1');
            
            u1c1(:,:,i,j,k) = SI(:,:,1);
            %save(['filtertest_CDManual.mat'], 'CI1', 'complex_strainimage', 'SI', 'U1', 'u1c1');
            u1c1_ = u1c1(2:(end-1),2:(end-1),i,j,k);
            %u1c1_ = u1c1(u1c1_>sc/8*0.05 & u1c1_<sc/8*0.35);
            u1c1std(i,j,k) = nanstd(u1c1_(:));%/nanmean(u1c1_(:))
            u1c1error(i,j,k) = 100*sqrt(sum(((u1c1_(:)-u1c1meanNN(i))/u1c1meanNN(i)).^2)/numel(u1c1_));
        end
        %figure(2)
        %plot(1:256,SI(:,48,1)); hold on;
        
        %% Automatic version
        for k = 1:length(3:2:19)
            filt = ['optimal' num2str(filtnum(k))];
            % [GPhij(:,:,1), GPhij(:,:,2)] = complexDivideN(CI0,filt);
            % GUij = -angle(GPhij);
            [GUij(:,:,1), GUij(:,:,2)] = complexDivideN(CI1,filt);
            
            %1 is reserved for Uli's method
            U1(:,:,i,j,k+1) = cumsum(GUij(:,:,1),1);
            u1c1(:,:,i,j,k+1) = GUij(:,:,1);
            
            hfilt = ceil(filtnum(k)/2);
            u1c1_ = u1c1(hfilt:(end-hfilt),hfilt:(end-hfilt),i,j,k+1);
            %u1c1_ = u1c1(u1c1_>sc/8*0.05 & u1c1_<sc/8*0.35);
            u1c1std(i,j,k+1) = nanstd(u1c1_(:));%/u0mean
            u1c1mean(i,j,k+1) = mean(u1c1_(:));
            u1c1error(i,j,k+1) = 100*sqrt(sum(((u1c1_(:)-u1c1meanNN(i))/u1c1meanNN(i)).^2)/numel(u1c1_));
            
            %figure(i*100+j);
            %ksdensity(u1c1_(:)); hold on;
            %save(['filtertest_' filt '.mat'], 'CI1', 'GUij', 'U1', 'u1c1');
        end
        %save('stdU1c1.mat','stdU1c1','u0mean');
    end
end
toc

if saving
    save('Case1-filterCompare.mat','numwraps','SNRs','filtnum','x1o','y1o','u1c1mean','u1c1meanNN','u1c1std','u1c1stdNN','u1c1','u1c1NN','U1','u1c1error')
end

magmavals = colormap(magma);

for j=1:5 
figure; 
    for k=1:size(u1c1error,3) 
        plot(1:2.5:41,u1c1error(j,:,k),'LineWidth',2,'Color',magmavals(1+floor((size(u1c1error,3)-k)/size(u1c1error,3)*64),:)); hold on; 
    end
    set(gca,'ylim',[0 300]);
    line([21 21], [0 300],'LineWidth',1); 
    title(['Number of wraps: ' num2str(numwraps(j))]);
    xlabel('Signal-to-noise ratio [dB]'); ylabel('Noise-induced error in u_{1,1} [%]')
end
%% Manual version
%
% filt = 'optimal3';
% % [GPhij(:,:,1), GPhij(:,:,2)] = complexDivideN(CI0,filt);
% % GUij = -angle(GPhij);
% [GUij(:,:,1), GUij(:,:,2)] = complexDivideN(CI1,filt);
% figure(1); subplot(2,5,2);
% imagesc(GUij(:,:,1)); colormap cividis; axis image; caxis([0.10 0.30]);
% title(['complexDivideN - e11, ' filt]);
% U1 = cumsum(GUij(:,:,1),1);
% subplot(2,5,7)
% imagesc(U1(:,:,1));  colormap cividis; axis image;
% title(['complexDivideN - u1, ' filt]);
% u1c1 = GUij(:,:,1);
% save(['filtertest_' filt '.mat'], 'CI1', 'GUij', 'U1', 'u1c1');
%
% figure(2)
% plot(1:256,GUij(:,48,1)); hold on;
% figure(3);
% %histogram(u1c1(:),0.05:0.005:0.35); hold on;
% u1c1_ = u1c1(u1c1>sc/8*0.05 & u1c1<sc/8*0.35);
% stdU1c1(2) = nanstd(u1c1_)%/u0mean
%
% ksdensity(u1c1_(:)); hold on;
%
% filt = 'optimal5';
% % [GPhij(:,:,1), GPhij(:,:,2)] = complexDivideN(CI0,filt);
% % GUij = -angle(GPhij);
% [GUij(:,:,1), GUij(:,:,2)] = complexDivideN(CI1,filt);
% figure(1); subplot(2,5,3);
% imagesc(GUij(:,:,1)); colormap cividis; axis image; caxis([0.10 0.30]);
% title(['complexDivideN - e11, ' filt]);
% U1 = cumsum(GUij(:,:,1),1);
% subplot(2,5,8)
% imagesc(U1(:,:,1));  colormap cividis; axis image;
% title(['complexDivideN - u1, ' filt]);
% u1c1 = GUij(:,:,1);
% save(['filtertest_' filt '.mat'], 'CI1', 'GUij' ,'U1', 'u1c1');
%
% figure(2)
% plot(1:256,GUij(:,48,1)); hold on;
% figure(3);
% %histogram(u1c1(:),0.05:0.005:0.35); hold on;
% u1c1_ = u1c1(u1c1>sc/8*0.05 & u1c1<sc/8*0.35);
% stdU1c1(3) = nanstd(u1c1_)%/u0mean
% ksdensity(u1c1_(:));
%
% filt = 'optimal7';
% % [GPhij(:,:,1), GPhij(:,:,2)] = complexDivideN(CI0,filt);
% % GUij = -angle(GPhij);
% [GUij(:,:,1), GUij(:,:,2)] = complexDivideN(CI1,filt);
% figure(1); subplot(2,5,4);
% imagesc(GUij(:,:,1)); colormap cividis; axis image; caxis([0.10 0.30]);
% title(['complexDivideN - e11, ' filt]);
% U1 = cumsum(GUij(:,:,1),1);
% subplot(2,5,9)
% imagesc(U1(:,:,1));  colormap cividis; axis image;
% title(['complexDivideN - u1, ' filt]);
% u1c1 = GUij(:,:,1);
% save(['filtertest_' filt '.mat'], 'CI1', 'GUij' ,'U1', 'u1c1');
%
% figure(2)
% plot(1:256,GUij(:,48,1)); hold on;
% figure(3);
% %histogram(u1c1(:),0.05:0.005:0.35); hold on;
% u1c1_ = u1c1(u1c1>sc/8*0.05 & u1c1<sc/8*0.35);
% stdU1c1(4) = nanstd(u1c1_)%/u0mean
% ksdensity(u1c1_(:));
%
% filt = 'optimal9';
% % [GPhij(:,:,1), GPhij(:,:,2)] = complexDivideN(CI0,filt);
% % GUij = -angle(GPhij);
% [GUij(:,:,1), GUij(:,:,2)] = complexDivideN(CI1,filt);
% figure(1); subplot(2,5,5);
% imagesc(GUij(:,:,1)); colormap cividis; axis image; caxis([0.10 0.30]);
% title(['complexDivideN - e11, ' filt]);
% U1 = cumsum(GUij(:,:,1),1);
% subplot(2,5,10)
% imagesc(U1(:,:,1));  colormap cividis; axis image;
% title(['complexDivideN - u1, ' filt]);
% u1c1 = GUij(:,:,1);
% save(['filtertest_' filt '.mat'], 'CI1', 'GUij', 'U1', 'u1c1');
%
% figure(2)
% plot(1:256,GUij(:,48,1)); hold on;
% figure(3);
% %histogram(u1c1(:),0.05:0.005:0.35); hold on;
% u1c1_ = u1c1(u1c1>sc/8*0.05 & u1c1<sc/8*0.35);
% stdU1c1(5) = nanstd(u1c1_)%/u0mean
% ksdensity(u1c1_(:));
% set(gca,'XLim',[sc/8*0.05 sc/8*0.35]);

%% Case 2: linear plus increasing freq sinusoid displacements,

% x1_ = [zeros(1,16), sin(linspace(0,2*pi,256-32))*1000, zeros(1,16)]
% y1_ = [zeros(1,16), cos(linspace(0,2*pi,256-32))*1000, zeros(1,16)];
% x1 = repmat(x1_,[96 1])';
% y1 = repmat(y1_,[96 1])';
%
%
% ph = angle(y1+1i*x1);
% figure; imagesc(ph); axis image;
