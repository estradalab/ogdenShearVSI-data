%Run this after runComplexUnwrap
clear all; close all;

mu = 1E6;
a = 0.5;
%mu2 = 1E6;
K = 1E7;

%% Calculate experimental output sensitivity fields for a single timepoint

dfil = glob('*1656*');
load(dfil{1});

mfil = glob('mechvars*');
load(mfil{1});
[hx,hy,hz] = meshgrid(HIRES.axis2,HIRES.axis1,HIRES.axis3);

validIdxs = find(~isnan(Eij{1,1}) & ~isnan(Eij{2,2}) & ~isnan(Eij{3,3}));
%experimental coordinates
eCOORD = [hx(validIdxs) hy(validIdxs) hz(validIdxs)];

elam123 = zeros(size(eCOORD));
eU11 = 2*Eij{1,1}(validIdxs) + 1;
eU22 = 2*Eij{2,2}(validIdxs) + 1;
eU33 = 2*Eij{3,3}(validIdxs) + 1;
eU12 = 2*Eij{1,2}(validIdxs);
eU13 = 2*Eij{1,3}(validIdxs);
eU23 = 2*Eij{2,3}(validIdxs);
eU = [eU11, eU22, eU33, eU12, eU13, eU23];


ed2PidKdF = zeros(length(eU),3,3,3,3);
ed2PidMudF = zeros(length(eU),3,3,3,3);
ed2PidadF = zeros(length(eU),3,3,3,3);
ed2PieqdKdF_deltaF = zeros(length(eU),3,3);
ed2PieqdMudF_deltaF = zeros(length(eU),3,3);
ed2PieqdadF_deltaF = zeros(length(eU),3,3);
for k=1:length(eU)
    U_ = [eU(k,1) eU(k,4) eU(k,5); eU(k,4) eU(k,2) eU(k,6); eU(k,5) eU(k,6) eU(k,3)];
    lam123(k,:) = sort(eig(U_)','descend');
    J_ = det(U_);
    
    
    ed2PidKdF(k,:,:,:,:) = d2PidKdF_NH(U_,mu,K);
    ed2PidMudF(k,:,:,:,:) = d2PidmudF_MR(U_,mu,a,K);
    %ed2PidMu2dF(k,:,:,:,:) = d2PidMu2dF_Yeoh2(U_,mu,mu2,K);
    ed2PidadF(k,:,:,:,:) = d2PidadF_MR(U_,mu,a,K);

    for n=1:3
        for m=1:3
            ed2PieqdKdF_deltaF(k,:,:) = ed2PieqdKdF_deltaF(k,:,:)+squeeze(ed2PidKdF(k,:,:,n,m))*J_^(-1/3)*(0.00001*(U_(n,m)-dij(n,m))+dij(n,m));%*dij(n,m));%
        end
    end
    
    
    %csum = zeros(3);
    for n=1:3
        for m=1:3
            %csum = csum + squeeze(d2PidMudF(:,:,n,m,k)*0.00001*U_(n,m));
            ed2PieqdMudF_deltaF(k,:,:) = ed2PieqdMudF_deltaF(k,:,:)+squeeze(ed2PidMudF(k,:,:,n,m))*J_^(-1/3)*(0.00001*(U_(n,m)-dij(n,m))+dij(n,m));
        end
    end
    %d2PieqdMudF_deltaF(k) = sqrt(sum(sum(csum.*csum)));
    
    %csum = zeros(3);
    for n=1:3
        for m=1:3
            %csum = csum + squeeze(d2PidMu2dF(:,:,n,m,k)*0.00001*U_(n,m));
            ed2PieqdadF_deltaF(k,:,:) = ed2PieqdadF_deltaF(k,:,:)+squeeze(ed2PidadF(k,:,:,n,m))*J_^(-1/3)*(0.00001*(U_(n,m)-dij(n,m))+dij(n,m));
        end
    end
    
end

figure; histogram(ed2PieqdMudF_deltaF(:,1,1),'EdgeColor','none'); hold on;
histogram(ed2PieqdMudF_deltaF(:,2,2),'EdgeColor','none'); 
histogram(ed2PieqdMudF_deltaF(:,3,3),'EdgeColor','none');
histogram(ed2PieqdMudF_deltaF(:,1,2),'EdgeColor','none'); 
histogram(ed2PieqdMudF_deltaF(:,1,3),'EdgeColor','none'); 
histogram(ed2PieqdMudF_deltaF(:,2,3),'EdgeColor','none');
title(['Stress components sensitivity to \mu, histogram, \alpha = ' num2str(a)])
L = [-2.5e-5 2.5E-5];
xlim(L)

%Plot of stress sensitivity d^2Pi_ij/dmu dF_km*deltaF_km
%L = [-2.5e-5 2.5E-5];
figure; subplot(3,3,1);
scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)),  4, ed2PieqdMudF_deltaF(:,1,1),'filled','Marker','o'); axis image; caxis(L);
subplot(3,3,5);
scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)),  4, ed2PieqdMudF_deltaF(:,2,2),'filled','Marker','o'); axis image; caxis(L);
subplot(3,3,9);
scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)),  4, ed2PieqdMudF_deltaF(:,3,3),'filled','Marker','o'); axis image; caxis(L);

subplot(3,3,2);
scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)),  4, ed2PieqdMudF_deltaF(:,1,2),'filled','Marker','o'); axis image; caxis(L);
subplot(3,3,3);
scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)),  4, ed2PieqdMudF_deltaF(:,1,3),'filled','Marker','o'); axis image; caxis(L);
subplot(3,3,6);
scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)),  4, ed2PieqdMudF_deltaF(:,2,3),'filled','Marker','o'); axis image; caxis(L);

%Log10 of the absolute value of above
L_ = [log10(L(2))-2 log10(L(2))];
figure; subplot(3,3,1);
scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)),  4, log10(abs(ed2PieqdMudF_deltaF(:,1,1))),'filled','Marker','o'); axis image; caxis(L_);
subplot(3,3,5);
scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)),  4, log10(abs(ed2PieqdMudF_deltaF(:,2,2))),'filled','Marker','o'); axis image; caxis(L_);
subplot(3,3,9);
scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)),  4, log10(abs(ed2PieqdMudF_deltaF(:,3,3))),'filled','Marker','o'); axis image; caxis(L_);

subplot(3,3,2);
scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)),  4, log10(abs(ed2PieqdMudF_deltaF(:,1,2))),'filled','Marker','o'); axis image; caxis(L_);
subplot(3,3,3);
scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)),  4, log10(abs(ed2PieqdMudF_deltaF(:,1,3))),'filled','Marker','o'); axis image; caxis(L_);
subplot(3,3,6);
scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)),  4, log10(abs(ed2PieqdMudF_deltaF(:,2,3))),'filled','Marker','o'); axis image; caxis(L_);


%Plot of stress sensitivity d^2Pi_ij/dK dF_km*deltaF_km
figure; subplot(3,3,1);
scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)),  4, ed2PieqdKdF_deltaF(:,1,1),'filled','Marker','o'); axis image; caxis(L);
subplot(3,3,5);
scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)),  4, ed2PieqdKdF_deltaF(:,2,2),'filled','Marker','o'); axis image; caxis(L);
subplot(3,3,9);
scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)),  4, ed2PieqdKdF_deltaF(:,3,3),'filled','Marker','o'); axis image; caxis(L);

subplot(3,3,2);
scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)),  4, ed2PieqdKdF_deltaF(:,1,2),'filled','Marker','o'); axis image; caxis(L);
subplot(3,3,3);
scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)),  4, ed2PieqdKdF_deltaF(:,1,3),'filled','Marker','o'); axis image; caxis(L);
subplot(3,3,6);
scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)),  4, ed2PieqdKdF_deltaF(:,2,3),'filled','Marker','o'); axis image; caxis(L);

%Plot of stress sensitivity d^2Pi_ij/da dF_km*deltaF_km

La = L*mu;
figure; subplot(3,3,1);
scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)),  4, ed2PieqdadF_deltaF(:,1,1),'filled','Marker','o'); axis image; caxis(La);
subplot(3,3,5);
scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)),  4, ed2PieqdadF_deltaF(:,2,2),'filled','Marker','o'); axis image; caxis(La);
subplot(3,3,9);
scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)),  4, ed2PieqdadF_deltaF(:,3,3),'filled','Marker','o'); axis image; caxis(La);

subplot(3,3,2);
scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)),  4, ed2PieqdadF_deltaF(:,1,2),'filled','Marker','o'); axis image; caxis(La);
subplot(3,3,3);
scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)),  4, ed2PieqdadF_deltaF(:,1,3),'filled','Marker','o'); axis image; caxis(La);
subplot(3,3,6);
scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)),  4, ed2PieqdadF_deltaF(:,2,3),'filled','Marker','o'); axis image; caxis(La);

figure; histogram(ed2PieqdadF_deltaF(:,1,1),'EdgeColor','none'); hold on;
histogram(ed2PieqdadF_deltaF(:,2,2),'EdgeColor','none'); 
histogram(ed2PieqdadF_deltaF(:,3,3),'EdgeColor','none');
histogram(ed2PieqdadF_deltaF(:,1,2),'EdgeColor','none'); 
histogram(ed2PieqdadF_deltaF(:,1,3),'EdgeColor','none'); 
histogram(ed2PieqdadF_deltaF(:,2,3),'EdgeColor','none');
title(['Stress components sensitivity to \alpha, histogram']);
xlim(La);

figure; histogram(ed2PieqdKdF_deltaF(:,1,1),'EdgeColor','none'); hold on;
histogram(ed2PieqdKdF_deltaF(:,2,2),'EdgeColor','none'); 
histogram(ed2PieqdKdF_deltaF(:,3,3),'EdgeColor','none');
histogram(ed2PieqdKdF_deltaF(:,1,2),'EdgeColor','none'); 
histogram(ed2PieqdKdF_deltaF(:,1,3),'EdgeColor','none'); 
histogram(ed2PieqdKdF_deltaF(:,2,3),'EdgeColor','none');
title(['Stress components sensitivity to K, histogram']);

%% 3D sensitivity plots for experimental data

%%Old version, with the absolute magnitude of stress for plotting
% parulamap = colormap(parula);
% figure; subplot(3,1,1);
% eq = ed2PieqdMudF_deltaF;
% eqlim = [mean(eq)-0.2 mean(eq)+0.2];
% epimupct = (eq-eqlim(1))./(eqlim(2)-eqlim(1));
% %pimupct = (q-min(q))./(max(q)-min(q));
% emuidxs = ceil(epimupct*64);
% ecmu = zeros(length(epimupct),3);
% for j=1:length(epimupct)
%     try
%         ecmu(j,:) = parulamap(emuidxs(j),:);
%     catch
%         if emuidxs(j)>=64
%             emuidxs(j) = 64;
%         elseif emuidxs(j)<=0
%             emuidxs(j) = 1;
%         end
%         ecmu(j,:) = parulamap(emuidxs(j),:);
%     end
% end
% scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)),  4, ecmu,'filled','Marker','o'); axis image;
% 
% subplot(3,1,2);
% eq2 = log10(ed2PieqdKdF_deltaF);
% eq2lim = [mean(eq2)-0.2 mean(eq2)+0.2];
% epiKpct = (eq2-eq2lim(1))./(eq2lim(2)-eq2lim(1));
% %piKpct = (q2-min(q2))./(max(q2)-min(q2));
% eKidxs = ceil(epiKpct*64);
% ecK = zeros(length(epiKpct),3);
% for j=1:length(epiKpct)
%     try
%         ecK(j,:) = parulamap(eKidxs(j),:);
%     catch
%         if eKidxs(j)>=64
%             eKidxs(j) = 64;
%         elseif eKidxs(j)<=0
%             eKidxs(j) = 1;
%         end
%         ecK(j,:) = parulamap(eKidxs(j),:);
%     end
% end
% scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)), 4, ecK,'filled','Marker','o'); axis image;
% 
% subplot(3,1,3);
% eq3 = log10(ed2PieqdadF_deltaF);
% eq3lim = [mean(eq3)-0.2 mean(eq3)+0.2];
% epiMu2pct = (eq3-eq3lim(1))./(eq3lim(2)-eq3lim(1));
% %piMu2pct = (q3-min(q3))./(max(q3)-min(q3));
% eMu2idxs = ceil(epiMu2pct*64);
% ecMu2 = zeros(length(epiMu2pct),3);
% for j=1:length(epiMu2pct)
%     try
%         ecMu2(j,:) = parulamap(eMu2idxs(j),:);
%     catch
%         if eMu2idxs(j)>=64
%             eMu2idxs(j) = 64;
%         elseif eMu2idxs(j)<=0
%             eMu2idxs(j) = 1;
%         end
%         ecMu2(j,:) = parulamap(eMu2idxs(j),:);
%     end
% end
% scatter3(eCOORD(:,1)-median(eCOORD(:,1)), eCOORD(:,2)-median(eCOORD(:,2)), eCOORD(:,3)-median(eCOORD(:,3)), 4, ecMu2,'filled','Marker','o'); axis image;




%% Calculate Abaqus output sensitivity fields for a single timepoint


load('abaqus_181212_Wrench_3p0.mat');
lam123 = zeros(size(COORD));
U = zeros(size(LE));
d2PidKdF = zeros(3,3,3,3,length(LE));
d2PidMudF = zeros(3,3,3,3,length(LE));
d2PidMu2dF = zeros(3,3,3,3,length(LE));
%d2PieqdKdF_deltaF = zeros(length(LE),1);
%d2PieqdMudF_deltaF = zeros(length(LE),1);
%d2PieqdMu2dF_deltaF = zeros(length(LE),1);
d2PieqdKdF_deltaF = zeros(3,3,length(LE));
d2PieqdMudF_deltaF = zeros(3,3,length(LE));
d2PieqdMu2dF_deltaF = zeros(3,3,length(LE));

for k=1:length(LE)
    U_ = expm([LE(k,1) LE(k,4) LE(k,5); LE(k,4) LE(k,2) LE(k,6); LE(k,5) LE(k,6) LE(k,3)]);
    U(k,:) = putMatIntoCol(U_);
    %U(k,1) = U_(1,1); U(k,2) = U_(2,2); U(k,3) = U_(3,3);
    %U(k,4) = U_(1,2); U(k,5) = U_(1,3); U(k,6) = U_(2,3);
    lam123(k,:) = sort(eig(U_)','descend');
    
    B_ = U_*U_;
    B(k,:) = putMatIntoCol(B_);
    J(k) = det(U_);
    
    I1(k) = trace(B_);
    I2(k) = 0.5*(trace(B_)^2-trace(B_*B_));
    
    I1bar(k) = I1(k)*J(k)^(-2/3);
    I2bar(k) = I2(k)*J(k)^(-4/3);
    
    
    d2PidKdF(:,:,:,:,k) = d2PidKdF_NH(U_,mu,K);
    d2PidMudF(:,:,:,:,k) = d2PidMudF_NH(U_,mu,K);
    d2PidMu2dF(:,:,:,:,k) = d2PidMu2dF_Yeoh2(U_,mu,mu2,K);
    
    %csum = zeros(3);
    %for n=1:3
    %    for m=1:3
    %        csum = csum + squeeze(d2PidKdF(:,:,n,m,k)*0.00001*U_(n,m));%*dij(n,m));%
    %    end
    %end
    %d2PieqdKdF_deltaF(k) = sqrt(sum(sum(csum.*csum)));
    for n=1:3
        for m=1:3
            d2PieqdKdF_deltaF(:,:,k) = d2PidKdF(:,:,n,m,k)*0.00001*U_(n,m);%*dij(n,m));%
        end
    end
    
    
    %csum = zeros(3);
    for n=1:3
        for m=1:3
            %csum = csum + squeeze(d2PidMudF(:,:,n,m,k)*0.00001*U_(n,m));
            d2PieqdMudF_deltaF(:,:,k) = d2PidMudF(:,:,n,m,k)*0.00001*U_(n,m);
        end
    end
    %d2PieqdMudF_deltaF(k) = sqrt(sum(sum(csum.*csum)));
    
    %csum = zeros(3);
    for n=1:3
        for m=1:3
            %csum = csum + squeeze(d2PidMu2dF(:,:,n,m,k)*0.00001*U_(n,m));
            d2PieqdMu2dF_deltaF(:,:,k) = d2PidMu2dF(:,:,n,m,k)*0.00001*U_(n,m);
        end
    end
    %d2PieqdMu2dF_deltaF(k) = sqrt(sum(sum(csum.*csum)));
end

figure(30);
plot(I1bar(:),I2bar(:),'.');
xlim([3 3.1]);
ylim([3 3.1]);
hold on;

lam = 0.5:0.01:3;
I1lam = lam.^2+2./lam;
I1uniax_lam = lam.^2+2./lam;
I2uniax_lam = 2*lam+(lam.^-2);
I1ps_lam = 1+lam.^2+(lam.^-2);
I2ps_lam = I1ps_lam;
I1biax_lam = 2*lam.^2+lam.^-4;
I2biax_lam = lam.^4+2./lam.^2;
plot(I1uniax_lam,I2uniax_lam,'-','color','r');
plot(I1ps_lam,I2ps_lam,'-','color','r');
plot(I1biax_lam,I2biax_lam,'-','color','r');




%% 3D sensitivity plots for ABAQUS data
parulamap = colormap(parula);
figure; subplot(3,1,1);
q = log10(d2PieqdMudF_deltaF);
qlim = [mean(q)-0.1 mean(q)+0.1];
pimupct = (q-qlim(1))./(qlim(2)-qlim(1));
%pimupct = (q-min(q))./(max(q)-min(q));
muidxs = ceil(pimupct*64);
cmu = zeros(length(pimupct),3);
for j=1:length(pimupct)
    try
        cmu(j,:) = parulamap(muidxs(j),:);
    catch
        if muidxs(j)>=64
            muidxs(j) = 64;
        elseif muidxs(j)<=0
            muidxs(j) = 1;
        end
        cmu(j,:) = parulamap(muidxs(j),:);
    end
end
scatter3(COORD(:,1)-median(COORD(:,1)), COORD(:,2)-median(COORD(:,2)), COORD(:,3)-median(COORD(:,3)), 20, cmu,'filled','Marker','o'); axis image;

subplot(3,1,2);
q2 = log10(d2PieqdKdF_deltaF);
q2lim = [mean(q2)-0.1 mean(q2)+0.1];
piKpct = (q2-q2lim(1))./(q2lim(2)-q2lim(1));
%piKpct = (q2-min(q2))./(max(q2)-min(q2));
Kidxs = ceil(piKpct*64);
cK = zeros(length(piKpct),3);
for j=1:length(piKpct)
    try
        cK(j,:) = parulamap(Kidxs(j),:);
    catch
        if Kidxs(j)>=64
            Kidxs(j) = 64;
        elseif Kidxs(j)<=0
            Kidxs(j) = 1;
        end
        cK(j,:) = parulamap(Kidxs(j),:);
    end
end
scatter3(COORD(:,1)-median(COORD(:,1)), COORD(:,2)-median(COORD(:,2)), COORD(:,3)-median(COORD(:,3)),  20, cK,'filled','Marker','o'); axis image;

subplot(3,1,3);
q3 = log10(d2PieqdMu2dF_deltaF);
q3lim = [mean(q3)-0.1 mean(q3)+0.1];
piMu2pct = (q3-q3lim(1))./(q3lim(2)-q3lim(1));
%piMu2pct = (q3-min(q3))./(max(q3)-min(q3));
Mu2idxs = ceil(piMu2pct*64);
cMu2 = zeros(length(piMu2pct),3);
for j=1:length(piMu2pct)
    try
        cMu2(j,:) = parulamap(Mu2idxs(j),:);
    catch
        if Mu2idxs(j)>=64
            Mu2idxs(j) = 64;
        elseif Mu2idxs(j)<=0
            Mu2idxs(j) = 1;
        end
        cMu2(j,:) = parulamap(Mu2idxs(j),:);
    end
end
scatter3(COORD(:,1)-median(COORD(:,1)), COORD(:,2)-median(COORD(:,2)), COORD(:,3)-median(COORD(:,3)),  20, cMu2,'filled','Marker','o'); axis image;

