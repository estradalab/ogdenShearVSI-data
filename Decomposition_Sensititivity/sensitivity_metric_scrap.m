clear all;
close all;


% lam = 1.01.^linspace(0,60,301);
lam = linspace(1,1.8,81);
alpha = linspace(-5,5,1001);

%Example parameters from Treloar's data:
%mu1 = 0.62 MPa, alpha1 = 1.3
%mu2 = 0.00118 MPa, alpha2 = 5
%mu3 = -0.00981 MPa, alpha3 = -2

mu = 1;


%Salpha_ps = mu*bsxfun(@times,log(lam'),(bsxfun(@power,lam',alpha-1)-bsxfun(@power,lam',-alpha-1)));
Salpha_ps =  mu*bsxfun(@times,log(lam'),(bsxfun(@power,lam',alpha-1)-bsxfun(@power,lam',-alpha-1)));

%Salpha_ua = mu*bsxfun(@times,log(lam'),(bsxfun(@power,lam',alpha-1)-0.5*bsxfun(@power,lam',-0.5*alpha-1)));
Salpha_ua = mu*bsxfun(@times,log(lam'),(bsxfun(@power,lam',alpha-1)-bsxfun(@power,lam',-0.5*alpha-1)));

%Salpha_ba = mu*bsxfun(@times,log(lam'),(0.5*bsxfun(@power,lam',0.5*(alpha-1))-bsxfun(@power,lam',-alpha-0.5)));
Salpha_ba = mu*bsxfun(@times,log(lam'),(bsxfun(@power,lam',0.5*alpha-1)-bsxfun(@power,lam',-alpha-1)));

%k = 0:0.01:1;
%Old, corrected on 220126
%Salpha_gen = @(k) mu*bsxfun(@times,log(lam'),((-k^2+1.5*k+0.5)*bsxfun(@power,lam',(-k^2+1.5*k+0.5)*alpha-1)...
%    +(-k+0.5)*bsxfun(@power,lam',(-k+0.5)*alpha-1)...
%    +(k^2-0.5*k-1)*bsxfun(@power,lam',(k^2-0.5*k-1)*alpha-1)));
Salpha_gen = @(k) mu*bsxfun(@times,log(lam'),((-k^2+1.5*k+0.5)^2*bsxfun(@power,lam',(-k^2+1.5*k+0.5)*alpha-1)...
    +(-k+0.5)^2*bsxfun(@power,lam',(-k+0.5)*alpha-1)...
    +(k^2-0.5*k-1)^2*bsxfun(@power,lam',(k^2-0.5*k-1)*alpha-1)));
Smu_gen = @(k) (-k^2+1.5*k+0.5)*bsxfun(@power,lam',(-k^2+1.5*k+0.5)*alpha-1)...
    +(-k+0.5)*bsxfun(@power,lam',(-k+0.5)*alpha-1)...
    +(k^2-0.5*k-1)*bsxfun(@power,lam',(k^2-0.5*k-1)*alpha-1);

W_gen = @(k) mu*bsxfun(@times,1./alpha,(bsxfun(@power,lam',(-k^2+1.5*k+0.5)*alpha)...
    +bsxfun(@power,lam',(-k+0.5)*alpha)...
    +bsxfun(@power,lam',(k^2-0.5*k-1)*alpha)-3));

figure;  subplot(1,3,1);
imagesc(log10(abs(Salpha_ba)))
set(gca,'YDir','normal'); colormap viridis;
set(gca,'XTick',[1,round(length(alpha)/4),round(length(alpha)/2),round(3*length(alpha)/4),length(alpha)],'XTickLabel',[alpha(1),alpha(round(end/4)),alpha(round(end/2)),alpha(round(3*end/4)),alpha(end)])
set(gca,'YTick',[1,round(length(lam)/4),round(length(lam)/2),round(3*length(lam)/4),length(lam)],'YTickLabel',[lam(1),lam(round(end/4)),lam(round(end/2)),lam(round(3*end/4)),lam(end)])

subplot(1,3,2);
imagesc(log10(abs(Salpha_ps)))
set(gca,'YDir','normal');colormap viridis;
set(gca,'XTick',[1,round(length(alpha)/4),round(length(alpha)/2),round(3*length(alpha)/4),length(alpha)],'XTickLabel',[alpha(1),alpha(round(end/4)),alpha(round(end/2)),alpha(round(3*end/4)),alpha(end)])
set(gca,'YTick',[1,round(length(lam)/4),round(length(lam)/2),round(3*length(lam)/4),length(lam)],'YTickLabel',[lam(1),lam(round(end/4)),lam(round(end/2)),lam(round(3*end/4)),lam(end)])

subplot(1,3,3);
imagesc(log10(abs(Salpha_ua)))
set(gca,'YDir','normal');colormap viridis;
set(gca,'XTick',[1,round(length(alpha)/4),round(length(alpha)/2),round(3*length(alpha)/4),length(alpha)],'XTickLabel',[alpha(1),alpha(round(end/4)),alpha(round(end/2)),alpha(round(3*end/4)),alpha(end)])
set(gca,'YTick',[1,round(length(lam)/4),round(length(lam)/2),round(3*length(lam)/4),length(lam)],'YTickLabel',[lam(1),lam(round(end/4)),lam(round(end/2)),lam(round(3*end/4)),lam(end)])

k=0:0.01:1;
for j = 1:length(k)
    Salpha_all(:,:,j) = Salpha_gen(k(j));
    Smu_all(:,:,j) = Smu_gen(k(j));
end

[maxS_for_alpha, maxks_for_alpha] = max(abs(Salpha_all),[],3);

ks_for_alpha =zeros(size(maxks_for_alpha));
for m=1:size(maxks_for_alpha,1)
    for n=1:size(maxks_for_alpha,2)
    ks_for_alpha(m,n) = k(maxks_for_alpha(m,n));
    end
end
    
figure(11); contourf(maxS_for_alpha,21);
title('maximum \alpha sensitivity as f(\alpha, \lambda)');
set(gca,'YDir','normal');colormap viridis;
set(gca,'XTick',[1,round(length(alpha)/4),round(length(alpha)/2),round(3*length(alpha)/4),length(alpha)],'XTickLabel',[alpha(1),alpha(round(end/4)),alpha(round(end/2)),alpha(round(3*end/4)),alpha(end)])
set(gca,'YTick',[1,round(length(lam)/4),round(length(lam)/2),round(3*length(lam)/4),length(lam)],'YTickLabel',[lam(1),lam(round(end/4)),lam(round(end/2)),lam(round(3*end/4)),lam(end)])

figure(12); contourf(ks_for_alpha,21,'LineStyle','none'); caxis([0 1])%caxis([0.25 0.75])
title('maximum k to maximize \alpha sensitivity as f(\alpha, \lambda)');
set(gca,'YDir','normal'); colormap viridis;
set(gca,'XTick',[1,round(length(alpha)/4),round(length(alpha)/2),round(3*length(alpha)/4),length(alpha)],'XTickLabel',[alpha(1),alpha(round(end/4)),alpha(round(end/2)),alpha(round(3*end/4)),alpha(end)])
set(gca,'YTick',[1,round(length(lam)/4),round(length(lam)/2),round(3*length(lam)/4),length(lam)],'YTickLabel',[lam(1),lam(round(end/4)),lam(round(end/2)),lam(round(3*end/4)),lam(end)])


[maxS_for_mu, maxks_for_mu] = max(abs(Smu_all),[],3);

ks_for_mu =zeros(size(maxks_for_mu));
for m=1:size(maxks_for_mu,1)
    for n=1:size(maxks_for_mu,2)
    ks_for_mu(m,n) = k(maxks_for_mu(m,n));
    end
end
    
figure(21); contourf(maxS_for_mu,21);
title('maximum \mu sensitivity as f(\alpha, \lambda)');
set(gca,'YDir','normal');colormap viridis;
set(gca,'XTick',[1,round(length(alpha)/4),round(length(alpha)/2),round(3*length(alpha)/4),length(alpha)],'XTickLabel',[alpha(1),alpha(round(end/4)),alpha(round(end/2)),alpha(round(3*end/4)),alpha(end)])
set(gca,'YTick',[1,round(length(lam)/4),round(length(lam)/2),round(3*length(lam)/4),length(lam)],'YTickLabel',[lam(1),lam(round(end/4)),lam(round(end/2)),lam(round(3*end/4)),lam(end)])

figure(22); contourf(ks_for_mu,21);
title('maximum k to maximize \mu sensitivity as f(\alpha, \lambda)');
set(gca,'YDir','normal');colormap viridis;
set(gca,'XTick',[1,round(length(alpha)/4),round(length(alpha)/2),round(3*length(alpha)/4),length(alpha)],'XTickLabel',[alpha(1),alpha(round(end/4)),alpha(round(end/2)),alpha(round(3*end/4)),alpha(end)])
set(gca,'YTick',[1,round(length(lam)/4),round(length(lam)/2),round(3*length(lam)/4),length(lam)],'YTickLabel',[lam(1),lam(round(end/4)),lam(round(end/2)),lam(round(3*end/4)),lam(end)])

Smu_ua = Smu_gen(1);
Smu_ba = Smu_gen(0);
Smu_ps = Smu_gen(0.5);
Smu_75 = Smu_gen(0.75);
Smu_25 = Smu_gen(0.25);

a=301;
figure; plot(lam,abs(Smu_ba(:,a)),'r'); hold on;  plot(lam,Smu_ps(:,a),'g'); plot(lam,Smu_ua(:,a),'b'); plot(lam,Smu_75(:,a),'cyan'); plot(lam,Smu_25(:,a),'Color',[0.8 0.8 0]) 
xlim([1 1.4])
xlabel('\lambda');
ylabel('S_\mu')
title('S_\mu at a value of \alpha = -2 for biaxial, pure shear, and uniaxial (RGB)');

a=631;
figure; plot(lam,abs(Smu_ba(:,a)),'r'); hold on;  plot(lam,Smu_ps(:,a),'g'); plot(lam,Smu_ua(:,a),'b'); plot(lam,Smu_75(:,a),'cyan'); plot(lam,Smu_25(:,a),'Color',[0.8 0.8 0]) 
xlim([1 1.4])
xlabel('\lambda');
ylabel('S_\mu')
title('S_\mu at a value of \alpha = 1.3 for biaxial, pure shear, and uniaxial (RGB)');

a=1001;
figure; plot(lam,abs(Smu_ba(:,a)),'r'); hold on;  plot(lam,Smu_ps(:,a),'g'); plot(lam,Smu_ua(:,a),'b'); plot(lam,Smu_75(:,a),'cyan'); plot(lam,Smu_25(:,a),'Color',[0.8 0.8 0]) 
xlim([1 1.4])
xlabel('\lambda');
ylabel('S_\mu')
title('S_\mu at a value of \alpha = 5 for biaxial, pure shear, and uniaxial (RGB)');


%%Next thing to do is put in the expression for energy and look at the
%%energy as a function of k, alpha, and lambda
figure; contourf(W_gen(0.5),21);

figure(51); 
% lamidx = 268; % 1.7 for exponential linspace
lamidx = 71;   % 1.7 for regular linspace
%NH
plot(k,squeeze(Salpha_all(lamidx,701,:))); hold on; 
%alpha=5
plot(k,squeeze(Salpha_all(lamidx,1001,:)));
%alpha=1.3
plot(k,squeeze(Salpha_all(lamidx,631,:)));
%alpha=-2
plot(k,-squeeze(Salpha_all(lamidx,301,:)));

figure(52);
%NH
plot(k,squeeze(Salpha_all(lamidx,701,:))/max(squeeze(Salpha_all(lamidx,701,:)))); hold on; 
%alpha=5
plot(k,squeeze(Salpha_all(lamidx,1001,:))/max(squeeze(Salpha_all(lamidx,1001,:))));
%alpha=1.3
plot(k,squeeze(Salpha_all(lamidx,631,:))/max(squeeze(Salpha_all(lamidx,631,:))));
%alpha=-2
plot(k,abs(squeeze(Salpha_all(lamidx,301,:)))/max(abs(squeeze(Salpha_all(lamidx,301,:)))));

