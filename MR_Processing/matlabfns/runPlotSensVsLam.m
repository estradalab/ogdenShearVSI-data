%runPlotSensVsLam

close all;
clear all;

K = 0;
mu = 1E6;
a = 1;
mu01 = mu*(1-a); mu10 = mu*a;
nrm = true;

k = 0.9;

power = (log(3)-log(2))/(log(2));

L = 1.01:0.2:5;
%zr = zeros(length(L),3,3);
zr = zeros(length(L),1);
sensToa_uni = zr; sensToa_ps = zr; sensToa_bi = zr; sensToa_arb = zr;
sensTomu_uni = zr; sensTomu_ps = zr; sensTomu_bi = zr; sensTomu_arb = zr;
sensTomu10_uni = zr; sensTomu10_ps = zr; sensTomu10_bi = zr; sensTomu10_arb = zr;
sensTomu01_uni = zr; sensTomu01_ps = zr; sensTomu01_bi = zr; sensTomu01_arb = zr;

for idx = 1:length(L)
    
    Lam = L(idx);
    F_uni = [Lam 0 0; 0 1/sqrt(Lam) 0; 0 0 1/sqrt(Lam)];
    dF_uni = [1 0 0; 0 -Lam^(-3/2)/2 0; 0 0 -Lam^(-3/2)/2];
    

    
    %%Normalization: Make all deformations contribute the same energy?
    W_MR_uni = W_MR_alt(F_uni,mu,a,0);
    if nrm
        syms x;
        vals = solve(W_MR_alt([x 0 0; 0 1 0; 0 0 1/x],mu,a,0) == W_MR_uni, x>0);
        x = eval(vals);
        if abs(x(1)-Lam)<abs(x(2)-Lam)
            lam_ps(idx) = x(1);
        else
            lam_ps(idx) = x(2);
        end
        
        syms x;
        vals = vpasolve(W_MR_alt([(x)^0.5 0 0; 0 (x)^0.5 0; 0 0 1/x],mu,a,0) == W_MR_uni);
        x = abs(eval(vals));
        lam_bi(idx) = (x);
        
        syms x
        vals = vpasolve(W_MR_alt([x^(fofk(k) + 0.5*fofk(1-k)) 0 0; 0 (x)^(-0.5*fofk(k) + 0.5*fofk(1-k)) 0; 0 0 (x)^(-0.5*fofk(k) - fofk(1-k))],mu,a,0) == W_MR_uni);
        x = abs(eval(vals));
        lam_arb(idx) = x;%^(fofk(k) + 0.5*fofk(1-k));
        %     if abs(x(1)-lam)<abs(x(2)-lam)
        %         lam_bi(idx) = x(1);
        %     else
        %         lam_bi(idx) = x(2);
        %     end
    else
        lam_ps(idx) = Lam;
        lam_bi(idx) = Lam;
        lam_arb(idx) = Lam;
    end
    
    F_ps = [lam_ps(idx) 0 0; 0 1 0; 0 0 1/lam_ps(idx)];
    dF_ps = [1 0 0; 0 0 0; 0 0 -1/lam_ps(idx)^2];
    W_MR_alt(F_ps,mu,a,0);
    F_bi = [sqrt(lam_bi(idx)) 0 0; 0 sqrt(lam_bi(idx)) 0; 0 0 1/lam_bi(idx)];
    dF_bi = [1/2*lam_bi(idx)^(-1/2) 0 0; 0 1/2*lam_bi(idx)^(-1/2) 0; 0 0 -lam_bi(idx)^-2];
    W_MR_alt(F_bi,mu,a,0);
%    F_arb = [lam_arb(idx)^(0.5*(1+k)) 0 0; 0 lam_arb(idx)^(0.5*(1-2*k)) 0; 0 0 lam_arb(idx)^(0.5*(-2+k))];
%    F_arb = [lam_arb(idx)^(0.5*(1+k)) 0 0; 0 lam_arb(idx)^(0.5*(1-2*k)) 0; 0 0 lam_arb(idx)^(0.5*(-2+k))];
    F_arb = [lam_arb(idx)^(fofk(k) + 0.5*fofk(1-k)) 0 0; 0 lam_arb(idx)^(-0.5*fofk(k) + 0.5*fofk(1-k)) 0; 0 0 lam_arb(idx)^(-0.5*fofk(k) - fofk(1-k))];
    %dF_arb = lam_arb(idx)^((-1-k)/2)/(0.5*(1+k))*0.5*[1+k 0 0; 0 1-2*k 0; 0 0 -2+k].*F_arb;
    %dF_arb = (lam_arb(idx)^(k^2-3/2*k-1/2))/(-k^2+3/2*k+1/2)*[(fofk(k) + 0.5*fofk(1-k)) 0 0; 0 (-0.5*fofk(k) + 0.5*fofk(1-k)) 0; 0 0 (-0.5*fofk(k) - fofk(1-k))].*F_arb;
    dF_arb = (1/lam_arb(idx))*[(fofk(k) + 0.5*fofk(1-k)) 0 0; 0 (-0.5*fofk(k) + 0.5*fofk(1-k)) 0; 0 0 (-0.5*fofk(k) - fofk(1-k))].*F_arb;
    W_MR_alt(F_arb,mu,a,0);
    
    %lam = Lam^(0.5*(1+k));
    %F_arb = [lam^(0.5*(1+k)) 0 0; 0 lam^(0.5*(1-2*k)) 0; 0 0 lam^(0.5*(-2+k))];
    %hold on; plot(I1(F_arb), I2(F_arb), '.', 'Color','green')
    
    % %     deltaF_uni_ = (0.00001*(F_uni-eye(3,3))+eye(3,3));
    % %     deltaF_ps_ = (0.00001*(F_ps-eye(3,3))+eye(3,3));
    % %     deltaF_bi_ = (0.00001*(F_bi-eye(3,3))+eye(3,3));
    %
    %   %%The Alternate form of the Mooney-Rivlin model (mu*alpha and mu*(1-alpha))
    %   %%Analytical dependence of sensitivity of a on lambda
    % %     %Calculate 4th rank d2Pi/dadF
    % %     d2PidadF_uni_ = d2PidadF_MR(F_uni, mu, a, K);
    % %     d2PidadF_ps_ = d2PidadF_MR(F_ps, mu, a, K);
    % %     d2PidadF_bi_ = d2PidadF_MR(F_bi, mu, a, K);
    %
    %     %Calculate the stress first derivatives also w.r.t a
    %     dPida_uni_ = dPida_MR_alt(F_uni, mu, a, K);
    %     dPida_ps_ = dPida_MR_alt(F_ps, mu, a, K);
    %     dPida_bi_ = dPida_MR_alt(F_bi, mu, a, K);
    %
    %     %Multiply that by 2nd rank deltaF to get sensitivity
    % %     sensToa_uni_ = zeros(3,3);
    % %     for n=1:3
    % %         for m=1:3
    % %             sensToa_uni_ = sensToa_uni_+d2PidadF_uni_(:,:,n,m)*dF_uni(n,m);%deltaF_uni_(n,m);%*dij(n,m));%
    % %         end
    % %     end
    %     sensToa_uni(idx) = trace(dPida_uni_.*dF_uni);
    % %     sensToa_ss_ = zeros(3,3);
    % %     for n=1:3
    % %         for m=1:3
    % %             sensToa_ss_ = sensToa_ss_+d2PidadF_ps_(:,:,n,m)*dF_ps(n,m);%deltaF_ps_(n,m);%*dij(n,m));%
    % %         end
    % %     end
    %     sensToa_ps(idx) = trace(dPida_ps_.*dF_ps);
    % %     sensToa_bi_ = zeros(3,3);
    % %     for n=1:3
    % %         for m=1:3
    % %             sensToa_bi_ = sensToa_bi_+d2PidadF_bi_(:,:,n,m)*dF_bi(n,m);%deltaF_bi_(n,m);%*dij(n,m));%
    % %         end
    % %     end
    %     sensToa_bi(idx) = trace(dPida_bi_.*dF_bi);
    %     %sensToa_uni(idx,:,:) = sensToa_uni_;
    %     %sensToa_ss(idx,:,:) = sensToa_ss_;
    %     %sensToa_bi(idx,:,:) = sensToa_bi_;
    %
    % %    %%Analytical dependence of sensitivity of mu on lambda
    % %     %Calculate 4th rank d2Pi/dadF
    % %     d2PidmudF_uni_ = d2PidmudF_MR(F_uni, mu, a, K);
    % %     d2PidmudF_ps_ = d2PidmudF_MR(F_ps, mu, a, K);
    % %     d2PidmudF_psNH_ = d2PidMudF_NH(F_ps, mu, K);
    % %     d2PidmudF_bi_ = d2PidmudF_MR(F_bi, mu, a, K);
    %     %Calculate the stress first derivatives also w.r.t a
    %     dPidmu_uni_ = dPidMu_MR_alt(F_uni, mu, a, K);
    %     %dPidmu_uniNH = dPidMu_NH(F_uni,mu,K);
    %     dPidmu_ps_ = dPidMu_MR_alt(F_ps, mu, a, K);
    %     dPidmu_bi_ = dPidMu_MR_alt(F_bi, mu, a, K);
    %     %Multiply that by 2nd rank deltaF to get sensitivity
    % %     sensTomu_uni_ = zeros(3,3);
    % %     for n=1:3
    % %         for m=1:3
    % %             sensTomu_uni_ = sensTomu_uni_+d2PidmudF_uni_(:,:,n,m)*dF_uni(n,m);%deltaF_uni_(n,m);%*dij(n,m));%
    % %         end
    % %     end
    % %     sensTomu_ss_ = zeros(3,3);
    % %     for n=1:3
    % %         for m=1:3
    % %             sensTomu_ss_ = sensTomu_ss_+d2PidmudF_ps_(:,:,n,m)*dF_ps(n,m);%*deltaF_ps_(n,m);%*dij(n,m));%
    % %         end
    % %     end
    % %     sensTomu_bi_ = zeros(3,3);
    % %     for n=1:3
    % %         for m=1:3
    % %             sensTomu_bi_ = sensTomu_bi_+d2PidmudF_bi_(:,:,n,m)*dF_bi(n,m);%*deltaF_bi_(n,m);%*dij(n,m));%
    % %         end
    % %     end
    % %     sensTomu_uni(idx,:,:) = sensTomu_uni_;
    % %     sensTomu_ss(idx,:,:) = sensTomu_ss_;
    % %     sensTomu_bi(idx,:,:) = sensTomu_bi_;
    % sensTomu_uni(idx) = trace(dPidmu_uni_.*dF_uni);
    % sensTomu_ps(idx) = trace(dPidmu_ps_.*dF_ps);
    % sensTomu_bi(idx) = trace(dPidmu_bi_.*dF_bi);

    dPida_uni_ = dPida_MR_alt(F_uni, mu, a, K);
    dPida_ps_ = dPida_MR_alt(F_ps, mu, a, K);
    dPida_bi_ = dPida_MR_alt(F_bi, mu, a, K);
    dPida_arb_ = dPida_MR_alt(F_arb, mu, a, K);
    
    dPidmu_uni_ = dPidMu_MR_alt(F_uni, mu, a, K);
    dPidmu_ps_ = dPidMu_MR_alt(F_ps, mu, a, K);
    dPidmu_bi_ = dPidMu_MR_alt(F_bi, mu, a, K);
    dPidmu_arb_ = dPidMu_MR_alt(F_arb, mu, a, K);
    
    dPidmu01_uni_ = dPidmu01_MR(F_uni, mu10, mu01, K);
    dPidmu01_ps_ = dPidmu01_MR(F_ps, mu10, mu01, K);
    dPidmu01_bi_ = dPidmu01_MR(F_bi, mu10, mu01, K);
    dPidmu01_arb_ = dPidmu01_MR(F_arb, mu10, mu01, K);
    
    dPidmu10_uni_ = dPidmu10_MR(F_uni, mu10, mu01, K);
    dPidmu10_ps_ = dPidmu10_MR(F_ps, mu10, mu01, K);
    dPidmu10_bi_ = dPidmu10_MR(F_bi, mu10, mu01, K);
    dPidmu10_arb_ = dPidmu10_MR(F_arb, mu10, mu01, K);
    
    % Second Order version
%     d2PidadF_uni_ = d2PidadF_MR(F_uni, mu, a, K);
%     d2PidadF_ps_ = d2PidadF_MR(F_ps, mu, a, K);
%     d2PidadF_bi_ = d2PidadF_MR(F_bi, mu, a, K);
%     
%     d2PidmudF_uni_ = d2PidmudF_MR(F_uni, mu, a, K);
%     d2PidmudF_ps_ = d2PidmudF_MR(F_ps, mu, a, K);
%     d2PidmudF_bi_ = d2PidmudF_MR(F_bi, mu, a, K);
%     
%     d2Pidmu01dF_uni_ = d2Pidmu01dF_MR(F_uni, mu10, mu01, K);
%     d2Pidmu01dF_ps_ = d2Pidmu01dF_MR(F_ps, mu10, mu01, K);
%     d2Pidmu01dF_bi_ = d2Pidmu01dF_MR(F_bi, mu10, mu01, K);
%     
%     d2Pidmu10dF_uni_ = d2PidMudF_NH(F_uni, mu10, K);
%     d2Pidmu10dF_ps_ = d2PidMudF_NH(F_ps, mu10, K);
%     d2Pidmu10dF_bi_ = d2PidMudF_NH(F_bi, mu10, K);
%     
    firstO = true;
    
    if firstO
        
        sensToa_uni(idx) = trace(dPida_uni_.*dF_uni);
        sensToa_ps(idx) = trace(dPida_ps_.*dF_ps);
        sensToa_bi(idx) = trace(dPida_bi_.*dF_bi);
        sensToa_arb(idx) = trace(dPida_arb_.*dF_arb);
        
        sensTomu_uni(idx) = trace(dPidmu_uni_.*dF_uni);
        sensTomu_ps(idx) = trace(dPidmu_ps_.*dF_ps);
        sensTomu_bi(idx) = trace(dPidmu_bi_.*dF_bi);
        sensTomu_arb(idx) = trace(dPidmu_arb_.*dF_arb);
        
        sensTomu01_uni(idx) = trace(dPidmu01_uni_.*dF_uni);
        sensTomu01_ps(idx) = trace(dPidmu01_ps_.*dF_ps);
        sensTomu01_bi(idx) = trace(dPidmu01_bi_.*dF_bi);
        sensTomu01_arb(idx) = trace(dPidmu01_arb_.*dF_arb);
        
        sensTomu10_uni(idx) = trace(dPidmu10_uni_.*dF_uni);
        sensTomu10_ps(idx) = trace(dPidmu10_ps_.*dF_ps);
        sensTomu10_bi(idx) = trace(dPidmu10_bi_.*dF_bi);
        sensTomu10_arb(idx) = trace(dPidmu10_arb_.*dF_arb);
    else
        sensToa_uni_ = 0; sensToa_ps_ = 0; sensToa_bi_ = 0;
        sensTomu_uni_ = 0; sensTomu_ps_ = 0; sensTomu_bi_ = 0;
        for k=3:-1:1
            for m=3:-1:1
                sensToa_uni_ = sensToa_uni_+d2PidadF_uni_(:,:,k,m)*dF_uni(k,m);
                sensToa_ps_ = sensToa_ps_+d2PidadF_ps_(:,:,k,m)*dF_ps(k,m);
                sensToa_bi_ = sensToa_bi_+d2PidadF_bi_(:,:,k,m)*dF_bi(k,m);
                
                sensTomu_uni_ = sensTomu_uni_+d2PidmudF_uni_(:,:,k,m)*dF_uni(k,m);
                sensTomu_ps_ = sensTomu_ps_+d2PidmudF_ps_(:,:,k,m)*dF_ps(k,m);
                sensTomu_bi_ = sensTomu_bi_+d2PidmudF_bi_(:,:,k,m)*dF_bi(k,m);                
                
            end
        end
        
    end
    
    
    W(idx) = W_MR_uni;
    I1_lambi(idx) = I1(F_bi);
    I2_lambi(idx) = I2(F_bi);
    I1_lamps(idx) = I1(F_ps);
    I2_lamps(idx) = I2(F_ps);
    I1_lamuni(idx) = I1(F_uni);
    I2_lamuni(idx) = I2(F_uni);
    I1_lamarb(idx) = I1(F_arb);
    I2_lamarb(idx) = I2(F_arb);    
end

%%
%close all;
figure(1);%(a*1000+1);
% ind = 0;
% for i=1:3
%     for j=1:3
%         ind = ind+1;
%        subplot(3,3,ind);
%        plot(L,sensToa_uni(:,i,j)); hold on;
%        plot(L,sensToa_ss(:,i,j));
%        plot(L,sensToa_bi(:,i,j));
%     end
% end
subplot(1,2,1);
plot(L,abs(sensToa_uni)); hold on;
plot(lam_ps,abs(sensToa_ps)); hold on;
plot(lam_bi,abs(sensToa_bi)); hold on;
plot(lam_arb,abs(sensToa_arb)); hold on;
title('Stress sensitivity to \alpha');

subplot(1,2,2);
plot(L,sensTomu_uni); hold on;
plot(lam_ps,sensTomu_ps); hold on;
plot(lam_bi,sensTomu_bi); hold on;
plot(lam_arb,sensTomu_arb); hold on;
title('Stress sensitivity to \mu');


figure(2)%(a*1000+2);
subplot(1,2,1);
plot(L,abs(sensTomu10_uni)); hold on;
plot(lam_ps,abs(sensTomu10_ps)); hold on;
plot(lam_bi,abs(sensTomu10_bi)); hold on;
plot(lam_arb,abs(sensTomu10_arb)); hold on;
title('Stress sensitivity to \mu_{10}');

subplot(1,2,2);
plot(L,sensTomu01_uni); hold on;
plot(lam_ps,sensTomu01_ps); hold on;
plot(lam_bi,sensTomu01_bi); hold on;
plot(lam_arb,sensTomu01_arb); hold on;
title('Stress sensitivity to \mu_{01}');
% ind = 0;
% for i=1:3
%     for j=1:3
%         ind = ind+1;
%        subplot(3,3,ind);
%        plot(L,sensTomu_uni(:,i,j)); hold on;
%        plot(L,sensTomu_ss(:,i,j));
%        plot(L,sensTomu_bi(:,i,j));
%     end
% end


% figure;
% subplot(1,3,1);
% plot(L,sensTomu_uni(:,1,1),'Color','b'); hold on;
% plot(L,sensTomu_uni(:,2,2),'Color','b');
% plot(L,sensTomu_uni(:,3,3),'Color','b');
% subplot(1,3,2);
% plot(L,sensTomu_ss(:,1,1),'Color',[1 .5 0]); hold on;
% plot(L,sensTomu_ss(:,2,2),'Color',[1 .5 0]);
% plot(L,sensTomu_ss(:,3,3),'Color',[1 .5 0]);
% subplot(1,3,3);
% plot(L,sensTomu_bi(:,1,1),'Color','y'); hold on;
% plot(L,sensTomu_bi(:,2,2),'Color','y');
% plot(L,sensTomu_bi(:,3,3),'Color','y');



 figure (a*1000+3);
 plot(I1_lamuni,I2_lamuni,'.','MarkerSize',12); hold on;
 plot(I1_lamps,I2_lamps,'.','MarkerSize',12);
 plot(I1_lambi,I2_lambi,'.','MarkerSize',12);
 plot(I1_lamarb,I2_lamarb,'.','MarkerSize',12);
 title('Plot of I_2-bar vs. I_1-bar');
 
% for i=0:10
%    k = i/10;
%    F_arb_ = [lam_arb(idx)^(fofk(k) + 0.5*fofk(1-k)) 0 0; 0 lam_arb(idx)^(-0.5*fofk(k) + 0.5*fofk(1-k)) 0; 0 0 lam_arb(idx)^(-0.5*fofk(k) - fofk(1-k))];
%    hold on; plot(I1(F_arb_), I2(F_arb_), '.', 'Color','green','MarkerSize',12)
% end
