clear all; close all;

%load('NeoHook_pars.mat')
hdir = pwd;

sample = '10A-glossy';

switch sample
    case '10A-sanded'
        comm = '-10AcorrFijHamMRcorrNO-repVFs';
        matcomm = 'fm-newpad-ham-MRcorr-1last';
        fixedJ = 1.002065; %Offset value for 10A-Sanded (calculated as mean of all central J)
        date = '191204';
        nums = {'1821','1548','1805','1639','1454','1954','1859','2049','2146','2304'};
        forcefn = '\ADMET\disploads.mat';
    case '2A'
        comm = '-2AhamMRLoadJcorrNO-repVFs';
        %comm = '-2ANoFcorrection';
        matcomm = 'fm-newpad-ham-MRcorr-1last';%
        fixedJ = 1.002742; %Offset value for 2A-Rect1to5 (calculated as mean of all central J)
        date = '191209';
        nums={'1820','1617','1416','1718','1517'};
        forcefn = '\ADMET\disploads.mat';
    case '10A-glossy'
        comm = '-10AcorrFijHamMRcorrNO-repVFs';
        %comm = '-2ANoFcorrection';
        matcomm = 'fm-newpad-ham-MRcorr-1last';%
        fixedJ = 1.002451; %Offset value for 2A-Rect1to5 (calculated as mean of all central J)
        date = '191209';        
        nums={'1605','1718','1813'};
        forcefn = '\ADMET\disploads-10ANoFcorrection.mat';
end



type = 'normVFnormPhirange';
blurfilt = '000';
constit = 'GenNeoHook1NoJ';
load([hdir '\mat-vfm\' type '\' constit '-fmin-filt' blurfilt '-' date comm '.mat']);
load([hdir  forcefn]);

offset = false;
shearonly = false;
%options: 'true','fixedJ','false'
isIncompressible = 'fixedJ'; 


checkLamFij = 'Lam';

midpct = [0.5 0.8 0.8]; %[0.50 0.75 0.75];

%%

%middle = [50:90];%50:90

for i=1:length(nums)
    fn = glob([hdir '\mat\DESTE_strains_' nums{i} '*']);
    load(fn{1},'wiggle','origpars');
    X1 = origpars.axis1;
    selstr(i) = wiggle(1);
    
    Ffn = glob([hdir '\mat\Fij_optimal3h' blurfilt '_' nums{i} matcomm '.mat']);
    load(Ffn{1},'Fij','mask');
            
    [detF{i},Fbarij] = calculateJ(Fij);
    switch isIncompressible
        case 'true'
            detF{i} = calculateJ(Fbarij);
            tFij{1} = Fbarij;
            Finv_ij = calculateInvFij(Fbarij);
            Fij = tFij{1};
        case 'false'
            tFij{1} = Fij;
            Finv_ij = calculateInvFij(Fij);
            Fij = tFij{1};
        case 'fixedJ'
            J_ = fixedJ; J_div = J_^(1/3);
            tFij{1} = Fij;
            for iii = 1:3
                tFij{1}{iii,iii} = Fij{iii,iii}./J_div; 
            end
            Finv_ij = calculateInvFij(tFij{1});
            Fij = tFij{1}; 
            [detF{i}] = calculateJ(Fij);
    end        
    
    linmask = nansum(nansum(mask,3),2);
    inds = find(~isnan(linmask)&linmask>20);
    rnginds = max(inds)-min(inds);
    midinds = (max(inds)+min(inds))/2;
    middleX = round(midinds-rnginds*midpct(1)/2):round(midinds+rnginds*midpct(1)/2);
    
    linmask = nansum(nansum(mask,3),1);
    inds = find(~isnan(linmask)&linmask>40);
    rnginds = max(inds)-min(inds);
    midinds = (max(inds)+min(inds))/2;
    middleY = round(midinds-rnginds*midpct(2)/2):round(midinds+rnginds*midpct(2)/2);
    
    linmask = nansum(nansum(mask,2),1);
    inds = find(~isnan(linmask)&linmask>40);
    rnginds = max(inds)-min(inds);
    midinds = (max(inds)+min(inds))/2;
    middleZ = round(midinds-rnginds*midpct(3)/2):round(midinds+rnginds*midpct(3)/2);
    
    for p=3:-1:1
        for q=3:-1:1
        tmp = (Fij{p,q}(middleX,middleY,middleZ));
        F(p,q,i) = nanmedian(tmp(:));
        end
    end

end

errF = F(:,:,1)-eye(3,3);

if offset
    if shearonly
        errF(1,1) = 0; errF(2,2) = 0; errF(3,3) = 0;
    end
    F = F-repmat(errF,[1 1 length(nums)]);
end

for i=1:length(nums)
    fn = glob([hdir '\mat\DESTE_strains_' nums{i} '*']);
    load(fn{1},'wiggle');
    selstr(i) = wiggle(1);
    
    Ffn = glob([hdir '\mat\Fij_optimal3h' blurfilt '_' nums{i} matcomm '.mat']);
    load(Ffn{1},'Fij','mask');
    
    if offset
        for p=1:3
            for q=1:3
                Fij{p,q} = Fij{p,q}-errF(p,q);
            end
        end
    end
    
    [detF{i},Fbarij] = calculateJ(Fij);
    switch isIncompressible
        case 'true'
            detF{i} = calculateJ(Fbarij);
            tFij{1} = Fbarij;
            Finv_ij = calculateInvFij(Fbarij);
            Fij = tFij{1};
        case 'false'
            tFij{1} = Fij;
            Finv_ij = calculateInvFij(Fij);
            Fij = tFij{1};
        case 'fixedJ'
            J_ = fixedJ; J_div = J_^(1/3);
            tFij{1} = Fij;
            for iii = 1:3
                tFij{1}{iii,iii} = Fij{iii,iii}./J_div; 
            end
            Finv_ij = calculateInvFij(tFij{1});
            Fij = tFij{1}; 
            [detF{i}] = calculateJ(Fij);
    end 
    
    lilF11{i} = Fij{1,1}(middleX,middleY,middleZ);
    F11(i) = nanmedian(lilF11{i}(:));
    stdF11(i) = nanstd(lilF11{i}(:));
    
    lilF12{i} = Fij{1,2}(middleX,middleY,middleZ);
    F12(i) = nanmedian(lilF12{i}(:));
    stdF12(i) = nanstd(lilF12{i}(:));
    
    lilF13{i} = Fij{1,3}(middleX,middleY,middleZ);
    F13(i) = nanmedian(lilF13{i}(:));
    stdF13(i) = nanstd(lilF13{i}(:));
    
    lilF21{i} = Fij{2,1}(middleX,middleY,middleZ);
    F21(i) = nanmedian(lilF21{i}(:));
    stdF21(i) = nanstd(lilF21{i}(:));
    
    lilF22{i} = Fij{2,2}(middleX,middleY,middleZ);
    F22(i) = nanmedian(lilF22{i}(:));
    stdF22(i) = nanstd(lilF22{i}(:));
    
    lilF23{i} = Fij{2,3}(middleX,middleY,middleZ);
    F23(i) = nanmedian(lilF23{i}(:));
    stdF23(i) = nanstd(lilF23{i}(:));
    
    lilF31{i} = Fij{3,1}(middleX,middleY,middleZ);
    F31(i) = nanmedian(lilF31{i}(:));
    stdF31(i) = nanstd(lilF31{i}(:));
    
    lilF32{i} = Fij{3,2}(middleX,middleY,middleZ);
    F32(i) = nanmedian(lilF32{i}(:));
    stdF32(i) = nanstd(lilF32{i}(:));
    
    lilF33{i} = Fij{3,3}(middleX,middleY,middleZ);
    F33(i) = nanmedian(lilF33{i}(:));
    stdF33(i) = nanstd(lilF33{i}(:));
    
    tUij = calculateUij(tFij);
    Uij{i} = tUij{1};
    
    tLami = calculateLami(tUij);
    Lami{i} = tLami{1};
    
    detF_(i) = det(squeeze(F(:,:,i)));
    
    stdF_(:,:,i) = [stdF11(i) stdF12(i) stdF13(i);...
        stdF21(i) stdF22(i) stdF23(i);...
        stdF31(i) stdF32(i) stdF33(i)];
    
    lilF11inv{i} = Finv_ij{1,1}(middleX,middleY,middleZ);
    invF11(i) = nanmean(lilF11inv{i}(:));
    
    lilLam1{i} = Lami{i}{1}(middleX,middleY,middleZ);
    Lam1(i) = nanmedian(lilLam1{i}(:));
    stdLam1(i) = nanstd(lilLam1{i}(:));
    lilLam2{i} = Lami{i}{2}(middleX,middleY,middleZ);
    Lam2(i) = nanmedian(lilLam2{i}(:));
    stdLam2(i) = nanstd(lilLam2{i}(:));
    lilLam3{i} = Lami{i}{3}(middleX,middleY,middleZ);
    Lam3(i) = nanmedian(lilLam3{i}(:));
    stdLam3(i) = nanstd(lilLam3{i}(:));
    
    lilJ{i} = detF{i}(middleX,middleY,middleZ);
    volratio_medianpxwise(i) = nanmedian(lilJ{i}(:));
    volratio_meanpxwise(i) = nanmean(lilJ{i}(:));
    stdvolratio(i) = nanstd(lilJ{i}(:));
    
    
end
%selstr = [0 0.5 1 2 3 4 4.5 5 6 7];


L=42/1e3; % mm to m
w=7.4/1e3;
t=8./1e3;
%262 px average
lilmask = mask(middleX,:,:); numpxInXsec = nansum(lilmask(:))/length(middleX);
A=(numpxInXsec/4)/1e6;%24.6/1e6;%w*t; %from mean number of cross-section px on x=101:150

A*1e6/8;
mask3 = nansum(mask,3); mean(mask3(mask3>0))*0.5
mask2 = nansum(mask,2); mean(mask2(mask2>0))*0.5
mask1 = nansum(mask,1); mean(mask1(mask1>0))*0.4375
%correction factor for MRI HIRES - known dim of 8 is scaled to 8.5
%Note: only need to correct for this afterward on the load/area nominal
%stress because the du/dX scaling cancels out for F.
A = A;%*(8/8.5)*(7.4/8);

Pos= smooth(loaddisp.Positionmm);
Loads = smooth(loaddisp.LoadN);%-0.01889; %%%%%Temp offset for 2A

force = interpn(Pos,Loads,selstr);

delta=selstr*1e-3;
exp_stretch=(L+delta)/L;
%exp_stretch = squeeze(F_(1,1,:));
%exp_stretch2 = squeeze(F_(2,2,:));
%exp_stretch3 = squeeze(F_(3,3,:));
exp_stretch = Lam1;
exp_stretch2 = Lam2;
exp_stretch3 = Lam3;
exp_J = Lam1.*Lam2.*Lam3;

%force=[0 0.15 0.28 0.41 0.52 0.63 0.725];%[0 0.33 0.64 0.92 1.17 1.4 1.61]; % N
exp_nom_stress=force/A;%-force(1)/A;

%par(1) = par(1)/1.13;

par = fminpars;%/1.2;%/1.3;
%Override
par(2) = 1e11;
%par = [33496.5439157828,1374041.97501252];
%par = [2.665537357034300e+04,1E9];
%par(2) = 1.385E6; %override bulk for plot
%par = [37.2e+03,2.05e+06];
%par(2) = 1E6;

%Uniaxial model
stretch=1:0.002:max(exp_stretch);
for i=1:length(stretch)
    syms Js
    assume(Js>=1)
    Fs=[stretch(i) 0 0; 0 sqrt(Js/stretch(i)) 0; 0 0 sqrt(Js/stretch(i))];
    Bs=Fs*Fs';
    I1s=trace(Bs);
    model_cauchy=2*par(1)*Js^(-5/3)*(Bs-1/3*I1s*eye(3))+par(2)*(Js-1)*eye(3);
    J=double(vpasolve(model_cauchy(2,2)==0,Js));
    %     J=1;
    F=[stretch(i) 0 0; 0 sqrt(J/stretch(i)) 0; 0 0 sqrt(J/stretch(i))];
    B=F*F';
    I1=trace(B);
    model_true=2*par(1)*J^(-5/3)*(B-1/3*I1*eye(3))+par(2)*(J-1)*eye(3);
    model_nominal=J*F^(-1)*model_true;
    model_S11(i)=model_nominal(1,1);
    trans_stretch(i) = sqrt(J/stretch(i));
    
    tr_sig(i) = trace(model_true);
    model_detF(i) = J;
end

%plot(model_S11,stretch,'--','LineWidth',2,'Color','red')

%whole mask J
%medJ = [1.00015779095399,1.00057238282810,1.00110953434783,1.00161555422045,1.00211083577762,1.00236027484481,1.00290870767779,1.00336023631818,1.00394088013618];
%medJ = [ 1.0004    1.0010    1.0019    1.0028    1.0036    1.0039    1.0046    1.0055    1.0062];


for j=1:length(nums)
    thirdtrsig(j) = (1/3)*(exp_nom_stress(j)/(volratio_medianpxwise(j)*invF11(j)));
end

%%

magmaVals = colormap(magma);

color = [1 0 0];
switch blurfilt(3)
    case '0'
        color = magmaVals(1,:);
    case '1'
        color = magmaVals(round(64*2/9),:);
    case '2'
        color = magmaVals(round(64*4/9),:);
    case '3'
        color = magmaVals(round(64*5/9),:);
    case '4'
        color = magmaVals(round(64*6/9),:);
    case '8'
        color = magmaVals(round(64*8/9),:);
end

%%
figure(100+10*1+str2double(blurfilt(3)))
% subplot(1,3,1)
% plot(exp_stretch,exp_nom_stress,'.','MarkerSize',15,'Color',color)
% hold on
% plot(stretch,model_S11,'--','LineWidth',2,'Color',color)
% xlabel('Stretch')
% ylabel('Nominal stress')
% set(gca,'FontName','Arial','FontSize',14)
% legend('Experimental data','Compressible NeoHookean UT simulation')
% title('Axial nominal stress-strain');
plot(exp_nom_stress,exp_stretch,'.','MarkerSize',15,'Color',color)
hold on;
Y = [];
for i=1:length(lilLam1)
   switch checkLamFij
       case 'Fij'
        Y = [Y, lilF11{i}(:)];
       case 'Lam'
        Y = [Y, lilLam1{i}(:)];   
   end
end
X = round(exp_nom_stress);
iosr.statistics.boxPlot(X,Y,'showViolin',true,'boxColor','none','boxWidth',300,'symbolColor',color,'LineColor',color);
hold on;
plot(model_S11,stretch,'--','LineWidth',2,'Color','k')
xlim([0 ceil(max(model_S11)/1E4)*1e4]);
ylabel('Stretch')
xlabel('Nominal stress')
set(gca,'FontName','Arial','FontSize',12)
%legend('Experimental data','Compressible NeoHookean UT simulation')
title('Axial nominal stress-strain');
set(gca,'View',[90 -90])
lft = fittype({'2*(x-1/x^2)'})
fit(exp_stretch(2:end)',exp_nom_stress(2:end)',lft)

figure(100+10*2+str2double(blurfilt(3)))
% subplot(1,3,2)
plot(exp_nom_stress,exp_stretch2,'.','MarkerSize',15,'Color','r'); hold on
hold on;
Y = [];
for i=1:length(lilLam1)
   switch checkLamFij
       case 'Fij'
        Y = [Y, lilF22{i}(:)];
       case 'Lam'
        Y = [Y, lilLam2{i}(:)];   
   end
end
X = round(exp_nom_stress);
iosr.statistics.boxPlot(X,Y,'showViolin',true,'boxColor','none','boxWidth',300,'symbolColor','red','LineColor','red');
hold on;
plot(model_S11,trans_stretch,'--','LineWidth',2,'Color','k')
plot(exp_nom_stress,exp_stretch3,'.','MarkerSize',15,'Color','b')
Y = [];
for i=1:length(lilLam1)
   switch checkLamFij
       case 'Fij'
        Y = [Y, lilF33{i}(:)];
       case 'Lam'
        Y = [Y, lilLam3{i}(:)];   
   end
end
iosr.statistics.boxPlot(X,Y,'showViolin',true,'boxColor','none','boxWidth',300,'symbolColor','blue','LineColor','blue');
hold on;
xlim([0 ceil(max(model_S11)/1E4)*1e4]);
ylabel('Transverse Stretch')
xlabel('Nominal stress')
set(gca,'FontName','Arial','FontSize',12)
%legend('Experimental data','Compressible NeoHookean UT simulation')
title('Transverse nominal stress-strain');
set(gca,'View',[90 -90])


%plot(exp_stretch2-(exp_stretch2(1)-1),exp_nom_stress,'.','MarkerSize',15,'Color','r'); hold on
%plot(exp_stretch3-(exp_stretch3(1)-1),exp_nom_stress,'.','MarkerSize',15,'Color','b')
% 
% plot(trans_stretch,model_S11,'--','LineWidth',2,'Color',color)
% xlabel('Stretch')
% ylabel('Nominal stress')
% set(gca,'FontName','Arial','FontSize',14)
% legend('Experimental data','Compressible NeoHookean UT simulation')

figure(100+10*3+str2double(blurfilt(3)))
% subplot(1,3,3)
Y = [];
for i=1:length(lilLam1)
Y = [Y, lilJ{i}(:)-1];
end
X = round(thirdtrsig);
iosr.statistics.boxPlot(X,Y,'showViolin',true,'boxColor','none','boxWidth',300,'symbolColor','k','LineColor','k');
hold on;
plot(thirdtrsig,volratio_medianpxwise-1,'.','MarkerSize',15,'Color','k')
plot(tr_sig/3,model_detF-1,'--','LineWidth',2,'Color','k')
hold on;
%plot(volratio_meanpxwise-1,thirdtrsig,'.','MarkerSize',15,'Color','k')
plot(thirdtrsig,exp_J-1,'.','MarkerSize',15,'Color','green');
%plot((exp_stretch3-(exp_stretch3(1)-1)).*(exp_stretch2-(exp_stretch2(1)-1)).*(exp_stretch-(exp_stretch(1)-1))-1,thirdtrsig,'.','MarkerSize',15,'Color','green');
xlim([0 ceil(max(thirdtrsig)/1E3)*1e3]);
xlabel('1/3*trace(Cauchy Stress)')
ylabel('J-1')
set(gca,'FontName','Arial','FontSize',12)
%legend('Experimental data','Compressible NeoHookean UT simulation')
title('Volumetric stress-strain');
set(gca,'View',[90 -90])

%%

figure(200+str2double(blurfilt(3)))
subplot(1,3,1)
histogram(lilLam1{end}); hold on;
histogram(lilF11{end})
subplot(1,3,2)
histogram(lilLam2{end}); hold on;
histogram(lilF22{end})
subplot(1,3,3)
histogram(lilLam3{end}); hold on;
histogram(lilF33{end})


magmaVals = colormap(magma);
figure(300+str2double(blurfilt(3)))
for i=1:length(lilJ)
   h{i} = histogram(lilJ{i}); hold on;
  h{i}.EdgeColor = 'none';
  h{i}.LineStyle = 'none';
  h{i}.FaceColor = magmaVals(round(i/length(lilJ)*64),:);
end


figure(400);
%Average J along the cross-section
for j=1:length(detF)    

    hold on;
    volsz = size(detF{j});
    Y=reshape(detF{j},volsz(1),volsz(2)*volsz(3))';
    X = X1;
    col = magmaVals(round(j/length(detF)*64),:);
%     iosr.statistics.functionalPlot(X,Y,'mainLineColor',magmaVals(round(j/length(detF)*64),:))
    iosr.statistics.functionalSpreadPlot(X,Y,'mainLineColor',col,'outerLineColor','none','spreadColor',col,'spreadBorderLineColor',col, ...
        'spreadAlpha',0.5,'outlierEdgeColor',col);
    hold on;
    
end

figure(500);
%Average Lambdas along the cross-section
for j=1:length(detF)    
    subplot(1,3,1);
    avgLam1_xA{j} = nansum(nansum(Lami{j}{1},3),2)./sum(sum(~isnan(Lami{j}{1}),3),2); hold on;
    plot((1:128)*0.43-64*.43,avgLam1_xA{j},'Color',magmaVals(round(j/length(detF)*64),:),'LineWidth',2);
    xlim([1, 128]*0.43-64*.43);
    
    subplot(1,3,2);
    avgLam2_xA{j} = nansum(nansum(Lami{j}{2},3),2)./sum(sum(~isnan(Lami{j}{2}),3),2); hold on;
    plot((1:128)*0.43-64*.43,avgLam2_xA{j},'Color',magmaVals(round(j/length(detF)*64),:),'LineWidth',2);
    xlim([1, 128]*0.43-64*.43);
    
    subplot(1,3,3);
    avgLam3_xA{j} = nansum(nansum(Lami{j}{3},3),2)./sum(sum(~isnan(Lami{j}{3}),3),2); hold on;
    plot((1:128)*0.43-64*.43,avgLam3_xA{j},'Color',magmaVals(round(j/length(detF)*64),:),'LineWidth',2);
    xlim([1, 128]*0.43-64*.43);
    
end

figure(501);
for j=1:length(detF)    
    subplot(1,3,1);
    hold on;
    volsz = size(Lami{j}{1});
    Y=reshape(Lami{j}{1},volsz(1),volsz(2)*volsz(3))';
    X = X1;
    col = magmaVals(round(j/length(detF)*64),:);
%     iosr.statistics.functionalPlot(X,Y,'mainLineColor',magmaVals(round(j/length(detF)*64),:))
    iosr.statistics.functionalSpreadPlot(X,Y,'mainLineColor',col,'outerLineColor','none','spreadColor',col,'spreadBorderLineColor',col, ...
        'spreadAlpha',0.5,'outlierEdgeColor',col);
    hold on;
    
    subplot(1,3,2);
    hold on;
    Y=reshape(Lami{j}{2},volsz(1),volsz(2)*volsz(3))';
    X = X1;
    col = magmaVals(round(j/length(detF)*64),:);
%     iosr.statistics.functionalPlot(X,Y,'mainLineColor',magmaVals(round(j/length(detF)*64),:))
    iosr.statistics.functionalSpreadPlot(X,Y,'mainLineColor',col,'outerLineColor','none','spreadColor',col,'spreadBorderLineColor',col, ...
        'spreadAlpha',0.5,'outlierEdgeColor',col);
    hold on;
    
    subplot(1,3,3);
    hold on;
    Y=reshape(Lami{j}{3},volsz(1),volsz(2)*volsz(3))';
    X = X1;
    col = magmaVals(round(j/length(detF)*64),:);
%     iosr.statistics.functionalPlot(X,Y,'mainLineColor',magmaVals(round(j/length(detF)*64),:))
    iosr.statistics.functionalSpreadPlot(X,Y,'mainLineColor',col,'outerLineColor','none','spreadColor',col,'spreadBorderLineColor',col, ...
        'spreadAlpha',0.5,'outlierEdgeColor',col);
    hold on;
    
    
%     subplot(1,3,2);
%     Y=reshape(Lami{j}{2},volsz(1),volsz(2)*volsz(3))';
%     X = X1;
%     iosr.statistics.functionalSpreadPlot(X,Y,'LineColor',magmaVals(round(j/length(detF)*64),:))
%     
%     subplot(1,3,3);
%     Y=reshape(Lami{j}{3},volsz(1),volsz(2)*volsz(3))';
%     X = X1;
%     iosr.statistics.functionalSpreadPlot(X,Y,'LineColor',magmaVals(round(j/length(detF)*64),:))    
end


figure(600);
for i=1:length(lilLam1)
    subplot(3,1,1);
    histogram(lilLam1{i}); hold on; 
    xlim([0.95 1.25]);
    subplot(3,1,2);
    histogram(lilLam2{i}); hold on; 
    xlim([0.90 1.05]);
    subplot(3,1,3);
    histogram(lilLam3{i}); hold on; 
    xlim([0.90 1.05]);
end

figure(601);
for i=1:length(lilLam1)
    subplot(3,1,1);
    histogram(lilF11{i}); hold on; 
    xlim([0.95 1.25]);
    subplot(3,1,2);
    histogram(lilF22{i}); hold on; 
    xlim([0.90 1.05]);
    subplot(3,1,3);
    histogram(lilF33{i}); hold on; 
    xlim([0.90 1.05]);
end


% figure(4)
% % Y=[rand(1000,1),gamrnd(1,2,1000,1),normrnd(10,2,1000,1),gamrnd(10,0.1,1000,1)];
% % violin(Y,'x',[-1 .7 3.4 8.8],'facecolor',[1 1 0;0 1 0;.3 .3 .3;0 0.3 0.1],'edgecolor','none',...
% % 'bw',0.3,'mc','k','medc','r-.')
% % axis([-2 10 -0.5 20])
% % ylabel('\Delta [yesno^{-2}]','FontSize',14)
% Y=[];
% for k=1:length(lilJ)
%     temp = lilJ{k}(:);
%     temp(isnan(temp))=[];
%     Y = cat(2,Y,temp);
% end
% X = thirdtrsig/max(thirdtrsig)*15;
% violin(Y,'x',X,'edgecolor','none','bw',0.3,'mc','k','medc','r-.')
