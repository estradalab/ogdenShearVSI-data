%prpocess DESTE wkspcs
clear
close all

%lambda over rox =2, variable stretch ampounts
timestampvec={'1134','1142','1150','1158','1206','1214','1222','1230','1237'};

%vary RO #, lambda over rox =2, stretch is 4 and 0, 4 acquisitions each
%timestampvec={'1448','1456','1504','1512','1520','1526','1533','1539'};

for jj=1:numel(timestampvec);
    vs(jj)=load(['DESTE_1d_' timestampvec{jj} '.mat']);  %vs = var stretch
    display([num2str(jj) ': ' vs(jj).comment]);
end

%%
close all;

sliceind=[9:14];
numsl=numel(sliceind);

% in the sample
rowvec=30:100;
colvec=10:24;

%plot it all;
for jj=1:numel(vs);
    
    figure('Position',[100 100 1800 800],'name',vs(jj).comment);
    
    for kk=1:numsl;
        
        sli=sliceind(kk);
        
        subplot(2,numsl+1,kk);
        title(['slice ' num2str(sli)]);
        imagesc(angle(vs(jj).origdata(:,:,sli)));
        
        subplot(2,numsl+1,kk+numsl+1);
        title(['slice ' num2str(sli)]);
        imagesc(log(abs(vs(jj).origdata(:,:,sli))),[5 10]);
    end
end


%% Figure to explain the method, using dataset indexed by ds, and slice sl
close all;
for ds=1:numel(vs);
    sl=10;
    si=size(vs(ds).origdata(:,:,sl));
    
    rowvec=41:90;
    colvec=7:28;
    
    %initialize the rotated data
    outmatrix=zeros(si);
    rotmatrix=zeros(si);
    strerr=zeros(si);
    
    test=vs(ds).origdata(:,:,sl);                   %the full slice
    
    %fit phase of the slice;
    
    dummy=test(rowvec,colvec);       %the slice in the center, rest is nan's
    [outmatrix(rowvec,colvec),rotmatrix(rowvec,colvec),aux]=gradshim(dummy);      %shim on the center
    
    %estimate the error
    [REFX,REFY]=gradient(aux.bowl);                 %aux
    [FX,FY]=gradient(angle(outmatrix(rowvec,colvec)));
    [dummyFX,dummyFY]=gradient(angle(dummy));
    strainerror=(FY+eps)./(REFY+eps);
    
    
    figure('position', [ 100 100 1200 800],'DefaultRectangleEdgeColor',[1 1 1],...
        'DefaultRectangleLinewidth',2,'Name', vs(ds).comment);
    
    subplot(3,4,1);
    imagesc(log(abs(test)));
    hold on;
    rectangle('position', [colvec(1),rowvec(1), colvec(end)-colvec(1), rowvec(end)-rowvec(1)]);
    title('log(magnitude)');
    
    subplot(3,4,2);
    imagesc(angle(test));
    hold on;
    h=rectangle('position', [colvec(1),rowvec(1), colvec(end)-colvec(1), rowvec(end)-rowvec(1)]);
    title('phase');
    
    subplot(3,4,3);
    imagesc(colvec,rowvec,diff(angle(dummy)),[-pi pi]);
    title('grad_y(phase data), \pm 180 deg');
    hold on;
    
    subplot(3,4,4);
    imagesc(colvec,rowvec,diff(angle(dummy),1,2),[-pi pi]/18)
    title('grad_x(phase data) \pm 10 deg');
    hold on;
    
    
    
    subplot(3,4,5);
    imagesc(colvec,rowvec,angle(outmatrix(rowvec,colvec)),[-pi pi]/4);
    title('phase of residual, \pm 45 deg');
    %set(gca,'Position', [ 0.4   0.1100    0.1566    0.3412]);
    %hold on; h=rectangle('position', [colvec(1),rowvec(1), colvec(end)-colvec(1), rowvec(end)-rowvec(1)]);
    
    subplot(3,4,10);
    imagesc(colvec,rowvec,FY,[-pi pi]/18);
    title('grad_y(phase residual), scale \pm 10 deg');
    %hold on; h=rectangle('position', [colvec(1),rowvec(1), colvec(end)-colvec(1), rowvec(end)-rowvec(1)]);
    
    
    subplot(3,4,6);
    dummy=rotmatrix;
    dummy(isnan(outmatrix))=NaN;
    imagesc(angle(dummy));
    title('fit to exp phase, 3rd order taylor.');
    %set(gca,'Position', [ 0.45    0.4096    0.1566    0.2157]);
    %hold on; h=rectangle('position', [colvec(1),rowvec(1), colvec(end)-colvec(1), rowvec(end)-rowvec(1)]);
    
    subplot(3,4,8);
    dummy=REFY;
    dummy(isnan(outmatrix))=NaN;
    imagesc(dummy,[-pi pi]);
    title('grad_y(fitted phase), \pm 180 deg');
    %set(gca,'Position', [ 0.45    0.4096    0.1566    0.2157]);
    %hold on; h=rectangle('position', [colvec(1),rowvec(1), colvec(end)-colvec(1), rowvec(end)-rowvec(1)]);
    
    subplot(3,4,12);
    imagesc(colvec,rowvec,strainerror,[-0.1 0.1]);
    hold on;
    title('grad_y(residual)/grad_y(fit)');
    %set(gca,'Position', [0.7    0.1100    0.1491    0.2157])
    
    
    %text output
    subplot(3,4,11);
    set(gca,'Visible','off');
    
    relerror=std(strainerror(~isnan(strainerror(:))));
    h=text(0.1,1, ['rel. strain error =' num2str(round(relerror*100)/100)],'color',[0 0 0]);
    
    abserror=std(FY(:))/(2*pi) *vs(ds).lambda /(vs(ds).axis1(1)-vs(ds).axis1(2));
    h=text(0.1,0.9, ['abs. strain error =' num2str(round(abserror*1000)/1000)],'color',[0 0 0]);
    
    %meshplot of the unwrapped phase and the fit
    subplot(3,4,9);
    mesh(colvec,rowvec,angle(outmatrix(rowvec,colvec)));
    view([-70 50]);
    title('phase residual');
    set(gca,'zlim',[-0.5 0.5]);
    
    subplot(3,4,7)
    mesh(colvec,rowvec,aux.bowl);
    view([-70 50]);
    title('fit to exp. phase');
    set(gca,'zlim',[-50 50]);
end


