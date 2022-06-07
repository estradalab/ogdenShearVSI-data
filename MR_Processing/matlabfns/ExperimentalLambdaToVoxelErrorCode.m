


for m = 1:size(Eij{3,3},1)
    eps11plane = Eij{3,3}(m,:,:);
    std_eps11(m) = nanstd(eps11plane(:));
    med_eps11(m) = nanmedian(eps11plane(:));
end
%figure; plot(med_eps11,'-o','LineWidth',2);
%hold on; errorbar(smooth(med_eps11,0.1,'sgolay'),std_eps11);

%figure; errorbar((med_eps11-smooth(med_eps11,0.1,'sgolay')')./(smooth(med_eps11,0.1,'sgolay')'),std_eps11./smooth(med_eps11,0.1,'sgolay')');
%set(gca,'ylim',[-1 1]);

useprev = false;
if ~useprev
    smoothedcurve = smooth(med_eps11,0.15,'sgolay')';
end


lims = round([50 90]/128*size(data,1)); %[50 90] for 0906, [60 150] for 0831 
errorVals = (med_eps11-smoothedcurve);%./smoothedcurve;
stdevVals = std(med_eps11(lims(1):lims(2))-smoothedcurve(lims(1):lims(2)))%std_eps11./smoothedcurve;
rmsVals = rms(med_eps11(lims(1):lims(2))-smoothedcurve(lims(1):lims(2)))
ratio = 2*lambda(1)/abs(axis1(2)-axis1(1));%*30/(30+wiggle(3)) %should be lambda_e from the MR/voxel sz
medError = median(abs(errorVals(lims(1):lims(2))));
medStdev = median(stdevVals)
medRms = median(rmsVals)

if isempty(wiggle)
    wiggle = [0 0 0];
end

magmavals = magma(64);

 figure(10); 
h = shadedErrorBar(1:length(med_eps11),med_eps11,std_eps11,'lineProps',{'-','Color',magmavals(round(63*wiggle(3)/7+1),:)}); hold on;
% plot(med_eps11,'-o','LineWidth',2,'lineProps',{'-','markerfacecolor',magmavals(round(63*wiggle(3)/7+1))});
 %hold on; errorbar(smoothedcurve,std_eps11); %hold off;
 
 xlim([40 220]); ylim([-0.1 0.3]);%ylim([-0.02 0.02]);
%figure; errorbar((med_eps11-smoothedcurve)./smoothedcurve,std_eps11./smoothedcurve);
%set(gca,'ylim',[-1 1]);

%%



figure(1);
hold on;

if size(data,1) == 128
    str = 'd';
elseif size(data,1) == 192
    str = '+';
elseif size(data,1) == 96
    str = 'o';
elseif size(data,1) == 256
    str = 'v';
else
    str = '^';
end

try
plot(ratio,medRms,str,'Color',magmavals(round(63*wiggle(3)/7+1),:),'LineWidth',2);
catch
    plot(ratio,medRms,str,'Color',magmavals(1,:),'LineWidth',2);
end
xlabel('Position [px]');
ylabel('Displacement gradient, u_{1,1} [ ]');
title('Non-nicked sample, Hamming filter off and 5-point optimal tap CD');
xlim([10 240]);

figure(2);
hold on;
str = 'd';
plot(wiggle(3),medRms,str,'Color',magmavals(round(63*wiggle(3)/7+1),:),'LineWidth',2);
xlabel('End-to-end displacement [mm]');
ylabel('Displacement gradient error magnitude, |\Delta u_{1,1}| [ ]');



clear all;