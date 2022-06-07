%clear all;

rsize = 1000;
p0 = rand(rsize,1)*2*pi;

%make 16 wraps of 2pi, starting at a random phase
x = 0:(2*pi/7200):32*pi;
dx = x(2)-x(1);

lambda = zeros(length(x),length(p0));
Lambda = zeros(length(x),length(p0));
for i=1:length(p0)
    lambda(:,i) = abs(sin(x-p0(i)))+abs(cos(x-p0(i)));
%    Lambda(:,i) = sin(x-p0(i)).*sign(cos(x-p0(i)))-cos(x-p0(i)).*sign(sin(x-p0(i)));
end

%always want to take 16 voxels, with a maximum ratio of vx/lambda = 1

Ints = zeros(16,rsize);

ratios = [1:0.01:1.33 1.35:0.05:3, 3.1:0.1:10 10.5:0.5:20;]
for rIdx = 1:length(ratios)
vx2lam = 1/ratios(rIdx);
    
for d = 1:16
    AngLims = [(d-1)*2*pi d*2*pi]*vx2lam;
    inIdx = x>AngLims(1) & x<=AngLims(2);
    xlimits = [find(x>AngLims(1),1), find(x<=AngLims(2),1,'last')];
    Ints(d,:) = sum(lambda(xlimits(1):xlimits(2),:))*dx;
    
end
rangeInts(rIdx) = (max(Ints(:))-min(Ints(:)));
stdInts(rIdx) = std(Ints(:));
rmsInts(rIdx) = rms(Ints(:)-mean(Ints(:)));
normstdInts(rIdx) = std(Ints(:))/median(Ints(:));
relI(rIdx) = rangeInts(rIdx)/median(Ints(:));

end

plot(ratios,relI,'-o','LineWidth',2); hold on;
plot(ratios,normstdInts,'-o','LineWidth',2,'Color','blue');
ylabel('St.Dev., MR Signal Acquisition');
xlabel('Ratio of \lambda_{enc}/voxel length [ ]');
%legend('Max range','StDev of MR Acq');

