function uwPh = WeightedUnwrap3(Ph,w,mask)

d2Ph = secondOrderDiffs(Ph);
Rijk = voxelReliability(d2Ph,mask).*w;
Edgeijk = edgeReliability(Rijk);

uwPh = Unwrap3(Rijk,Edgeijk,Ph);

end

function d2Ph = secondOrderDiffs(Ph)

ijk = makeijk();
d2Ph = zeros([size(Ph), length(ijk)]);
PhPad = padarray(Ph,[1 1 1],nan);
for n = 1:length(ijk)
    d1L = wrapToPi(PhPad((2:(end-1))-ijk(n,1),(2:(end-1))-ijk(n,2),(2:(end-1))-ijk(n,3))-Ph);
    %d1L(d1L>pi) = d1L(d1L>pi)-2*pi; d1L(d1L<-pi) = d1L(d1L<-pi)+2*pi;
    d1R = wrapToPi(Ph-PhPad((2:(end-1))+ijk(n,1),(2:(end-1))+ijk(n,2),(2:(end-1))+ijk(n,3)));
    %d1R(d1R>pi) = d1R(d1R>pi)-2*pi; d1R(d1R<-pi) = d1R(d1R<-pi)+2*pi;
    d2Ph(:,:,:,n) = d1L-d1R;
    %     imagesc(d1L(:,:,60)-d1R(:,:,60)); axis image;
end

end

function Rijk = voxelReliability(d2Ph,mask)

edgeon = false;

Rijk = sum(d2Ph.^2,4).^(-1/2);
Rijk = Rijk./max(Rijk(:));

if edgeon
    Rijk(isnan(Rijk)&~isnan(mask))=0.00001;
    Rijk = Rijk.*sum(~isnan(d2Ph),4)/26;
else
    Rijk(isnan(Rijk)&~isnan(mask))=0;
end

% figure; imagesc(Rijk(:,:,60)); axis image; set(gca,'YDir','normal');
end

function Edgeijk = edgeReliability(Rijk)

edgeX = 2*movmean(Rijk,2,1);
edgeX = edgeX(2:end,:,:);

edgeY = 2*movmean(Rijk,2,2);
edgeY = edgeY(:,2:end,:);

edgeZ = 2*movmean(Rijk,2,3);
edgeZ = edgeZ(:,:,2:end);

% figure(2); subplot(1,3,1);
% imagesc(edgeX(:,:,60)); axis image; set(gca,'YDir','normal');
% subplot(1,3,2);
% imagesc(edgeY(:,:,60)); axis image; set(gca,'YDir','normal');
% subplot(1,3,3);
% imagesc(edgeZ(:,:,60)); axis image; set(gca,'YDir','normal');


Edgeijk{1} = edgeX;
Edgeijk{2} = edgeY;
Edgeijk{3} = edgeZ;

end

function uwPh = Unwrap3(Rijk,Edgeijk,Ph)
gifon = false;

%[Rsorted,sortIdx] = sort(Rijk(:));
AllEdges = [Edgeijk{1}(:); Edgeijk{2}(:); Edgeijk{3}(:)];
[EdgeSorted,sortIdx] = sort(AllEdges,'descend','MissingPlacement','last');
goodidx = sum(~isnan(EdgeSorted));

sizePh = size(Ph);
[mx,my,mz] = ndgrid(1:sizePh(1),1:sizePh(2),1:sizePh(3));
%sizeEdges = sizeU-[1 1 1];
exmx = movmean(mx,2,1); exmx = exmx(2:end,:,:); exmy = my(2:end,:,:); exmz = mz(2:end,:,:);
eymy = movmean(my,2,2); eymy = eymy(:,2:end,:); eymx = mx(:,2:end,:); eymz = mz(:,2:end,:);
ezmz = movmean(mz,2,3); ezmz = ezmz(:,:,2:end); ezmx = mx(:,:,2:end); ezmy = my(:,:,2:end);

%XYorZind = [ones(length(Edgeijk{1}(:)),1); 2*ones(length(Edgeijk{2}(:)),1); 3*ones(length(Edgeijk{3}(:)),1)];
ex = [exmx(:), exmy(:), exmz(:)];
ey = [eymx(:), eymy(:), eymz(:)];
ez = [ezmx(:), ezmy(:), ezmz(:)];
AllEijk = [ex; ey; ez];

clear AllEdges EdgeSorted mx my mz ex ey ez exmx exmy exmz eymx eymy eymz ezmx ezmy ezmz;

groupMx = zeros(sizePh);
uwPh = nan(sizePh);
for n = 1:goodidx
    %     if n==110651%158284
    %         pause;
    %     end
    
    edge_ijk = AllEijk(sortIdx(n),:);
    pt1 = floor(edge_ijk); pt2 = ceil(edge_ijk);
    
    
    groupMx1 = groupMx(pt1(1),pt1(2),pt1(3));
    groupMx2 = groupMx(pt2(1),pt2(2),pt2(3));
    
    Ph1 = Ph(pt1(1),pt1(2),pt1(3));
    Ph2 = Ph(pt2(1),pt2(2),pt2(3));
    
    % if sum((pt1==[193,52,70]&pt2==[193,52,71])| (pt2==[193,52,70]&pt1==[193,52,71]))==3
    %      [groupMx1 groupMx2]
    %      [Ph1 Ph2;]
    %      figure; imagesc(squeeze(uwPh(:,52,:)))
    % end
    
    %If neither voxel is assigned to a group
    if groupMx1+groupMx2 == 0
        %Add them both to a new group
        groupMx(pt1(1),pt1(2),pt1(3)) = max(groupMx(:))+1;
        groupMx(pt2(1),pt2(2),pt2(3)) = groupMx(pt1(1),pt1(2),pt1(3));
        uwPh(pt1(1),pt1(2),pt1(3)) = Ph1;
        
        dPh = Ph(pt2(1),pt2(2),pt2(3))-uwPh(pt1(1),pt1(2),pt1(3));
        %determine if any unwrap needs to be done (e.g. shift ~= 0)
        shift = (dPh-wrapToPi(dPh));
        %shift the phase of point 2 w.r.t. point 1
        uwPh(pt2(1),pt2(2),pt2(3)) = Ph(pt2(1),pt2(2),pt2(3))-shift;
        %combine the groups into point 1's group
        groupMx(pt2(1),pt2(2),pt2(3)) = groupMx1;
        
        %If one voxel is assigned and the other isn't
    elseif groupMx1==0 && groupMx2>0
        %Unwrap and add the unassigned voxel to the same group as the previously assigned one
        %calculate phase difference
        dPh = Ph(pt1(1),pt1(2),pt1(3))-uwPh(pt2(1),pt2(2),pt2(3));
        %determine if any unwrap needs to be done (e.g. shift ~= 0)
        shift = (dPh-wrapToPi(dPh));
        %shift the phase of point 1 w.r.t. point 2
        uwPh(pt1(1),pt1(2),pt1(3)) = Ph(pt1(1),pt1(2),pt1(3))-shift;
        %combine the groups into point 2's group
        groupMx(pt1(1),pt1(2),pt1(3)) = groupMx2;
    elseif groupMx2==0 && groupMx1>0
        %Unwrap and add the unassigned voxel to the same group as the previously assigned one
        %calculate phase difference
        dPh = Ph(pt2(1),pt2(2),pt2(3))-uwPh(pt1(1),pt1(2),pt1(3));
        %determine if any unwrap needs to be done (e.g. shift ~= 0)
        shift = (dPh-wrapToPi(dPh));
        %shift the phase of point 2 w.r.t. point 1
        uwPh(pt2(1),pt2(2),pt2(3)) = Ph(pt2(1),pt2(2),pt2(3))-shift;
        %combine the groups into point 1's group
        groupMx(pt2(1),pt2(2),pt2(3)) = groupMx1;
        
        %If both are assigned and are in different groups
    elseif (groupMx1>0 && groupMx2>0) && (groupMx1~=groupMx2)
        %Unwrap the smaller group with respect to the larger group
        %If a)group 1 is larger, or
        %   b)the groups are of equal size and group 1 has a lower index/is
        %   more reliable
        if (sum(groupMx(:)==groupMx1)>sum(groupMx(:)==groupMx2))...
                || (sum(groupMx(:)==groupMx1)==sum(groupMx(:)==groupMx2) && groupMx1<groupMx2)
            %calculate phase difference
            dPh = uwPh(pt2(1),pt2(2),pt2(3))-uwPh(pt1(1),pt1(2),pt1(3));
            %determine if any unwrap needs to be done (e.g. shift ~= 0)
            shift = (dPh-wrapToPi(dPh));
            %shift the phase of point 2's group w.r.t. point 1
            uwPh(groupMx==groupMx2) = uwPh(groupMx == groupMx2)-shift;
            %combine the groups into point 1's group
            groupMx(groupMx==groupMx2) = groupMx1;
        else %sum(groupMx==groupMx(pt1))<sum(groupMx==groupMx(pt2))
            %calculate phase difference
            dPh = uwPh(pt1(1),pt1(2),pt1(3))-uwPh(pt2(1),pt2(2),pt2(3));
            %determine if any unwrap needs to be done (e.g. shift ~= 0)
            shift = (dPh-wrapToPi(dPh));
            %shift the phase of point 1's group w.r.t. point 2
            uwPh(groupMx==groupMx1) = uwPh(groupMx == groupMx1)-shift;
            %combine the groups into point 2's group
            groupMx(groupMx==groupMx1) = groupMx2;
        end
    end
    
    
    % if groupMx1==3348 || groupMx2==3348
    %     [Ph1 Ph2; uwPh(pt1(1),pt1(2),pt1(3)) uwPh(pt2(1),pt2(2),pt2(3))]
    % end
    %
    
    if gifon
        z=round((96-63)/96*20);
        if n==1
            edge1 = AllEijk(sortIdx(1),:);
            pt1 = floor(edge1); pt2 = ceil(edge1);
            %figure(101); imagesc(uwPh(:,:,48)-uwPh(pt1(1),pt1(2),pt1(3))); truesize; set(gca,'ydir','normal'); colormap magma;
            figure(102); imagesc(uwPh(:,:,z)-nanmedian(uwPh(:))); truesize; set(gca,'ydir','normal'); colormap magma;
            caxis([-4*pi 4*pi]);
            fn = '1515mshift_wDESTEu.gif';
            pause(0.2);
            gif(fn,'LoopCount',Inf,'DelayTime',1/7);
            %imwrite(rgb2ind(frame2im(getframe),256),fn,'LoopCount',Inf,'DelayTime',1/10);
        end
        
        freq = 100; %was 10000
        if mod(n,freq)==0
            edge1 = AllEijk(sortIdx(1),:);
            pt1 = floor(edge1); pt2 = ceil(edge1);
            %figure(n/1000);
            %figure(101); imagesc(uwPh(:,:,48)-uwPh(pt1(1),pt1(2),pt1(3))); truesize; set(gca,'ydir','normal'); colormap magma;
            figure(102); imagesc(uwPh(:,:,z)-nanmedian(uwPh(:))); truesize; set(gca,'ydir','normal'); colormap magma;
            caxis([-4*pi 4*pi]);
            gif
            %pause(0.2);
            %imwrite(rgb2ind(frame2im(getframe),256),fn,'WriteMode','append','DelayTime',1/10);
        end
    
    edge1 = AllEijk(sortIdx(1),:);
    pt1 = floor(edge1); pt2 = ceil(edge1);
    %figure(101); imagesc(uwPh(:,:,48)-uwPh(pt1(1),pt1(2),pt1(3))); axis image; set(gca,'ydir','normal'); colormap magma; truesize
    figure(102); imagesc(uwPh(:,:,z)-nanmedian(uwPh(:))); truesize; set(gca,'ydir','normal'); colormap magma;
    caxis([-4*pi 4*pi]);
    gif
    %pause(0.2);
    %imwrite(rgb2ind(frame2im(getframe),256),fn,'WriteMode','append','DelayTime',1/10);

    end
end
end

function ijk = makeijk()
ijk = [1 0 0; 0 1 0; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1; ...
    1 1 0; 1 0 1; 0 1 1; 1 -1 0; 1 0 -1; 0 1 -1;...
    -1 -1 0; -1 0 -1; 0 -1 -1; -1 1 0; -1 0 1; 0 -1 1;...
    1 1 1; -1 1 1; 1 -1 1; 1 1 -1;...
    -1 -1 -1; 1 -1 -1; -1 1 -1; -1 -1 1;];
end