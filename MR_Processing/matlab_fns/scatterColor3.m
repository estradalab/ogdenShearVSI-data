function f = scatterColor3(varargin)
%ScatterColor3 - for a gridded 3D scatterplot binned into colors, with
%transparency highlights for high values

%Syntax: f = scatterColor3(XYZC, numbins, colormap, gamma, maxalpha, markerSz, camerapos);
% XYZC = three position coordinates and data.
% size(X) = size(Y) = size(Z) = size(C)
% numbins = number of total color bins for the plot
% colormap = the chosen color scheme for the plot e.g. parula
% gamma = gamma value for shifting transparency values, with 0 being flat,
%   1 being linear increase, -1 being linear decrease
% maxalpha = (0,1] transparency maximum for the extreme data value
% markerSz = size of the point markers for the scatterplot
% camerapos = position of the camera
% climits = colorbar limits if user wants them defined, (min/max default)
% 

% Copyright 2018 Jonathan B. Estrada, Ph.D.
% Research Fellow, Mechanical Engineering, University of Michigan


[XYZC, numbins, colormap, gamma, maxalpha, markerSz, camerapos, climits] = parseInputs(varargin{:});

X = XYZC{1}; Y = XYZC{2}; Z = XYZC{3}; C = XYZC{4};

%If there are nans, take them out before processing to save time
C_ = ~isnan(C);

if isempty(climits)
    C_sc = ceil((C-min(C(:)))/(max(C(:))-min(C(:)))*(size(colormap,1)-1)+1);
else
    C_sc = ceil((C-min(climits(:)))/(max(climits(:))-min(climits(:)))*(size(colormap,1)-1)+1);
end

C1 = nan(size(C_sc)); C2 = C1; C3 = C1;

figure; hold on;
for nt = 1:numbins
    Cbin =  C_sc>((nt-1)/numbins*size(colormap,1)) & C_sc<=(nt/numbins*size(colormap,1));
    idxs = find(C_ & Cbin);
    
    for cc = 1:sum(Cbin(:))
        C1(idxs(cc)) = colormap(C_sc(idxs(cc)),1);
        C2(idxs(cc)) = colormap(C_sc(idxs(cc)),2);
        C3(idxs(cc)) = colormap(C_sc(idxs(cc)),3);
    end
    
    f{nt} = scatter3(X(idxs),Y(idxs),Z(idxs),markerSz,cat(2,C1(idxs),C2(idxs),C3(idxs)),'filled');
    if gamma>=0
        f{nt}.MarkerEdgeAlpha = (nt/numbins)^gamma*maxalpha;
        f{nt}.MarkerFaceAlpha = (nt/numbins)^gamma*maxalpha;
    else
        f{nt}.MarkerEdgeAlpha = (1-(nt/numbins))^abs(gamma)*maxalpha;
        f{nt}.MarkerFaceAlpha = (1-(nt/numbins))^abs(gamma)*maxalpha;
    end
end
daspect([1 1 1]);
campos(camerapos);
hold off;

end

function [data, numbins, colormap, gamma, maxalpha, markerSz, camerapos, climits] = parseInputs(varargin)

data = varargin{1};

if length(varargin) < 2, numbins = 11; else, numbins = varargin{2}; end

if length(varargin) < 3, colormap = parula; else, colormap = varargin{3}; end

if length(varargin) < 4, gamma = 1; else, gamma = varargin{4}; end

if length(varargin) < 5, maxalpha = 0.5; else, maxalpha = varargin{5}; end

if length(varargin) < 6, markerSz = 2; else, markerSz = varargin{6}; end

if length(varargin) < 7, camerapos =  [90 180 105]; else, camerapos = varargin{7}; end

if length(varargin) < 8, climits =  []; else, climits = varargin{8}; end

end