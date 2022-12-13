% v2 - fully parameterizable geometry
%      - implement ellipitical  notches (both symmetric pairs and single)

function [model_3D] = createGeometry_v2(l,w,h,hole,graph,mesh_ref,theta,db_lim,b) % line_res >3

model = createpde(1);

R1 = [3;4;-l/2;l/2;l/2;-l/2;-w/2;-w/2;w/2;w/2];

if w/2 - b*tan(theta)/2 < db_lim
    % Dogbone
    x = abs(w/2-db_lim)/tan(theta);
    Crack1 = [2;4;-b;-b+x;b-x;b;w/2;db_lim;db_lim;w/2];
    Crack2 = [2;4;-b;-b+x;b-x;b;-w/2;-db_lim;-db_lim;-w/2];
    gm = [R1,Crack1,Crack2];
    sf = 'R1-Crack1-Crack2'; ns = char('R1','Crack1','Crack2'); ns = ns';
elseif b*tan(theta)/2 > 0
    % Triangle cracks
    Crack1 = [2;3;-b;0;b;w/2;w/2-b*tan(theta)/2;w/2];
    Crack1 = [Crack1;zeros(length(R1) - length(Crack1),1)];
    Crack2 = [2;3;-b;0;b;-w/2;-w/2+b*tan(theta)/2;-w/2];
    Crack2 = [Crack2;zeros(length(R1) - length(Crack2),1)];
    gm = [R1,Crack1,Crack2];
    sf = 'R1-Crack1-Crack2'; ns = char('R1','Crack1','Crack2'); ns = ns';
else
    % Rectangular
    gm = R1; sf = 'R1'; ns = char('R1'); ns = ns';
end

% switch hole.yesno
%     case 'yes'
%         switch hole.type
%             case 'two_central_5mm'
%                 C1 = [1;-(l/2)+(l/3);0;2.5];
%                 C2 = [1;-(l/2)+(2*l/3);0;2.5];
%                 C1 = [C1;zeros(length(R1) - length(C1),1)];
%                 C2 = [C2;zeros(length(R1) - length(C2),1)];
%                 gm = [R1,C1,C2];
%                 sf = 'R1-C1-C2'; ns = char('R1','C1','C2'); ns = ns';
%         end
%     case 'no'
%         % gm = R1; sf = 'R1'; ns = char('R1'); ns = ns';
% 
%         % Notch test:
%         E1 = [4;0;-w/2;w/40;w/4;0];
%         E1 = [E1;zeros(length(R1) - length(E1),1)];
%         gm = [R1,E1];
%         sf = 'R1-E1'; ns = char('R1','E1'); ns = ns';
% 
% end

g = decsg(gm,sf,ns);
geometryFromEdges(model,g);
if graph
subplot(1,3,1)
pdegplot(model,"EdgeLabels","off")
axis([-l/1.3 l/1.3 -w/1.3 w/1.3])
axis off
end

model_3D = model;
model_3D.Geometry = extrude(model_3D.Geometry,h);
if graph
subplot(1,3,2)
pdegplot(model_3D,"FaceLabels","off","FaceAlpha",0.5)
axis off
end


generateMesh(model_3D);
hmax = model_3D.Mesh.MaxElementSize;
generateMesh(model_3D,'Hmax',hmax*mesh_ref);
if graph
subplot(1,3,3)
pdeplot3D(model_3D)
axis off
end


%% MATLAB In-line Commands
 % Dogbone
% l = 40; w = 8; h = 8.5; hole = hole; graph = true; mesh_ref = 0.7; theta = pi/4; db_lim = w/3; b = l/10;
% [model_3D] = createGeometry_v2(l,w,h,hole,graph,mesh_ref,theta,db_lim,b);

% Rectangular
% l = 40; w = 8; h = 8.5; hole = hole; graph = true; mesh_ref = 0.7; theta = pi/4; db_lim = w/3; b = 0;
% [model_3D] = createGeometry_v2(l,w,h,hole,graph,mesh_ref,theta,db_lim,b);

% Crack
% l = 40; w = 8; h = 8.5; hole = hole; graph = true; mesh_ref = 0.7; theta = pi/4; db_lim = w/6; b = l/10;
% [model_3D] = createGeometry_v2(l,w,h,hole,graph,mesh_ref,theta,db_lim,b);