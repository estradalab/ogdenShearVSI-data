function  [x,y,z,TRI] = genMesh(fileName,mesh,params,write,pres_disp,mesh_ref)

fN = erase(fileName,'_tet');

% Using meshgrid
% switch fN
%     case 'Uniaxial_7MMDisp'
%         X = -20:40/93:20; Y = -4:0.5:4; Z = 0:0.5:8.5; % Uniaxial 40X8X8.5 mm^3
% end
% 
% switch params
%     case 'ogden-treloar'
%         coef.model = 'Og_3';
%         coef.val = [0.4017 1.3 0.00295 5 0.00981 -2 0 0 0];
%                       % [mu_1 a_1 mu_2 a_2 mu_3 a_3 k_1 k_2 k_3]
% end
% 
% [xx, yy,zz]=meshgrid(X,Y,Z);
% x=reshape(xx,[length(X)*length(Y)*length(Z) 1]);
% y=reshape(yy,[length(X)*length(Y)*length(Z) 1]);
% z=reshape(zz,[length(X)*length(Y)*length(Z) 1]);
% TRI = delaunay(x,y,z);

% Using createGeometry
switch fN
    case 'Uniaxial_7MMDisp'
        l = 40; w = 8; h = 8.5;
        hole.yesno = 'no';
    case {'5MMHoles_2_5MMDisp','5MMHoles_5MMDisp'}
        l = 25.4; w = 19.1; h = 6.4;
        hole.yesno = 'yes'; hole.type = 'two_central_5mm';
    case {'NoHoles_2_5MMDisp','NoHoles_5MMDisp','NoHoles_7MMDisp'}
        l = 25.4; w = 19.1; h = 6.4;
        hole.yesno = 'no';
end

edge.shape{1} = 'line'; edge.shape{2} = 'line'; edge.shape{3} = 'line'; edge.shape{4} = 'line';
line_res = 4; graph = false;

model_3D = createGeometry(l,w,h,line_res,hole,edge,graph,mesh_ref);
x = model_3D.Mesh.Nodes(1,:)'; y = model_3D.Mesh.Nodes(2,:)'; 
z = model_3D.Mesh.Nodes(3,:)';
TRI = model_3D.Mesh.Elements';

switch params
    case 'ogden-treloar'
        coef.model = 'Og_3';
        coef.val = [0.4017 1.3 0.00295 5 0.00981 -2 0 0 0];
                      % [mu_1 a_1 mu_2 a_2 mu_3 a_3 k_1 k_2 k_3]
end

Nodes.gen=[x y z];                                   %(N*2) Nodes Coordinates 
% Different surfaces for boundary conditions
Surf.x1 = find(and(x<=max(x)+10^-9,x>=max(x)-10^-9));
Surf.x2 = find(and(x<=min(x)+10^-9,x>=min(x)-10^-9));
Surf.y1 = find(and(y<=max(y)+10^-9,y>=max(y)-10^-9));
Surf.y2 = find(and(y<=min(y)+10^-9,y>=min(y)-10^-9));
Surf.z1 = find(and(z<=max(z)+10^-9,z>=max(z)-10^-9));
Surf.z2 = find(and(z<=min(z)+10^-9,z>=min(z)-10^-9));

% Boundary conditions and surfaces
switch fN
    case 'Uniaxial_7MMDisp'
        Nodes.bc1 = Surf.x1;
        Nodes.bc2 = Surf.x2;
        Nodes.presDisp.dir = 'x';
    case {'5MMHoles_2_5MMDisp','5MMHoles_5MMDisp','NoHoles_2_5MMDisp','NoHoles_5MMDisp','NoHoles_7MMDisp'}
        Nodes.bc1 = Surf.z1;
        Nodes.bc2 = Surf.z2;
        Nodes.presDisp.dir = 'x';
end
Nodes.presDisp.mag = pres_disp;


for i=1:1:size(TRI,1) 
    Elements{i}=TRI(i,:);                          %Nodes indices vector for each Elements{i}              
end

Elements_Sets{1}.Name='Set-1';                  %Element set name
switch mesh
    case 'hex'
        Elements_Sets{1}.Elements_Type='C3D8RH';
    case 'tet'
        Elements_Sets{1}.Elements_Type='C3D4H';
end
Elements_Sets{1}.Elements=1:size(TRI,1);               %Elements indices vectors in the element set

switch write
    case 'Write'
        writeInp(Nodes,Elements,Elements_Sets,fN,coef)
    otherwise
end