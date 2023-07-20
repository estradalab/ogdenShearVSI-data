function [model_3D] = createGeometry(l,w,h,line_res,edge,mesh_ref,elementType) % line_res >3

maxelsize = mesh_ref.maxelsize;

model = createpde(1);

funcs = createEdgeFuncs(edge,l,w);

temp = linspace(-l/2,l/2,line_res(1))';
edges{1} = temp; % Edge 1, x-coords
edges{2} = funcs{1}(temp); % Edge 1, y-coords

temp = -linspace((-w/2)+(w/(line_res(2)-1)),(w/2)-(w/(line_res(2)-1)),line_res(2)-2)';
edges{2} = [edges{2};temp];
edges{1} = [edges{1};funcs{2}(temp)];

temp = -linspace(-l/2,l/2,line_res(3))';
edges{1} = [edges{1};temp];
edges{2} = [edges{2};funcs{3}(temp)];

temp = linspace((-w/2)+(w/(line_res(4)-1)),(w/2)-(w/(line_res(4)-1)),line_res(4)-2)';
edges{2} = [edges{2};temp];
edges{1} = [edges{1};funcs{4}(temp)];

R1 = [2;4*(mean(line_res)-1);edges{1};edges{2}];

gm = R1; sf = 'R1'; ns = char('R1'); ns = ns';

g = decsg(gm,sf,ns);
geometryFromEdges(model,g);

model_3D = model;
model_3D.Geometry = extrude(model_3D.Geometry,h);

switch elementType
    case 'C3D4H'
        generateMesh(model_3D,'GeometricOrder','linear');
    case 'C3D10H'
        generateMesh(model_3D);
end
hmax = model_3D.Mesh.MaxElementSize;

switch elementType
    case 'C3D4H'
        if exist('maxelsize','var') == 1
            generateMesh(model_3D,'Hmax',maxelsize,'Hmin',maxelsize,'GeometricOrder','linear');
        else
            generateMesh(model_3D,'Hmax',hmax*mesh_ref.defsize,'Hmin',maxelsize,'GeometricOrder','linear');
        end
    case 'C3D10H'
        if exist('maxelsize','var') == 1
            generateMesh(model_3D,'Hmax',maxelsize,'Hmin',maxelsize);
        else
            generateMesh(model_3D,'Hmax',hmax*mesh_ref.defsize,'Hmin',maxelsize);
        end
end