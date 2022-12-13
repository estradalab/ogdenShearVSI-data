function [model_3D] = createGeometry(l,w,h,line_res,hole,edge,graph,mesh_ref) % line_res >3

model = createpde(1);

funcs = createEdgeFuncs(edge,l,w);

temp = linspace(-l/2,l/2,line_res)';
edges{1} = temp; % Edge 1, x-coords
edges{2} = funcs{1}(temp); % Edge 1, y-coords

temp = -linspace((-w/2)+(w/(line_res-1)),(w/2)-(w/(line_res-1)),line_res-2)';
edges{2} = [edges{2};temp];
edges{1} = [edges{1};funcs{2}(temp)];

temp = -linspace(-l/2,l/2,line_res)';
edges{1} = [edges{1};temp];
edges{2} = [edges{2};funcs{3}(temp)];

temp = linspace((-w/2)+(w/(line_res-1)),(w/2)-(w/(line_res-1)),line_res-2)';
edges{2} = [edges{2};temp];
edges{1} = [edges{1};funcs{4}(temp)];

R1 = [2;4*(line_res-1);edges{1};edges{2}];
% R1 = [3;4;-l/2;l/2;l/2;-l/2;-w/2;-w/2;w/2;w/2];

switch hole.yesno
    case 'yes'
        switch hole.type
            case 'two_central_5mm'
                C1 = [1;-(l/2)+(l/3);0;2.5];
                C2 = [1;-(l/2)+(2*l/3);0;2.5];
                C1 = [C1;zeros(length(R1) - length(C1),1)];
                C2 = [C2;zeros(length(R1) - length(C2),1)];
                gm = [R1,C1,C2];
                sf = 'R1-C1-C2'; ns = char('R1','C1','C2'); ns = ns';
        end
    case 'no'
        gm = R1; sf = 'R1'; ns = char('R1'); ns = ns';
end

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