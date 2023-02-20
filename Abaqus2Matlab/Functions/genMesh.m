function  [x,y,z,TRI] = genMesh(curdir,fileName,mesh,params,write,pres_disp,mesh_ref)

fN = erase(fileName,'_tet');
graph = false;

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
switch true
    case (startsWith(fN,'Uniaxial')==1)
        l = 40; w = 8; h = 8.5;
        hole.yesno = 'no';
        edge.shape{1} = 'line'; edge.shape{2} = 'line'; edge.shape{3} = 'line'; edge.shape{4} = 'line';
        line_res = [4 4 4 4]; 
    case (startsWith(fN,'5MMHoles')==1)
        l = 25.4; w = 19.1; h = 6.4;
        hole.yesno = 'yes'; hole.type = 'two_central_5mm';
        edge.shape{1} = 'line'; edge.shape{2} = 'line'; edge.shape{3} = 'line'; edge.shape{4} = 'line';
        line_res = [4 4 4 4]; 
    case (startsWith(fN,'NoHoles')==1)
        l = 25.4; w = 19.1; h = 6.4;
        hole.yesno = 'no';
        edge.shape{1} = 'line'; edge.shape{2} = 'line'; edge.shape{3} = 'line'; edge.shape{4} = 'line';
        line_res = [4 4 4 4]; 
    case (startsWith(fN,'ShearWavy')==1)
        if startsWith(fN,'ShearWavyAmp2') || startsWith(fN,'ShearWavyPrd2')
            l = 40; w = 8; h = 30;
        elseif startsWith(fN,'ShearWavyWidth1')
            l = 40; w = 8;
            h = str2double(erase(fN,'ShearWavyWidth1_6.25MMDisp_Amp_0.6MM_Prd_4_Width_'));
        else
            l = 40; w = 8; h = 7.5;
        end
        hole.yesno = 'no';
        if startsWith(fN,'ShearWavyPrd1') || startsWith(fN,'ShearWavyPrd2')
            edge.coef{1} = -0.6;
            edge.coef{3} = 0.6;
            edge.shape{1} = 'sin'; edge.shape{2} = 'line'; edge.shape{3} = 'sin'; edge.shape{4} = 'line';
            line_res = [101 4 101 4]; 
            edge.period = str2double(erase(fN,'ShearWavyPrd1_6.25MMDisp_Amp_0.6MM_Prd_'));
            if isnan(edge.period) == 1
                edge.period = str2double(erase(fN,'ShearWavyPrd2_6.25MMDisp_Amp_0.6MM_Prd_'));
            end
        elseif startsWith(fN,'ShearWavyAmp2')
            edge.coef{1} = -str2double(erase(fN,'ShearWavyAmp2_6.25MMDisp_Amp_'));
            edge.coef{3} = str2double(erase(fN,'ShearWavyAmp2_6.25MMDisp_Amp_'));
            if str2double(erase(fN,'ShearWavyAmp2_6.25MMDisp_Amp_')) == 0
                edge.shape{1} = 'line'; edge.shape{2} = 'line'; edge.shape{3} = 'line'; edge.shape{4} = 'line';
                line_res = [4 4 4 4];
            else
                edge.shape{1} = 'sin'; edge.shape{2} = 'line'; edge.shape{3} = 'sin'; edge.shape{4} = 'line';
                line_res = [101 4 101 4]; 
            end
            edge.period = 5; % edge.amptol = 10^-3;
        elseif startsWith(fN,'ShearWavyWidth1')
            edge.coef{1} = -0.6;
            edge.coef{3} = 0.6;
            edge.period = 4;
            edge.shape{1} = 'sin'; edge.shape{2} = 'line'; edge.shape{3} = 'sin'; edge.shape{4} = 'line';
            line_res = [101 4 101 4]; 
        else
            edge.coef{1} = -str2double(erase(fN,'ShearWavy_6.25MMDisp_Amp_'));
            edge.coef{3} = str2double(erase(fN,'ShearWavy_6.25MMDisp_Amp_'));
            if str2double(erase(fN,'ShearWavy_6.25MMDisp_Amp_')) == 0
                edge.shape{1} = 'line'; edge.shape{2} = 'line'; edge.shape{3} = 'line'; edge.shape{4} = 'line';
                line_res = [4 4 4 4];
            else
                edge.shape{1} = 'sin'; edge.shape{2} = 'line'; edge.shape{3} = 'sin'; edge.shape{4} = 'line';
                line_res = [101 4 101 4]; 
            end
            edge.period = 5; % edge.amptol = 10^-3;
        end
        edge.func{1} = @(x) edge.coef{1}*sin(2*pi*edge.period*x/abs(l)) + w/2;
        edge.func{3} = @(x) edge.coef{3}*sin(2*pi*edge.period*x/abs(l)) - w/2;
end

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
switch true
    case {(startsWith(fN,'Uniaxial')==1),(startsWith(fN,'5MMHoles')==1),(startsWith(fN,'NoHoles')==1)}
        Surf.x1 = find(and(x<=max(x)+10^-9,x>=max(x)-10^-9));
        Surf.x2 = find(and(x<=min(x)+10^-9,x>=min(x)-10^-9));
        Surf.y1 = find(and(y<=max(y)+10^-9,y>=max(y)-10^-9));
        Surf.y2 = find(and(y<=min(y)+10^-9,y>=min(y)-10^-9));
        Surf.z1 = find(and(z<=max(z)+10^-9,z>=max(z)-10^-9));
        Surf.z2 = find(and(z<=min(z)+10^-9,z>=min(z)-10^-9));
    case (startsWith(fN,'ShearWavy')==1)
        Surf.x1 = find(and(x<=max(x)+10^-9,x>=max(x)-10^-9));
        Surf.x2 = find(and(x<=min(x)+10^-9,x>=min(x)-10^-9));
        Surf.z1 = find(and(z<=max(z)+10^-9,z>=max(z)-10^-9));
        Surf.z2 = find(and(z<=min(z)+10^-9,z>=min(z)-10^-9));
        if str2double(erase(fN,'ShearWavy_6.25MMDisp_Amp_')) == 0 || str2double(erase(fN,'ShearWavyAmp2_6.25MMDisp_Amp_')) == 0'
            Surf.y1 = find(and(y<=max(y)+10^-9,y>=max(y)-10^-9));
            Surf.y2 = find(and(y<=min(y)+10^-9,y>=min(y)-10^-9));
        else
            Surf.y1 = []; Surf.y2 = [];
            for j = 1:length(Nodes.gen)
                if ismembertol(y(j),edge.func{1}(x(j)))
                    Surf.y1 = [Surf.y1;j];
                elseif ismembertol(y(j),edge.func{3}(x(j)))
                    Surf.y2 = [Surf.y2;j];
                end
            end
        end
end

% Boundary conditions and surfaces
switch true
    case (startsWith(fN,'Uniaxial')==1)
        Nodes.bc1 = Surf.x1;
        Nodes.bc2 = Surf.x2;
        Nodes.presDisp.dir = 'x';
    case {(startsWith(fN,'5MMHoles')==1),(startsWith(fN,'NoHoles')==1)}
        Nodes.bc1 = Surf.z1;
        Nodes.bc2 = Surf.z2;
        Nodes.presDisp.dir = 'x';
    case (startsWith(fN,'ShearWavy')==1)
        Nodes.bc1 = Surf.y1;
        Nodes.bc2 = Surf.y2;
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

% Organizes optimization 'sweeps' into subfolders
jN = fN;
if contains(curdir,'sweep','IgnoreCase',true)
    fN = [curdir '\' fN];
    if ~exist(['InputFiles\' curdir], 'dir')
        mkdir(['InputFiles\' curdir])
        mkdir(['Data\' curdir])
    end
end

switch write
    case 'Write'
        writeInp(Nodes,Elements,Elements_Sets,fN,coef,jN)
    otherwise
end