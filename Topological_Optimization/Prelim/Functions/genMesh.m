function  [x,y,z,TRI] = genMesh(fileName,mesh,params,pres_disp,...
    mesh_ref,abaqus_ver,elementType,d,N,A,l)

fN = erase(fileName,['_' mesh]);

% Length parameters
w = d; h = d;

% Determine maximum element size
mesh_ref.maxelsize = nthroot(l*w*h*6*sqrt(2)/mesh_ref.num_of_el,3);
sin_ref_scale = 3;
tol = 0.01;

% Sinusoidal parameters
edge.coef{1} = -A; % Amplitude
edge.coef{3} = A;
edge.period = N; % Period

% Edge resolution
if A == 0 || N == 0
    edge.shape{1} = 'line'; edge.shape{2} = 'line'; edge.shape{3} = 'line'; edge.shape{4} = 'line';
    line_res = [4 4 4 4];
else
    edge.shape{1} = 'sin'; edge.shape{2} = 'line'; edge.shape{3} = 'sin'; edge.shape{4} = 'line';
    line_res = [sin_ref_scale*round(l/mesh_ref.maxelsize)+4 4 sin_ref_scale*round(l/mesh_ref.maxelsize)+4 4]; 
end
edge.func{1} = @(x) edge.coef{1}*sin(2*pi*edge.period*x/abs(l)) + w/2;
edge.func{3} = @(x) edge.coef{3}*sin(2*pi*edge.period*x/abs(l)) - w/2;

% Create the mesh given geometric parameters
if mesh_ref.exact
    model_3D = createGeometry(l,w,h,line_res,edge,mesh_ref);

%     Percent change by volume of element
%     change_vol = size(model_3D.Mesh.Elements,2)/mesh_ref.num_of_el-1;
%     mesh_ref.maxelsize = (1+change_vol)*model_3D.Mesh.MaxElementSize;
%     model_3D = createGeometry(l,w,h,line_res,edge,mesh_ref);

%     Percent change by elements (works the best but is not perfect)
%     mesh_ref.maxelsize = model_3D.Mesh.MaxElementSize + (size(model_3D.Mesh.Elements,2)-mesh_ref.num_of_el)*model_3D.Mesh.MaxElementSize/mesh_ref.num_of_el;
%     model_3D = createGeometry(l,w,h,line_res,edge,mesh_ref);
%       To do: include Criscione
%       decomposition pipeline (both FtoK2andK3 and galaxy plots)

%     Fminunc implementation (sometimes does not converge)
%     mesh_ref = optimize_mesh(model_3D,mesh_ref,l,w,h,line_res,edge);
    temp = [size(model_3D.Mesh.Elements,2) model_3D.Mesh.MaxElementSize 0];
    i = 1; j = false;
    init = 0.025; % Initial percentage decrement
    breakFlag = false;

    % Step 1 mesh optimization: while loop iterations (until it's within 0.5% of desired mesh elements)
    while (size(model_3D.Mesh.Elements,2) < mesh_ref.num_of_el*0.995 || size(model_3D.Mesh.Elements,2) > mesh_ref.num_of_el*1.005)
        if i == 1
            incr = init;
        else
            if ~j
                if ((temp(i,1) > mesh_ref.num_of_el) == (temp(i-1,1) > mesh_ref.num_of_el)) && (abs(temp(i,1)-mesh_ref.num_of_el) > 0.03*mesh_ref.num_of_el)
                    incr = temp(i,3);
                elseif ~((temp(i,1) > mesh_ref.num_of_el) == (temp(i-1,1) > mesh_ref.num_of_el)) && (abs(temp(i,1)-mesh_ref.num_of_el) > abs(temp(i-1,1)-mesh_ref.num_of_el))
                    incr = temp(i,3)/2;
                else
                    incr = temp(i,3)/2;
                    j = true;
                end
            else
                incr = temp(i,3)/2;
            end
        end
        if size(model_3D.Mesh.Elements,2) > mesh_ref.num_of_el
            mesh_ref.maxelsize = mesh_ref.maxelsize*(1+incr);
        else
            mesh_ref.maxelsize = mesh_ref.maxelsize*(1-incr);
        end
        model_3D = createGeometry(l,w,h,line_res,edge,mesh_ref);
        
        i = i + 1;
        temp = [temp;size(model_3D.Mesh.Elements,2) model_3D.Mesh.MaxElementSize incr];
        if incr < init/2^15
            breakFlag = true;
            break
        end
    end
    
    % Step 2 optimization: Sweep in range centered around best point in Step 1
    if breakFlag
        [~,idx]=min(abs(temp(:,1)-mesh_ref.num_of_el));
        minVal=temp(idx,:);
    else
        minVal=temp(end,:);
    end
    sweep = linspace(minVal(2)-minVal(2)*init/8,minVal(2)+minVal(2)*init/8,30);
    for i = 1:length(sweep)
        mesh_ref.maxelsize = sweep(i);
        model_3D = createGeometry(l,w,h,line_res,edge,mesh_ref);
        temp = [temp;size(model_3D.Mesh.Elements,2) model_3D.Mesh.MaxElementSize init/8];
    end
    [~,idx]=min(abs(temp(:,1)-mesh_ref.num_of_el));
    minVal=temp(idx,:);
    mesh_ref.maxelsize = minVal(2);
    model_3D = createGeometry(l,w,h,line_res,edge,mesh_ref);
else
    model_3D = createGeometry(l,w,h,line_res,edge,mesh_ref);
end
x = model_3D.Mesh.Nodes(1,:)'; y = model_3D.Mesh.Nodes(2,:)'; 
z = model_3D.Mesh.Nodes(3,:)';
TRI = model_3D.Mesh.Elements';

% To view the 3D model, set breakpoint here and type the following commands
% pdeplot3D(model_3D)
% axis off

% Set the material parameters
switch params
    case 'ogden-treloar'
        coef.model = 'Og_3';
        coef.val = [0.4017 1.3 0.00295 5 0.00981 -2 0 0 0];
                      % [mu_1 a_1 mu_2 a_2 mu_3 a_3 k_1 k_2 k_3]
                      % Units are in MPa
                      % See https://solidmechanics.org/text/Chapter3_5/Chapter3_5.htm
                      % Note the correction from the original Treloar data
                      % factors. This is due to the original paper
                      % utilizing the original formulation of the Ogden
                      % model. For the Abaqus formulation, all mu_i values
                      % should be positive (https://polymerfem.com/4-things-you-didnt-know-about-the-ogden-model/).
                      % Correction factor: mu_i =a_i*mu_orig_i/2;
end

% Set all node coordinates as a matrix
Nodes.gen=[x y z];

% Identify boundary surfaces
Surf.x1 = find(and(x<=max(x)+10^-9,x>=max(x)-10^-9));
Surf.x2 = find(and(x<=min(x)+10^-9,x>=min(x)-10^-9));
Surf.z1 = find(and(z<=max(z)+10^-9,z>=max(z)-10^-9));
Surf.z2 = find(and(z<=min(z)+10^-9,z>=min(z)-10^-9));
if A == 0 || N == 0
    Surf.y1 = find(and(y<=max(y)+10^-9,y>=max(y)-10^-9));
    Surf.y2 = find(and(y<=min(y)+10^-9,y>=min(y)-10^-9));
else % For sinusoidal surfaces
    Surf.y1 = []; Surf.y2 = [];
    for j = 1:length(Nodes.gen)
        if ismembertol(y(j),edge.func{1}(x(j)),tol*mesh_ref.maxelsize)
            Surf.y1 = [Surf.y1;j];
        elseif ismembertol(y(j),edge.func{3}(x(j)),tol*mesh_ref.maxelsize)
            Surf.y2 = [Surf.y2;j];
        end
    end
end

% Declare boundary conditions (may be parameterized)
Nodes.bc1 = Surf.y1;
Nodes.bc2 = Surf.y2;
Nodes.presDisp.dir = 'x';
Nodes.presDisp.mag = pres_disp;

% Creating the element set associated with boundary conditions
idx1 = [];
for k = 1:length(Nodes.bc1)
    [row1,~] = find(TRI==Nodes.bc1(k));
    idx1 = [idx1;row1];
end
Elements_Sets{1}.bc1 = unique(idx1);

idx2 = [];
for k = 1:length(Nodes.bc2)
    [row2,~] = find(TRI==Nodes.bc2(k));
    idx2 = [idx2;row2];
end
Elements_Sets{1}.bc2 = unique(idx2);

 % Nodes indices vector for each Elements{i}
for i=1:1:size(TRI,1) 
    Elements{i}=TRI(i,:);              
end

% Element set name
Elements_Sets{1}.Name='Set-1';

% Set element type
switch mesh
    case 'tet'
        Elements_Sets{1}.Elements_Type=elementType;
end

% Elements indices vectors in the element set
Elements_Sets{1}.Elements=1:size(TRI,1);

% Organizes optimization 'sweeps' into subfolders
writeInp(Nodes,Elements,Elements_Sets,fN,coef,abaqus_ver,mesh);
