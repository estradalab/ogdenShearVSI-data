function writeInp(Nodes,Elements,Elements_Sets,Filename,coef,Jobname)

fileID = fopen(['InputFiles\' Filename '_tet.inp'], 'w');

%Heading
fprintf(fileID,'*Heading\n');
fprintf(fileID,['** Job name: ' Jobname ' Model name: Model-1\n']);
fprintf(fileID,['** Generated by: Abaqus/CAE 2021\n*Preprint, echo=NO, model=NO, ' ...
    'history=NO, contact=NO\n**\n** PARTS\n**\n*Part, name=Part-1\n*End Part\n' ...
    '**\n**\n** ASSEMBLY\n**\n*Assembly, name=Assembly\n**\n*Instance, name' ...
    '=Part-1-1, part=Part-1\n']);

%Generate Nodes in Input File

fprintf(fileID,'*Node\n');
[NNode, ND]=size(Nodes.gen);

if ND==2  %2D                               
    for i=1:1:NNode
        fprintf(fileID,[num2str(i) ', ' num2str(Nodes.gen(i,1)) ', ' num2str(Nodes.gen(i,2)) '\n']);
    end

elseif ND==3  %3D
    for i=1:1:NNode
        fprintf(fileID,[num2str(i) ', ' num2str(Nodes.gen(i,1)) ', ' num2str(Nodes.gen(i,2)) ', ' num2str(Nodes.gen(i,3)) '\n']);
    end    
end

%Generate Elements in Input File
fprintf(fileID,strcat('*Element, type=',Elements_Sets{1}.Elements_Type,'\n'));
for i=1:1:length(Elements_Sets)
    for j=1:1:length(Elements_Sets{i}.Elements) %Loop for the elements in the elements set
        IE=Elements_Sets{i}.Elements(j); %Elements indices in elements sets
        NNN=[num2str(IE) ', '];
        for k=1:1:length(Elements{IE})    
            NNN=[NNN num2str(Elements{IE}(k)) ', '];
        end
        NNN=NNN(1:end-2);
        fprintf(fileID,[NNN '\n']);
    end
end

fprintf(fileID,['*Nset, nset=' Elements_Sets{1}.Name ', generate\n1, ' num2str(NNode) ...
    ', 1\n*Elset, elset=' Elements_Sets{1}.Name ', generate\n1, ' num2str(length(Elements_Sets{i}.Elements)) ...
    ', 1\n** Section: Section-1\n*Solid Section, ' ...
    'elset=' Elements_Sets{1}.Name ', controls=EC-1, material=Material-1\n,\n' ...
    '*End Instance\n**\n']);

% Prescribed boundary condition
fprintf(fileID,'*Nset, nset=bc-set-1, instance=Part-1-1\n');
bc(fileID,Nodes.bc1)

% Fixed boundary condition
fprintf(fileID,'*Nset, nset=bc-set-2, instance=Part-1-1\n');
bc(fileID,Nodes.bc2)

fprintf(fileID,['*End Assembly\n**\n** ELEMENT CONTROLS\n**\n*Section' ...
    ' Controls, name=EC-1, hourglass=ENHANCED\n1., 1., 1.\n**\n** MATERIALS\n' ...
    '**\n*Material, name=Material-1\n']);

% Material model coefficients
switch coef.model
    case 'Og_3'
        fprintf(fileID,'*Hyperelastic, n=3, ogden\n');
        for i = 1:floor(length(coef.val)/8)+1
            if i == floor(length(coef.val)/8)+1
                fprintf(fileID,[regexprep(num2str(coef.val(1+(8*(i-1)):end)), ' +', ', ') ',\n']);
            else
                fprintf(fileID,[regexprep(num2str(coef.val(1+(8*(i-1)):8+(8*(i-1)))), ' +', ', ') '\n']);
            end
        end
end

fprintf(fileID,['** ----------------------------------------------------------------\n' ...
    '**\n** STEP: Step-1\n**\n*Step, name=Step-1, nlgeom=YES\n*Static\n' ...
    '1., 1., 1e-05, 1.\n**\n** BOUNDARY CONDITIONS\n**\n']);

% Embed boundary conditions by directions of prescribed displacement

fprintf(fileID,'** Name: BC-1 Type: Displacement/Rotation\n*Boundary\n');
switch Nodes.presDisp.dir
    case 'x'
        fprintf(fileID,['bc-set-1, 1, 1, ' num2str(Nodes.presDisp.mag) '\nbc-set-1, 2, 6\n']);
    case 'y'
        fprintf(fileID,['bc-set-1, 1, 1\nbc-set-1, 2, 2, ' num2str(Nodes.presDisp.mag) ...
            '\nbc-set-1, 3, 6\n']);
    case 'z'
        fprintf(fileID,['bc-set-1, 1, 1\nbc-set-1, 2, 2\nbc-set-1, 3, 3, ' num2str(Nodes.presDisp.mag) ...
            '\nbc-set-1, 4, 6\n']);
end
fprintf(fileID,['** Name: BC-2 Type: Symmetry/Antisymmetry/Encastre\n*Boundary\n' ...
    'bc-set-2, ENCASTRE\n**\n** OUTPUT REQUESTS\n**\n*Restart, write, frequency=0\n' ...
    '**\n** FIELD OUTPUT: F-Output-1\n**\n*Output, field, variable=PRESELECT\n**\n' ...
    '** HISTORY OUTPUT: H-Output-1\n**\n' ...
    '*Output, history, variable=PRESELECT\n*End Step\n']);

fclose(fileID);

end

function bc(fileID,Nodes)
for i = 1:16:length(Nodes)
    if i + 16 >= length(Nodes)
        indx = length(Nodes)-i;
    else
        indx = 15;
    end
    NNN = [];
    for j = 0:indx
        NNN = [NNN num2str(Nodes(i+j))];
        if j ~= indx
            NNN = [NNN ', '];
        else
            NNN = [NNN '\n'];
        end
    end
    fprintf(fileID,NNN);
end
end