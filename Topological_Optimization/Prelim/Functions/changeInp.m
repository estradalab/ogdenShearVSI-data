function changeInp(fileName,nodeLength,elLength)
% Stores .inp file into MatLab cell array, A
fid = fopen([fileName '.inp'],'r');
i=1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);

% Finds name of "Instance"
for i = 1:length(A)
    C = strsplit(A{i},', ');
    if strcmp(C{1},'*Instance') == 1
        D = strsplit(C{2},'name=');
        instance_name = D{2};
        break
    end
end

% Finds the job name
for i = 1:length(A)
    C = strsplit(A{i},': ');
    if strcmp(C{1},'** Job name') == 1
        A{i} = ['** Job name: ' fileName '_test Model name: Model-1'];
        break
    end
end

% Finds location to add a set of nodes, RP, and requests an output of
% displacements, U1, U2, U3
indx = find(strcmp(A,'*End Instance'));
A_temp = {A{1:indx+1},['*Nset, nset=RP, instance=' instance_name ', generate'],['     1,  ' num2str(nodeLength) ',      1']...
    ,'*Elset, elset=RP, instance=Part-1-1, generate',['     1,  ' num2str(elLength) ',      1'],A{indx+2:end}};
B = A_temp;
indx = find(strcmp(B,'*End Step')); % Only for 1-step files
A_temp = {B{1:indx-1},'*EL PRINT','DG',B{indx:end}};
B = A_temp;
A_temp = {B{1:indx-1},'*NODE PRINT','U',B{indx:end}};
B = A_temp;

% Overwrites the additions into a new file
fid = fopen([fileName '_test.inp'],'w');
for i = 1:numel(B)
    if B{i+1} == - 1
        fprintf(fid,'%s', B{i});
        break
    else
        fprintf(fid,'%s\n',B{i});
    end
end