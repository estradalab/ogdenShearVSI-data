function changeInp(fileName)
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

% Finds the job name
for i = 1:length(A)
    C = strsplit(A{i},': ');
    if strcmp(C{1},'** Job name') == 1
        A{i} = ['** Job name: ' fileName '_test Model name: Model-1'];
        break
    end
end

% Finds location to add a set of nodes, RP, and requests an output of
% displacements, U, DG, S12
B = A;
indx = find(strcmp(B,'*End Step')); % Only for 1-step files
A_temp = {B{1:indx-1},'*EL PRINT, position = centroidal','DG',B{indx:end}};
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