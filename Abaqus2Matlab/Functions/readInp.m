function [X,Y,Z,ElNode] = readInp(fileName)
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

% Finds the critical indices that list node reference positions and
% connectivity matrix
temp = find(strcmp(A,'*Node'));
indx_nd = temp(1);
for i = 1:length(A)
    C = strsplit(A{i},', ');
    if strcmp(C{1},'*Element') == 1
        indx_el = i;
        break
    end
end

% Stores node references positions into matrices, X, Y, and Z
for i = indx_nd+1:length(A)
    C = strsplit(A{i},', ');
    X(i-indx_nd,1)=str2double(C{2}); Y(i-indx_nd,1)=str2double(C{3}); Z(i-indx_nd,1)=str2double(C{4});
    D = strsplit(A{i+1},', ');
    if isempty(str2double(D{1})) || str2double(D{1})-str2double(C{1}) ~= 1
        break
    end
end

% Stores connectivity matrix into matrix, ElNode
ElNode = [];
for i = indx_el+1:length(A)
    C = strsplit(A{i},', ');
    ElNode = [ElNode; str2double(C)];
    D = strsplit(A{i+1},', ');
    if isempty(str2double(D{1})) || str2double(D{1})-str2double(C{1}) ~= 1
        break
    end
end
ElNode(:,1) = [];