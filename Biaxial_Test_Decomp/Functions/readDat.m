function [U,F] = readDat(fileName,nodeLength)
% Stores .dat file into MatLab cell array, A
fid = fopen([fileName '_test.dat'],'r');
i=1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);

% Finds the locations that lists all element displacement gradients
idx_search = [];
for j = 1:length(A)-1
    splitline = strsplit(A{j});
    if ~isempty(find(strcmp(splitline, 'INCREMENT'),1)) && ~isempty(find(strcmp(splitline, 'SUMMARY'),1))
        idx_search = [idx_search j];
    end
end
indx = max(idx_search);
ii = indx(end)+19;

% Store displacement gradient into F cell matrix
while length(strsplit(A{ii})) > 1
    temp = strsplit(A{ii});
    F{1,1}(str2double(temp{2})) = str2double(temp{3});
    F{2,2}(str2double(temp{2})) = str2double(temp{4});
    F{3,3}(str2double(temp{2})) = str2double(temp{5});
    F{1,2}(str2double(temp{2})) = str2double(temp{6});
    F{1,3}(str2double(temp{2})) = str2double(temp{7});
    F{2,3}(str2double(temp{2})) = str2double(temp{8});
    F{2,1}(str2double(temp{2})) = str2double(temp{9});
    F{3,1}(str2double(temp{2})) = str2double(temp{10});
    F{3,2}(str2double(temp{2})) = str2double(temp{11});
    ii = ii + 1;
end

% Finds the locations that lists all nodal displacements (only for final
% increment)
for j = indx:length(A)-1
    if ~isempty(find(strcmp(A{j},'                                       N O D E   O U T P U T'),1))
        ii = j+10;
        break
    end
end

% Stores displacements into temporary cell matrix
while length(strsplit(A{ii})) > 1
    temp = strsplit(A{ii});
    B{str2double(temp{2})} = [str2double(temp{3}) str2double(temp{4}) str2double(temp{5})];
    ii = ii + 1;
end

% Stores nodal displacements into matrix, U, of which all "empty" cells are
% those that are encastred (U = [0 0 0])
for ii = 1:nodeLength
    try
        U(ii,:) = B{ii};
    catch
        U(ii,:) = zeros(1,3);
    end
end
