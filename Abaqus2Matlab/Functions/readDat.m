function U = readDat(fileName,nodeLength)
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

% Finds the locations that lists all nodal displacements (only for final
% increment)
indx = find(strcmp(A,'                                       N O D E   O U T P U T'));
finl_incr = indx(end)+10;
ii = finl_incr;

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