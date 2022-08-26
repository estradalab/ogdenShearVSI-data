function runInp(fileName)
% Runs the input file through abaqus
system(['abaqus job=' fileName '_test cpus=4 interactive']);

% Makes a new directory to store ammended .inp file, .odb output database,
% and .dat data
if ~exist(['Data\' fileName], 'dir')
    mkdir(['Data\' fileName])
end

% Moves files into respective directory, and clears clutter
movefile([fileName '_test.dat'],['Data\' fileName],'f');
movefile([fileName '_test.odb'],['Data\' fileName],'f');
fclose('all');
movefile([fileName '_test.inp'],['Data\' fileName],'f');
delete *_test.*