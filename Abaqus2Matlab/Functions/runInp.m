function runInp(fileName,jobName)
% Runs the input file through abaqus
system(['abaqus job=' jobName '_test cpus=4 interactive']);

% Makes a new directory to store ammended .inp file, .odb output database,
% and .dat data
if ~exist(['Data\' fileName], 'dir')
    mkdir(['Data\' fileName])
end

% Moves files into respective directory, and clears clutter
movefile([jobName '_test.dat'],['Data\' fileName],'f');
movefile([jobName '_test.odb'],['Data\' fileName],'f');
fclose('all');
movefile([jobName '_test.inp'],['Data\' fileName],'f');
delete *_test.*