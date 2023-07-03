function runInp(fileName)
% Runs the input file through abaqus
system(['abaqus job=' fileName '_test cpus=4']);

while 1
    if any(size(dir([fileName '*.odb']),1)) && any(size(dir([fileName '*.dat']),1)) && ~any(size(dir([fileName '*.simdir']),1))
        break
    else
        pause(1)
    end
end

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