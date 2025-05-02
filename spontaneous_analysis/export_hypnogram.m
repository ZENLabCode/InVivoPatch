clc; clear;
addpath(genpath('X:\OptoLab_v4.1\function'))

files = dir(['X:\1 TIDIS Lab\Tiago\Ephys_PROCESSED\In_vivo\Awake_Sleep\PFC_Processed\*\*\aWake\'...
    '*.H']);
files([files.isdir]) = [];
files = selectionList(fullfile({files.folder},{files.name}));


for fil = 1:numel(files)
    rPath = fileparts(files{fil});
    while ~strcmpi(pwd,rPath)
        cd(rPath)
    end
    f_hypnoCorrect
end


