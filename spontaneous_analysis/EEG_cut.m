% EEG_cut
%------------------------------------------------------------------------
% Cut EEGs to VM data recordings. So VM and EEGs must be synchronized
% first.
% Data will be saved to sub-folders, using VM filenames as sub-folder
% names
%
%
% Thomas Rusterholz, 11 Apr 2022
%------------------------------------------------------------------------

clc; clear; close all
addpath(genpath('X:\OptoLab_v4.1\function')) %for selectionList.m


%PARAMETERS
%----------
%FILES NAMES (without paths)
%data files to cut (*.mat)
%  - Extension not needed (will be replaced with .mat anyway)
%  - Also cut TTL for to re-check synchronization!
Files.EEG = {'EEG1','EEG2','EMG','TTL'}; %cell, data to cut
%vm data files (e.g. *INO.mat for to find VM data files)
Files.VM  = '*IN0.mat'; % e.g. *INO.mat for to find VM data files
%openBCI file (for recording start data from filename)
% PS: reads start date from filename and start time within file
Files.BCI = 'OpenBCI-RAW*.txt'; %for recording start date (*.txt)
dateFormat = 'OpenBCI-RAW-yyyy-mm-dd_HH-MM-SS';

%READ PATHS
tmp = dir(fullfile('X:\1 TIDIS Lab\Tiago\Ephys_PROCESSED\In_vivo\Awake_Sleep\RC_Processed\20211119_TC_P046_Day1\Evoked_Spikes',...
    Files.VM));
rPaths = selectionList(unique({tmp.folder})); %unique paths!


%MAIN SCRIPT
%-----------
scriptName = mfilename;
fprintf('%s\n%s\n',scriptName,repmat('-',size(scriptName)))

%INIT
if ~iscell(Files.EEG)
    Files.EEG = {Files.EEG};
end
%number of ...
noPAT  = numel(rPaths);
noEEG = numel(Files.EEG);
%remove extensions
for k = 1:noEEG
    [~,Files.EEG{k}] = fileparts(Files.EEG{k});
end

%PATH LOOP
nnPAT  = numel(num2str(noPAT));
indent = blanks(2*nnPAT+2);
for pat = 1:noPAT
    rPath = rPaths{pat};
    fprintf('%*i/%i: %s\n',nnPAT,pat,noPAT,rPath)
    
    %VM FILES
    filesVM = dir(fullfile(rPath,Files.VM));
    if isempty(filesVM)
        fprintf(2,'%s No file %s found\n',indent,filesVM)
        continue
    end
    filesVM = {filesVM.name};
    
    %GET START DATE/TIME TTL
    tmp = dir(fullfile(rPath,Files.BCI));
    fname = fullfile(tmp.folder,tmp.name);
    [~,fileBCI] = fileparts(fname); %cut extension
    fid = fopen(fname,'r');
    str = strtrim(fgetl(fid));
    while str(1)=='%' || isempty(str)
        str = strtrim(fgetl(fid));
    end
    fclose(fid);
    str = regexp(str,' ','split');
    str = regexp(str{end},':','split');
    %create datevec
    [~,rFile] = fileparts(fileBCI);
    startDate = datevec(rFile,dateFormat);
    startDate(4) = str2double(str{1});
    startDate(5) = str2double(str{2});
    startDate(6) = str2double(str{3});
    
    %LOAD EEGs
    fprintf('%s LOAD DATA\n',indent)
    EEGS = cell(noEEG,1);
    clear SampRate
    for eeg = 1:noEEG
        [~,rFile,rExt] = fileparts(Files.EEG{eeg});
        fprintf('%s   %s\n',indent,rFile)
        fname = fullfile(rPath,[rFile,'.mat']);
        if exist(fname,'file')~=2
            fprintf('\b - [\bdoes NOT exist!]\b\n')
            continue
        end
        %load data
        tmp = load(fname);
        if ~strcmpi(fileBCI,tmp.id)
            error('Something is wrong')
        end
        if ~exist('SampRate','var') %do only once per path
            SampRate  = tmp.SampRate; 
        elseif ~isequal(SampRate,tmp.SampRate)
            error('Something is wrong')
        end
        %append data
        EEGS{eeg} = tmp.resampled_data_mV;
    end
    noSAM = min(cellfun(@numel,EEGS)); %number of samples
    
    %CUT BY SYNCHRONIZATION DATE
    fprintf('%s CUT DATA\n',indent)
    for k = 1:numel(filesVM)
        [~,rFile,rExt] = fileparts(filesVM{k});
        fprintf('%s   %s\n',indent,rFile)
        tmp = matfile(fullfile(rPath,[rFile,rExt]));
        try
            startDateShifted = tmp.startDateShifted;
        catch
            fprintf('\n - [\bnot vet synchronized!]\b\n')
            continue
        end
        %save path
        sPath = fullfile(rPath,rFile);
        
        %CUT INDEX (start & end)
        dur    = prod(size(tmp,'resampled_data_mV'))/tmp.SampRate;
        indSTA = round(etime(startDateShifted,startDate)*SampRate);
        if indSTA==0 %round up instead down (probably never happens)
            indSTA = 1;
        end
        indEND = indSTA + ceil(dur*SampRate) - 1;
        if indEND>noSAM
            if indEND-2<=noSAM %ok because might be rounded up twice
                indEND = noSAM;
            else
                fprintf(2,'%s     VM outside EEG\n',indent)
                continue
            end
        end
       
        %cut channels
        for eeg = 1:noEEG
            data = EEGS{eeg};
            if isempty(data)
                continue
            end
            %data to save
            info.fileVM  = rFile;
            info.fileBCI = fileBCI;
            info.cutIndex = [indSTA,indEND];
            channel = Files.EEG{eeg};
            resampled_data_mV = data(indSTA:indEND);
            %save
            sname = fullfile(sPath,[channel,'.mat']);
            if true %to test, prevent saving anything
                if ~exist(sPath,'dir')
                    mkdir(sPath)
                end
                save(sname,'info','SampRate','channel','resampled_data_mV')
            end
            fprintf('%s     Saved: ...%s\n',indent,...
                strrep(sname,rPath,''))
        end %files EEG
    end %files VM
end %path loop

