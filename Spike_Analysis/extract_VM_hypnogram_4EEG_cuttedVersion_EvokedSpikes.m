% extract_VM_hypnogram
%-------------------------------------------------------------------------
% Extracts VM-data per sleep stage bouts. Additionally saves stage after
% each bout to (bout transition to). Last bout will always have transition
% to NaN.
%
% VM-data might be downsampled, as original recording is huge.
% Exported bout data might be limited to specified duration, and a margin
% to crop start/end of each bout.
%
%
% Thomas Rusterholz, 17 Jun 2021
% Updates:
%  8 Aug 21: exporting indices as well (same sample rate) for to use later
%            with EEGs, EMG, TTL, ... (indices must be corrected to
%            corresponing sample rate, will also be saved in set fsDAT)
%------------------------------------------------------------------------
clc; clear; close all; fclose all;
addpath(genpath('X:\OptoLab_v4.1\function')) %for selectionList.m

%PARAMETER
%----------
%FILES
% Without paths! Searches for files using dir command, so wildcard * works
%  - fileHYP: hypnogram files (must be unique per read path)
%  - fileBCI: OpenBCI txt-file to find start date/time (must be unique)
%  - fileVM : mat-files exported from abf-files. Performs a loop across all files found from read path.
fileHYP = 'Hypnogram.mat';
fileBCI = 'OpenBCI-RAW-*.txt';
fileVM  = '*_IN0.mat';
fileEMG = 'EMG.mat';
%READ PATHS
% Automatically searche read paths using dir command. Opens selection list
% with paths found.
% Best is to search for paths having one of the upper files
startSearchPath = 'X:\1 TIDIS Lab\Tiago\Ephys_PROCESSED\In_vivo\Awake_Sleep\RC_Processed\20211119_TC_P046_Day1\Evoked_Spikes';

%OPTIONS
%sampling rates
%  - fsDAT: in [Hz], VM data to resample to. Original is 10^5 Hz
%  - fsHYP: in [Hz], scored hypnogram (to upsample to fsDAT)
fsDAT = 5000;
fsHYP = 1;
%limit exported stage bouts
%  - margin  : in [s], crops margin from start AND end of bouts
%              (in total crops 2*margin seconds)
%  - durLimit: in [s], bouts smaller will be excluded
%              PS: durLimit inclusive margin, margin cut afterwards.
%                  --> durLimit must be > 2*margin !!!
margin   = 1;
durLimit = 2*margin + 1; %[s]
%stages to export based on hypnogram
Stages = {... {number, label}
    1, 'Wake';...
    2, 'NREM';...
    3, 'REM';...
    };

%MAIN SCRIPT
%------------
scriptName = mfilename;
fprintf('%s\n%s\n',scriptName,repmat('-',size(scriptName)))

%FIND & SELECT FILES
tmp = dir(fullfile(startSearchPath,'**',fileVM));
filesVM = selectionList(fullfile({tmp.folder},{tmp.name}));
noFIL = numel(filesVM);

%info to append
info.info = char({...
    sprintf('Exported by %s.m, %s',scriptName,date);...
    '';...
    'Exports VM data for stage bouts.';...
    'Date cropped by margin & bout duration limited.';...
    'Also exports stage after each bout.';...
    '';...
    'Info fields:';...
    ' - margin   : Margin in [s], bouts cropped at start & end.';...
    '              Total crop is 2 times margin!';
    ' - minDur   : Minimal bout duration of cropped bouts';...
    '              Shorter ones excluded!';...
    ' - indBouts : EEG stage bouts start/end indices';...
    '              Margin alredy excluded!';...
    '              In sample rate fs!';...
    });
info.margin = margin;
info.minDur = durLimit-2*margin;

%FILE LOOP
nnFIL = numel(num2str(noFIL));
indent = blanks(2*nnFIL+2);
for fil = 1:numel(filesVM)
    tic
    [rPath,rFile] = fileparts(filesVM{fil});
    fprintf('%*i/%i: %s\n',nnFIL,fil,noFIL,rFile)    
    
    %LOAD DATA
    %hypnogram
    fname = fullfile(rPath,rFile,fileHYP);
    if exist(fname,'file')~=2
        fprintf(2,'%s %s NOT found\n',indent,fileHYP)
        continue
    end
    tmp = load(fname);
    hypnogram = tmp.Hypnogram(:);
    if fsDAT~=fsHYP %upsample
        hypnogram = repmat(hypnogram',fsDAT/fsHYP,1);
        hypnogram = hypnogram(:);
    end
    hypnogram(end+1) = NaN;
    %vm
    fname = fullfile(rPath,[rFile,'.mat']);
    tmp = load(fname);
    if fsDAT==tmp.SampRate
        dat = tmp.resampled_data_mV(:);
    else %resample data
        dat = resample(tmp.resampled_data_mV(:),fsDAT,tmp.SampRate);
    end

    %BOUTS
    clear res;
    res.info = info;
    res.fs = fsDAT;
    for sta = 1:size(Stages,1)
        [stage,label] = Stages{sta,:};
        indSTA = find(...
            [hypnogram;NaN]==stage & ...
            [NaN;hypnogram]~=stage);
        indEND = find(...
            [hypnogram;NaN]~=stage & ...
            [NaN;hypnogram]==stage)-1; %-1, end belongs to bout
        dur = (indEND-indSTA+1)/fsDAT;
        
        %test plot
        if false
            figure
            t = 1:numel(hypnogram);
            plot(t,hypnogram,'b','marker','.'); hold on
            plot(t(indSTA),hypnogram(indSTA),'go')
            plot(t(indEND),hypnogram(indEND),'ro')
            if ~isempty(indSTA)
                legend('Hypnogram','Bout Start','Bout End')
            end
            title(sprintf('%s Bouts',label))
            ylim([0.8 3.2]); xlim([0,t(end)]); zoom xon;
            %durations
            ind = dur<durLimit;
            cols = {'r','k'}; vpos = {'top','bottom'};
            for k=1:2
                x = t(round((indSTA(ind)+indEND(ind))/2));
                y = ones(size(x))*stage;
                text(x,y,cellfun(@(x)sprintf('%g s',x),...
                    num2cell(dur(ind)'),'uniformoutput',false),...
                    'color',cols{k},...
                    'horizontalalignment','center',...
                    'verticalalignment',vpos{k})
                ind = ~ind;
            end
            drawnow
        end
        
        %duration limit
        ind = dur<durLimit;
        indSTA(ind) = [];
        indEND(ind) = [];
        
        %APPEND BOUT DATA
        N = numel(indSTA);
        res.(label) = cell(N,1);
        res.info.indBouts.(label) = NaN(N,2);
        indMAR = margin*fsDAT;
        for c = 1:N
            ind = indSTA(c)+indMAR:indEND(c)-indMAR;
            res.(label){c} = dat(ind);
            res.info.indBouts.(label)(c,:) = ind([1,end]);
        end
        %stage after bouts
        res.(sprintf('stageAfter%s',label)) = hypnogram(indEND+1);
    end
    
    %test plot
    if false
        figure
        tmp = res.NREM;
        for q=1:numel(tmp); plot(tmp{q}); pause(2);end
        return
    end
    
    %SAVE
    tmp = regexp(rFile, '_','split');
    sFile = sprintf('%s_vmStage.mat',tmp{1});
    save(fullfile(rPath,rFile,sFile),'-struct','res');
    fprintf('%s Saved %s\n',indent,sFile)
end %file loop
fprintf('%s ',indent); toc