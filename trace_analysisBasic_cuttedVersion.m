% trace_analysisBasic
%------------------------------------------------------------------------
% Calulates basic Vm statistic, histogram of Vm and cumulative sum of
% Vm-min(Vm), NaN values excluded.
%
% Histogram and cumsum will be normalized:
%  - histograms given in percentage of the amount of data points
%    --> equal weight of traces by averaging histograms across traces
%  - cumulative sum is normalized by the sampling rate
%    --> independent results having variable sampling rates
%    --> unit mV/Hz = mV*s
%   PS: cumulative sum of Vm-min(Vm)
%
% NEW, Also exports:
%  - Also includes EMG
%  - From extracted video data (if available)
%    (Reads TTL for video offset)
%    - Pupile, normalized by max value
%    - Motion, mean absolute value across all lines (500)
%

% Note: Traces to be analyzed can be selected. Excluded traces will not be
%       exported!
%
%
% Thomas Rusterholz, 30 Jun 2021
%------------------------------------------------------------------------
clc; clear; close all
addpath(genpath('Z:\OptoLab_v4.1\function'))

%PARAMETERS
%----------
%FILES (*vmStage !!!)
% - Auto searching for files using command dir
%   E.g.: files.vmStage = C:\data\**\*_vmStage.mat;
%   (uses TTL & VM files from the same folder)
% - Auto matching video files using last sub-path from vmStage file and
%   recording timestamp from *IN0.mat or *IN1.mat
% - Opens selection list to select vmStage files from files found
%   Option mustHaveVideo: if true, only displays vmStage file having a
%                         video file

files.vmStage = ['X:\1 TIDIS Lab\Tiago\Ephys_PROCESSED\In_vivo\Awake_Sleep\PFC_Processed\Pharmacology\NMDAR_Blockers\*\*\*_vmStage_10sec.mat'];
files.video   = ['X:\1 TIDIS Lab\Tiago\Pupil_Tracking_RAW\In_vivo_Patch\20210923_TC_P032_Day1\Image 2021_09_23_120043,069_processed_proc.mat'];
mustHaveVideo = false;

%STAGES
Stages = {...{Stage, corresponding number, traces to exclude}
    'Wake', 1, [];...
    'NREM', 2, [];...:
    'REM' , 3, [];...
    };

%EMG QUANTIFICATION (total power of frequency band)
f = 10:0.1:100; %frequency vector for pwelch

%HISTOGRAM (does not plot, only to save data)
edges = -80:-20; %in [mV]

%MAIN SCRIPT
%-----------
scriptName = mfilename;
fprintf('%s\n%s\n',scriptName,repmat('-',size(scriptName)))

%FIND FILES
%vmStage
tmp = dir(files.vmStage);
tmp([tmp.isdir]) = []; %must be file, even when miserable set ;-)
if isempty(tmp)
    fprintf(2,'No files found for: %s\n',files.vmStage)
    return
end
Files = cell2struct(fullfile({tmp.folder},{tmp.name}),'vmStage');
%video files
tmp = dir(files.video);
if isempty(tmp) && mustHaveVideo
    fprintf(2,'No files found for: %s\n',files.video)
    return
end
filesVID = fullfile({tmp.folder},{tmp.name})';

%MATCH FILES
for fil = 1:numel(Files)
    [rPath,rFile] = fileparts(Files(fil).vmStage);
    Files(fil).TTL = fullfile(rPath,'TTL.mat');
    Files(fil).EMG = fullfile(rPath,'EMG.mat');
    Files(fil).HYP = fullfile(rPath,'Hypnogram.mat');
    %label
    tmp = regexp(rFile,'_','split');
    label = tmp{1};
    %find * _IN0.mat file path
    mPath = rPath; %main path having BCI-file
    file  = sprintf('%s_IN0.mat',label);
    while exist(fullfile(mPath,file),'file')~=2
        mPath = fileparts(mPath);
        %break if no sub-path anymore
        if numel(mPath)<=3 %e.g. Z:\
            error('grrhh')
        end
    end
    [~,subPath] = fileparts(mPath);
    
    %get start time from IN-file
    try
        Files(fil).VM = fullfile(mPath,sprintf('%s_IN0.mat',label));
        tmp = load(Files(fil).VM,'info');
        T0 = tmp.info.recStart;
    catch
        Files(fil).VM = fullfile(mPath,sprintf('%s_IN1.mat',label));
        tmp = load(Files(fil).VM,'info');
        T0 = tmp.info.recStart;
    end
    
    %match video files
    ind = contains(filesVID,[filesep,subPath,filesep]);
    fileVID = filesVID(ind);
    if isempty(fileVID)
        Files(fil).video = {};
        Files(fil).led   = {};
    else
        %match time
        dt = NaN(numel(fileVID),1);
        for k = 1:numel(fileVID)
            [~,str] = fileparts(fileVID{k});
            str = regexp(str,'[0-9_]*','match');
            T1 = datevec(str{1},'yyyy_mm_dd_HHMMSS');
            dt(k) = etime(T1,T0);
        end
        [~,ind] = min(abs(dt));
        fileVID = fileVID{ind};
        if dt(ind)>5*60 %that's probably to long time delay
            warning('Grrhhh, do something')
        end
        %append
        [rPath,rFile] = fileparts(fileVID);
        Files(fil).video = fileVID;
        Files(fil).led = fullfile(rPath,[rFile,'_led.mat']);
    end
    
    %append label
    Files(fil).label = label;
end
%ckeck for unique video files
filesVID = {Files.video};
filesVID(cellfun(@isempty,filesVID)) = [];
if numel(filesVID)~=numel(unique(lower(filesVID)))
    error('Video files did not uniquely match with vmStage files')
end
clear filesVID;

%SELECT FILES
if mustHaveVideo
    ind = cellfun(@isempty,{Files.video});
    Files(ind) = [];
end
tmp = selectionList({Files.vmStage});
if isempty(tmp)
    fprintf(2,'No File Selected!\n')
    return
end
Files(~ismember({Files.vmStage},tmp)) = [];

%INIT
centers = edges(1:end-1) + diff(edges)/2;
noFIL = numel(Files);
noSTA = size(Stages,1);
%test figure for video delay
cols = ceil(sqrt(noFIL)); rows = ceil(noFIL/cols);
hf = figure('Name','Estimated Video Delay');
ha = fig_createAxes(hf,[rows,cols],[50 20 20],[50 50 70],'pixel');
set(ha,'unit','normalized')
set(hf,'CurrentAxes',ha(1,1)); %for super title, needs string distance
title(' '); fig_superLabel(ha,'title','TTL')
set(ha(noFIL+1:end),'visible','off')

%LOOP FILES
nnFIL = numel(num2str(noFIL));
indent = blanks(2*nnFIL+2);
for fil = 1:noFIL
    files = Files(fil);
    fprintf('%*i/%i: %s\n',nnFIL,fil,noFIL,files.label)
    fprintf('%s Path data  : %s\n',indent,fileparts(files.vmStage))
    if isempty(files.video)
        fprintf('%s [\bNo video file!]\b\n',indent)
    else
        fprintf('%s Path video : %s\n',indent,fileparts(files.video))
    end
    fprintf('%s Load Data\n',indent)
    
    %LOAD vmStage (data already per bout, inclusive indices)
    fname = files.vmStage;
    [sPath,rFile,rExt] = fileparts(fname); %read path = save path
    fprintf('%s  - %s\n',indent,[rFile,rExt])
    tmp = load(fname);
    clear VM
    VM.fs       = tmp.fs;
    VM.indBouts = tmp.info.indBouts;
    VM.Data     = tmp.no_spike; %per stage bout already
    VM.minDur   = tmp.info.minDur;
    
    %LOAD EMG (can use indBouts directly, same as for EEG)
    fname = files.EMG;
    [~,rFile,rExt] = fileparts(fname);
    fprintf('%s  - %s\n',indent,[rFile,rExt])
    tmp = load(fname);
    clear EMG
    EMG.fs   = tmp.SampRate;
    EMG.data = tmp.resampled_data_mV(:);
    
    %LOAD HYPNOGRAM (for a check only)
    fname = files.HYP;
    [~,rFile,rExt] = fileparts(fname);
    fprintf('%s  - %s (for checking only)\n',indent,[rFile,rExt])
    tmp = load(fname);
    hypnogram = tmp.Hypnogram(:);
    %upsample hypnogram to EMG if needed
    fac  = EMG.fs;
    fac2 = floor(numel(EMG.data)/numel(hypnogram));
    if false && fac~=fac2 %set to false to ignore
        prompt = sprintf(['EMG sampling rate is %g Hz\n',...
            'However, estimated sampling rate based on data size is %g Hz\n',...
            '(Maybe not full scored to the end)\n',...
            'Upsampling factor is %g if hypnogram is scored for 1-s epochs,\n',...
            'else change manually (remaining will be padded with NaNs)'],...
            fac,fac2,fac);
        answer = inputdlg({prompt},'Set Hypnogram Upsamlping Factor',...
            [1 65],{num2str(fac)});
        if isempty(answer)
            error('Hypnogram upsampling factor not set!')
        end
        fac = str2double(answer);
    end
    if fac>1
        hypnogram = repmat(hypnogram',fac,1);
        hypnogram = hypnogram(:);
        fprintf('%s    Upsampled by factor %i\n',indent,fac)
    end
    hypnogram(numel(hypnogram)+1:numel(EMG.data)) = NaN;
    
    %EMG INDEX (relative to vm bouts)
    n = numel(EMG.data);
    for k = 1:noSTA
        [stage,num,~] = Stages{k,:};
         ind = round(VM.indBouts.(stage)/VM.fs*EMG.fs);
        if any(ind(:)<0 | ind(:)>n+1) || any(ind(1:end-1,2)>ind(2:end,1))
            error('uups'); %something is wrong
        end
        ind(ind==0)   = 1; %might happen because rounding
        ind(ind==n+1) = n; %might happen because rounding
        %hypnogram check
        for q = 1:size(ind,1)
            if ~all(hypnogram(ind(q,1):ind(q,2))==num)
                error('uups'); %must be an insex error
            end
        end
        %append
        EMG.ind.(stage) = ind;
    end
    
    %LOAD VIDEO DATA (if available)
    figure(hf); %activate this one
    set(hf,'CurrentAxes',ha(fil));
    title(files.label);
    if isempty(files.video)
        set(gca,'xtick',[],'ytick',[],'xlim',[0 1],'ylim',[0 1])
        text(.5,.5,'No Video Data',...
            'horizontalalignment','center','verticalalignment','middle')
    else
        %TTL
        fname = files.TTL;
        [~,rFile,rExt] = fileparts(fname);
        fprintf('%s  - %s\n',indent,[rFile,rExt])
        tmp = load(fname);
        clear TTL
        TTL.fs   = tmp.SampRate; %same name in all data structures
        TTL.data = tmp.resampled_data_mV;
        TTL.t    = (1:numel(TTL.data))/TTL.fs;
        %video delay (synchronize video data to vm)
        lim = mean(TTL.data)+5*std(TTL.data);
        ind = find(abs(TTL.data)>lim,1,'first')-1;
        dt  = TTL.t(ind);
        %test plot
        plot(TTL.t(1:ind),TTL.data(1:ind),'r'); hold on
        plot(TTL.t(ind+1:end),TTL.data(ind+1:end),'b');
        plot([0,TTL.t(end)],[lim,lim],'k:')
        set(gca,'xlim',[0,min([TTL.t(end),10*dt])])
        text(dt,TTL.data(ind),sprintf('Video Delay: %g s',dt),...
            'horizontalalignment','left','verticalalignment','bottom')
        text(max(xlim),lim,'Limit: Mean + 5*STD ',...
            'horizontalalignment','right','verticalalignment','bottom')
        title(files.label); xlabel('Time [s]'); zoom xon
        drawnow
        
        %VIDEO ()
        fname = files.video;
        [rPath,rFile,rExt] = fileparts(fname);
        fprintf('%s  - %s\n',indent,[rFile,rExt])
        tmp = load(fname,'pupil','filenames','motSVD_1');
        %movie (needed ???)
        fname = squeeze([tmp.filenames])';
        [~,rFile,rExt] = fileparts(fname);
        if exist(fname,'file')~=2 %if was moved to another folder
            fname = fullfile(rPath,[rFile,rExt]);
        end
        fprintf('%s  - %s\n',indent,[rFile,rExt])
        %append
        clear VID;
        VID.movie  = VideoReader(fname);
        VID.fs     = VID.movie.FrameRate;
        VID.pupil  = tmp.pupil{1};
        VID.pupil  = double(VID.pupil.area_smooth(:));
        VID.pupil  = VID.pupil/max(VID.pupil);
        VID.motion = double(mean(abs(tmp.motSVD_1),2));
        VID.t      = (1:numel(VID.pupil))/VID.fs;
        
        %% VIDEO INDEX
        n = numel(VID.t);
        for k = 1:noSTA
            [stage,num] = Stages{k,1:2};
            ind = VM.indBouts.(stage);
            ind = round((ind/VM.fs-dt)*VID.fs); % minus dt !!!
            VID.ind.(stage) = ind;
            
            %             %correction
            %             ind(ind<1) = 1;
            %             ind(ind>n) = n;
            %             dur = (diff(ind,[],2)+1)/VID.fs;
            %             ind(dur<VM.minDur,:) = NaN;
            %             %append
            %             VID.ind.(stage) = ind;
            
            %optional calc
            % t = VID.t + (VM.indEEG(1)-1)/VM.fs + dt;
            % indB = VM.indBouts.(stage);
            % IND = NaN(size(indB));
            % for q = 1:size(IND,1)
            %     [~,IND(q,1)] = min(abs(t - indB(q,1)/VM.fs));
            %     [~,IND(q,2)] = min(abs(t - indB(q,2)/VM.fs));
            % end
            % VID.ind.(stage) = IND;
        end
        clear stages %just in case, should not be used, use Stages
    end
    
    %% INIT
    clear RES
    RES.info = sprintf('Exported by %s.m, %s',scriptName,date);
    RES.label = files.label;
    
    %LOOP STAGES
    for sta = 1:noSTA
        [stage,num,traceExcl] = Stages{sta,:};
        Data = VM.Data.(stage); %analysis on data without spikes
        if isempty(Data)
            RES.(stage) = [];
            fprintf('%s [\bNo %s trace!]\b\n',indent,stage)
            continue
        end
        
        %TRACES (exclusions)
        Traces = 1:numel(Data);
        Traces(ismember(Traces,traceExcl)) = [];
        noTRA = numel(Traces);

        %LOOP TRACE
        clear Res
        for tra = 1:noTRA
            trace = Traces(tra);
            data  = Data(trace).data(:);
            t     = (1:numel(data))'/VM.fs;
            
            %INIT
            Res.VM(tra).trace = trace;
            Res.VM(tra).t     = t;
            Res.VM(tra).data  = data;
            %remove NaNs from rest of analysis
            data(isnan(data)) = [];
            dataS = data-min(data);
            
            %STATISTIC VALUES
            Res.VM(tra).mean = mean(data);
            Res.VM(tra).std  = std(data);
            Res.VM(tra).prctile15 = prctile(data,15);
            Res.VM(tra).prctile95 = prctile(data,90);
            Res.VM(tra).prctile95_subtracted = prctile(dataS,95);
            
            %HISTOGRAM (normalized)
            counts = histcounts(data,edges) / numel(data)*100;
            Res.VM(tra).histX = centers(:);
            Res.VM(tra).histY = counts(:);
            Res.VM(tra).histLabelX = '[mV]';
            Res.VM(tra).histLabelY = 'Counts [%]';
            
            %CUMSUM (normalized)
            Res.VM(tra).cumsumX = (1:numel(dataS))'/VM.fs;
            Res.VM(tra).cumsumY = cumsum(dataS)/VM.fs;
            % MOD - Tiago Campelo (31/01/2022)
            Res.VM(tra).maxcumsum = max (Res.VM(tra).cumsumY);
            % END_MOD
            Res.VM(tra).cumsumLabelX = 'Time [s]';
            Res.VM(tra).cumsumLabelY = '[mV\cdots]';
            
            %VIDEO VM
            if isempty(files.video)
                Res.pupil  = [];
                Res.motion = [];
            else
                ind = VID.ind.(stage)(trace,1):VID.ind.(stage)(trace,2);
                t = (1:numel(ind))'/VID.fs;
                N = numel(VID.t);
                ind2 = ind>=1 & ind<=N;
                fields = {'pupil','motion'};
                for k = 1:numel(fields)
                    field = fields{k};
                    data = NaN(size(t));
                    data(ind2) = VID.(field)(ind(ind2));
                    Res.(field)(tra).trace = trace;
                    Res.(field)(tra).t     = t;
                    Res.(field)(tra).data  = data;
                    Res.(field)(tra).mean  = nanmean(data);
                    Res.(field)(tra).std   = nanstd(data);
                    %durations trace & available data
                    % PS: durTOT might be slightly different from dur of
                    %     EMG, as index was rounded to fs for video
                    Res.(field)(tra).durTOT = t(end);
                    Res.(field)(tra).durDAT = sum(ind2)/VID.fs;
                end
            end
            
            %EMG
            ind  = EMG.ind.(stage)(trace,1):EMG.ind.(stage)(trace,2);
            data = EMG.data(ind);
            Res.EMG(tra).trace = trace;
            Res.EMG(tra).t     = (1:numel(ind))'/EMG.fs;
            Res.EMG(tra).data  = data;
            pxx = pwelch(detrend(data),[],0,f,EMG.fs);
            Res.EMG(tra).pxx = pxx;
            Res.EMG(tra).f   = f;
            Res.EMG(tra).totalPower = sum(pxx)*mean(diff(f));
            %not usable
            % Res.EMG(tra).mean  = mean(data);
            % Res.EMG(tra).std   = std(data);
        end %loop trace
        
        %APPEND
        if exist('Res','var')
            RES.(stage) = Res;
        else
            RES.(stage) = [];
        end
    end %loop stages
    
    %% PLOT
    %clc; close all
    figure;
    ha1 = NaN(4,noSTA);
    for sta = 1:noSTA
        [stage,num,~] = Stages{sta,:};
        res = RES.(stage);
        if isempty(res)
            continue
        end
        
        %VM
        ha1(1,sta) = subplot(4,noSTA,0*noSTA+sta);
        t0 = 0;
        for k = 1:numel(res.VM)
            t = res.VM(k).t+t0;
            d = res.VM(k).data;
            plot(t,d); hold on;
            t0 = t(end);
        end
        title(sprintf('VM - %s',stage))
        
        %EMG
        ha1(2,sta) = subplot(4,noSTA,1*noSTA+sta);
        t0 = 0;
        for k = 1:numel(res.EMG)
            t = res.EMG(k).t+t0;
            d = detrend(res.EMG(k).data);
            plot(t,d); hold on;
            t0 = t(end);
        end
        title(sprintf('Detrended EMG - %s',stage))
        
        
        %PUPIL
        ha1(3,sta) = subplot(4,noSTA,2*noSTA+sta);
        t0 = 0;
        for k = 1:numel(res.pupil)
            t = res.pupil(k).t+t0;
            d = res.pupil(k).data;
            plot(t,d); hold on;
            t0 = t(end);
        end
        title(sprintf('pupil - %s',stage))
        
        
        %MOTION
        ha1(4,sta) = subplot(4,noSTA,3*noSTA+sta);
        t0 = 0;
        for k = 1:numel(res.motion)
            t = res.motion(k).t+t0;
            d = res.motion(k).data;
            plot(t,d); hold on;
            t0 = t(end);
        end
        title(sprintf('motion - %s',stage))
        
        linkaxes(ha1(:,sta),'x')
    end
    
    %     return
    %     %%
    
%     SAVE VM
    [~,rFile] = fileparts(files.vmStage);
    sFile = sprintf('%s_basicAna_10sec.mat',rFile);
    save(fullfile(sPath,sFile),'-struct','RES')
    fprintf('%s Saved: %s\n',indent,sFile)
end %loop files