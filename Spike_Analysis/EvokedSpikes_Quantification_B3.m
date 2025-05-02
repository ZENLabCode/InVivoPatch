% Spontaneous Spikes Quantification - Tiago Campelo (01/09/2024)
clc; clear; close all
addpath(genpath('X:\OptoLab_v4.1\function'));

%% PARAMETERS
%----------
%FILES
files = 'X:\1 TIDIS Lab\Tiago\Ephys_PROCESSED\In_vivo\Awake_Sleep\PFC_Processed\Pharmacology\Diazepam\*\Evoked_Spikes\*\*_vmStage.mat';

%SAVE OUTPUT
savePath = 'X:\1 TIDIS Lab\Tiago\Campelo_Manuscript\Fig4_PFC_Diazepam\RawData';
saveFile = 'DZP_EvokedSpikes.mat';
savecondition ='DZP_';
fullSavePath = fullfile(savePath, saveFile);
savefile = true;

%SPIKE DETECTION - Based on spikeRemoval from T. Rusterholz (23 Jun 2021)
par.threshold = 0; %[mV]
par.fun = @(x,fs)x-movmedian(x,ceil(20/1000*fs/2)*2+1); %odd nbins
par.margin = [3,5]; %[ms]
%interpolation method (using interp1)
par.interpMethod = 'pchip';
%ARTIFACTS
% minimal difference between peaks (artifacts if smaller)
art.minDiff = 0.001; %[s], set to zero if not wanted

%STAGES (case sensitive! variable names in read files)
stages = {'Wake','NREM','REM'};

%Timestamps of stimulation: might need to be tuned for each cell
Evoked_timestamps = [21964, 50964, 79964];

%% MAIN SCRIPT
scriptName = mfilename;
fprintf('%s\n%s\n',scriptName,repmat('-',size(scriptName)))

%FIND/SELECT FILES
tmp = dir(files);
if isempty(tmp)
    fprintf(2,'No files found for: %s\n',files)
    return
end
Files = selectionList(fullfile({tmp.folder},{tmp.name}));
noFIL = numel(Files);
if noFIL==0
    fprintf(2,'No file selected\n')
    return
end

%FILE LOOP
nnFIL = numel(num2str(noFIL));
indent = blanks(2*nnFIL+2);
DATA = {};

for fil = 1:noFIL
    fname = Files{fil};
    [rPath,rFile] = fileparts(fname);
    fprintf('%*i/%i: %s\n',nnFIL,fil,noFIL,fname);
    %create TTL filename
    rPathT = fileparts(rPath);
    rFileT = strrep(rFile,'_vmStage','_IN1.mat');
    fnameT = fullfile(rPathT, rFileT);
    if exist(fnameT,'file')~=2
        fprintf(2,'%s TTL file not available\n',indent)
        continue
    end
    data_clusters = NaN(noFIL, 1); %Check the need of this!
    
    %Load Hypnogram
    rHypno = '\Hypnogram.mat';
    Hypnogram_path = load ([rPath,rHypno]);
    cell_ID = strrep(rFile, '_vmStage', '');
    cell_ID = ['f', cell_ID];
    
    %Load Evoked Spikes - Vm
    lastBackslashIndex = find(rPath == '\', 1, 'last');
    rPath2 = rPath(1:lastBackslashIndex);
    pattern = '(\d+)_vmStage';
    replacement = '$1_IN0';
    rFile2 = regexprep(rFile, pattern, replacement);
    fname2 = [rPath2, rFile2];
    Vm_fulldata = load (fname2);
    Vm_fs = Vm_fulldata.SampRate;
    Vm = Vm_fulldata.resampled_data_mV;
    Vm_ts = (1:numel(Vm))/Vm_fs;
    %Load TTL for artifact removal
    TTL = load(fnameT);
    tmp = [min(TTL.resampled_data_mV), max(TTL.resampled_data_mV)];
    thr = tmp(1)+0.1*diff(tmp);
    t1  = (find(TTL.resampled_data_mV>thr,1,'first')-1)/TTL.SampRate -0.05;
    t2  = (find(TTL.resampled_data_mV>thr,1,'last' )+1)/TTL.SampRate +0.05;
    art = false(size(Vm_ts));
    art(Vm_ts>t1 & Vm_ts<t2) = true;
    if false %test plot
        Vm_tmp = Vm;
        Vm(art) = NaN;
        figure;
        subplot(211);
        plot(Vm_ts, Vm_tmp); hold on
        plot(Vm_ts, Vm, 'linewidth',2)
        subplot(212);
        ttl = TTL.resampled_data_mV;
        plot(ttl); hold on
        ttl(art) = 0;
        plot(ttl, 'linewidth',2); zoom xon
        return
    else
        Vm(art) = NaN;
    end
    
    %Resample Hypnogram to match Vm
    Hypnogram = resample (Hypnogram_path.Hypnogram, Vm_fs,1); %Hypnogram fs = 1;
    Hypnogram = round (Hypnogram);
    Hypnogram_ts = (1:numel(Hypnogram))/Vm_fs;
    
    %Crop the Vm trace to the last numerical value of the Hypnogram
    nanID = find (isnan(Hypnogram));
    if ~isempty (nanID)
        First_NanID = nanID(1);
        Vm = Vm(1:First_NanID);
        Vm_ts = Vm_ts(1:First_NanID);
        Hypnogram = Hypnogram (1:First_NanID-1);
        Hypnogram_ts = Hypnogram_ts (1:First_NanID-1);  
    end
    
    %Define the number of bouts on the evoked recording protocol
    Bouts = floor((numel (Vm) / (10*Vm_fs)));
    Total_ts = [];
    
    %Extract timestamps across all the bouts
    for rep = 2:Bouts
        % Calculate the offset for the current repetition
        offset = (rep - 1) * 10*Vm_fs; 
        % Add the offset to the base timestamps and append to the list
        currentTimestamps = Evoked_timestamps + offset;
        Total_ts = [Total_ts; currentTimestamps];
    end
    
    %All TimeStamps vector initiation
    allTimestamps = [Evoked_timestamps; Total_ts];
    reshaped_allTimestamps = reshape (allTimestamps, [], 1);
    reshaped_allTimestamps = sort (reshaped_allTimestamps, 1);
    
    %Find indexes and number of spikes on the full Vm trace
    indSTA = find(...
        [Vm;-inf]>=par.threshold & ...
        [-inf;Vm]<par.threshold);
    indEND = find(...
        [Vm;-inf]<par.threshold  & ...
        [-inf;Vm]>=par.threshold)-1;
    
    %Extract indexes and timming of all spikes
    noSPK = numel(indSTA);
    indPKS = NaN(noSPK,1);
    t_PKS = NaN(noSPK,1);
    
    if noSPK > 0
        for spk = 1:noSPK %start/end correction
            %spike peak
            ind = indSTA(spk):indEND(spk);
            [~,indP] = max(Vm(ind));
            indPKS(spk) = indP+ind(1)-1;
            t_PKS(spk) = Vm_ts(indPKS(spk));
        end %spk loop
    end
     
    %Cluster spikes according to their occurrence in time
    clear Spikes
    cnt = 0;
    for k = 1:size(reshaped_allTimestamps,1)
        stim_idx = reshaped_allTimestamps(k);
        % Determine spikes wihin the evoked timestamps
        SpikesFlag1 = find (indPKS > stim_idx);
        SpikesFlag2 = find (indPKS < (stim_idx + 5000));
        Spikesoverlap = intersect(SpikesFlag2, SpikesFlag1);
        N = numel(Spikesoverlap);
        %append
        Spikes(k).cluster = k;
        
        %Spikes(k).cluster ones (numel(Spikesoverlap), 1)*k; 
        Spikes(k).N  = N;
        if N==0
            Spikes(k).indPKS = stim_idx + 500; %500: filter
        else
            Spikes(k).indPKS = indPKS (Spikesoverlap);
        end
        Spikes(k).t_PKS  = t_PKS  (Spikesoverlap);
    end
    %test plot
    if true
        figure
        plot(Vm,'k'); hold on;
        for k = 1:numel(Spikes)
            if Spikes(k).N>0
                ind = Spikes(k).indPKS;
                plot(ind,Vm(ind),'o')
            end
        end
        title(strrep(cell_ID,'_','\_'))
        zoom xon
    end
    
    
    %Extract transitions based on the Hypnogram
    Hypno_change = NaN(0,3); %[stage, start index, end index]
    ind1 = 1; %init
    while true
        stage = Hypnogram(ind1);
        ind2  = find(Hypnogram(ind1:end)~=stage,1,'first')+ind1-1;
        if isempty(ind2)
            Hypno_change(end+1,:) = [stage,ind1,numel(Hypnogram)];
            break
        end
        Hypno_change(end+1,:) = [stage,ind1,ind2];
        ind1 = ind2+1;
    end
    
    %%
    Vm2filt = Vm;
    for k = 1:size (reshaped_allTimestamps,1)
        Stim_STA = reshaped_allTimestamps(k) - 800;
        Stim_END = reshaped_allTimestamps(k) + 4000;
        
        Vm2filt(Stim_STA:Stim_END) = NaN;
    end %loop across stimulations
    
    %Vm Filtering
    fs = Vm_fulldata.SampRate;
    Vm_filt = movmean(Vm2filt - min(Vm2filt), fs/10);
    Vm_filt = Vm_filt(:); %force to column vector
    t = (1:numel(Vm_filt))/fs;
    N = numel(Vm_filt);
    
    if false %plot filtering of Vm
        close all; 
        figure (fil); 
        subplot (3,1,1); plot (Hypnogram); hold on;
        subplot (3,1,2); plot (Vm_filt);
        subplot (3,1,3); plot (Vm)
    end 
    
    Labels = {... {label, threshold, factor, fun valid peaks, style}
        'Down',prctile(Vm_filt,15),-1,@(x,thr)x<thr;...
        'Up'  ,(prctile(Vm_filt,15) + (prctile(Vm_filt,15)*0.3)),+1,@(x,thr)x>thr;...  %threshold: prctile(Vm_filt,40) or 50
        };
    cnt = 0;
    
    for k = 1:2
        [label, thres, fac, fun] = Labels{k,:};
        [pks,locs] = findpeaks(fac*Vm_filt);
        pks = pks*fac; %correct values again
        
        %LIMITS
        %threshold
        indOK = fun(pks,thres);
        locs(~indOK) = [];
        pks(~indOK)  = [];
        %combine peaks when not reaching threshold in between
        indOK = false(size(locs));
        ind1  = 1;
        while true
            ind2 = ind1+1;
            while ind2+1<numel(locs)
                ind = locs(ind1):locs(ind2);
                if all(fun(Vm_filt(ind), thres))
                    ind2 = ind2+1;
                else
                    break
                end
            end
            ind0 = ind1:ind2-1;
            [~,ind] = max(fac*pks(ind0)>fac*thres);
            indOK(ind0(ind(1))) = true;
            ind1 = ind2;
            if ind2+1>numel(locs)
                break
            end
        end
        locs(~indOK) = [];
        pks(~indOK)  = [];
        %check peaks
        locs = [locs;inf]; %add inf for time diff limit
        pks  = [pks ;NaN]; %add NaN (same size as locs)
        while true
            %combine peaks based on time limit
            ind1 = find(diff(locs)>0.25*fs, 1, 'first');
            [~,ind] = max(fac*pks(ind1));
            loc = locs(ind1(ind));
            pk  = pks(ind1(ind));
            
            %find start/end point
            tmp = fun(Vm_filt(1:locs(ind1(1))), thres);
            indSTA = find(~tmp, 1,'last');
            tmp = fun(Vm_filt(locs(ind1(end)):end), thres);
            indEND = find(~tmp, 1,'first')+locs(ind1(end))-1;
            
            %append data
            if ~isempty(indSTA) && ~isempty(indEND)
                cnt = cnt+1;
                peaks(cnt).label    = label;
                peaks(cnt).thres    = thres;
                peaks(cnt).loc      = loc;
                peaks(cnt).peak     = pk;
                peaks(cnt).indSTA   = indSTA;
                peaks(cnt).indEND   = indEND;
                peaks(cnt).duration = (indEND-indSTA+1)/fs;
            end
            
            %cut
            locs(1:ind1(end)) = [];
            pks(1:ind1(end)) = [];
            if numel(pks)==1 %the added NaN
                break
            end
        end
    end

    %Select UP/DOWN START & END indexes
    UP_idx = [];
    DOWN_idx = [];
    
    for Pkz = 1:size ({peaks.peak},2)
        if strcmpi (peaks(Pkz).label, 'Up')
            UP = peaks(Pkz).indSTA;
            UP_idx = [UP_idx; UP];
        elseif strcmpi (peaks(Pkz).label, 'Down')
            DOWN = peaks(Pkz).indSTA;
            DOWN_idx = [DOWN_idx; DOWN];
        end %UP/DOWN segmentation
    end %Pkz loop
    
    if false
        close all;
        UP_Y_idx = ones (numel(UP_idx),1)*-60;
        DOWN_Y_idx = ones (numel(DOWN_idx),1)*-70;
        figure (fil);
        subplot (2,1,1)
        plot (Hypnogram);
        subplot (2,1,2)
        plot (Vm); hold on;
        scatter (UP_idx, UP_Y_idx);
        scatter (DOWN_idx, DOWN_Y_idx); zoom xon
    end

%%  Filter all the spike data according to parameters to be analyzed  
    %1) Timestamps of evoked currents
    % Initialize a template structure
    Stim1 = allTimestamps(:, 1);
    Stim2 = allTimestamps(:, 2);
    Stim3 = allTimestamps(:, 3);
    for i = 1:numel(Spikes) 
        spikes_cluster = Spikes(i).cluster;
        spikes_N       = Spikes(i).N;
        spikes_ind     = Spikes(i).indPKS';
        cluster_id     = sprintf('cluster_%i',spikes_cluster);
        
        condition1 = mean(spikes_ind) > Stim1 & mean(spikes_ind) < Stim1 + 4000;
        condition2 = mean(spikes_ind) > Stim2 & mean(spikes_ind) < Stim2 + 4000;
        condition3 = mean(spikes_ind) > Stim3 & mean(spikes_ind) < Stim3 + 4000;
        
        for sta = 1:numel (stages)
            stage = stages{sta};
            hypno_change = Hypno_change(Hypno_change(:,1)==sta,:);
                        
            for a = 1:size(hypno_change,1) %loop across timestamps where evoked occured at stage
                %if belongs to stage
                if mean (spikes_ind) > hypno_change(a,2) && mean (spikes_ind) < hypno_change(a,3)
                    %type UP or DOWN
                    previousUP   = max(UP_idx(UP_idx<spikes_ind (1)));
                    previousDOWN = max(DOWN_idx (DOWN_idx <spikes_ind (1)));
                    if isempty(previousUP) && isempty(previousDOWN)
                        Cluster_type = 'UNKNOWN';
                        Cluster_ts   = NaN;
                    elseif isempty(previousDOWN) || previousUP>previousDOWN
                        Cluster_type = 'UP';
                        Cluster_ts   =previousUP/Vm_fs;
                    else
                        Cluster_type = 'DOWN';
                        Cluster_ts   =previousDOWN/Vm_fs;
                    end

                    if any (condition1)
                        DATA.(cell_ID).(stage).Evoked1.spike.(cluster_id) = spikes_N;
                        DATA.(cell_ID).(stage).Evoked1.cluster_type.(cluster_id) = Cluster_type;
                        DATA.(cell_ID).(stage).Evoked1.cluster_ts.(cluster_id) = Cluster_ts;
                    elseif any (condition2)
                        DATA.(cell_ID).(stage).Evoked2.spike.(cluster_id) = spikes_N;
                        DATA.(cell_ID).(stage).Evoked2.cluster_type.(cluster_id) = Cluster_type;
                        DATA.(cell_ID).(stage).Evoked2.cluster_ts.(cluster_id) = Cluster_ts;
                    elseif any (condition3)
                        DATA.(cell_ID).(stage).Evoked3.spike.(cluster_id) = spikes_N;
                        DATA.(cell_ID).(stage).Evoked3.cluster_type.(cluster_id) = Cluster_type;
                        DATA.(cell_ID).(stage).Evoked3.cluster_ts.(cluster_id) = Cluster_ts;
                    end %sort spikes across evoked timestamps
                end %spikes-stimulation loop
            end % row loop
        end %stage loop
    end %lcluster loop
end %FIL loop
%% Create STAT file
STAT = [];
MainFields = {'Evoked1','Evoked2', 'Evoked3'};
SubFields = {'spike','cluster_type', 'cluster_ts'};
fieldCells = fieldnames(DATA);

for main = 1:numel(MainFields)
    mainfield = MainFields{main};
    
    for sub = 1:numel (SubFields)
        subfield = SubFields{sub};
        
        for sta = 1:numel(stages)
            stage = stages{sta};
            average_spikes = [];
            all_clusters = [];
        
            for field = 1:numel (fieldCells)
                cell_ID = fieldCells{field};            
                data_clusters = [];
                
                if isfield (DATA.(cell_ID), stage)
                    if isfield (DATA.(cell_ID).(stage), mainfield)
                        clusterz_ID = fieldnames (DATA.(cell_ID).(stage).(mainfield).(subfield));
                        
                        for clust = 1:numel(clusterz_ID)
                            clusterz = clusterz_ID{clust};
                            data = DATA.(cell_ID).(stage).(mainfield).(subfield).(clusterz);
                            
                            if sub == 2
                                 data_clusters{end+1} = data;
                            else
                                data_clusters = [data_clusters, data];
                            end %loop data concatenation                           
                        end %loop clusters
                        
                    end %isfield DATA-MainField
                end %isfield DATA-Subfield
                
                all_clusters = [all_clusters, data_clusters];
                
                if sub == 1
                    average_SPR = nanmean (data_clusters);
                    average_spikes = [average_spikes; average_SPR];
                end %average per cell loop
                
            end %loop field (cell_ID)
            
            if sub == 1
                STAT.(mainfield).(stage).Average_nrSpikes = average_spikes';
            end %conditional saving
            STAT.(mainfield).(stage).(subfield) = all_clusters;
        end %loop stages
    end %loop subfields
end %loop mainfields

%% Save the output files
if savefile
    newDATA = [savecondition 'DATA'];
    newSTAT = [savecondition 'STAT'];
    assignin('base', newDATA, evalin('base', 'DATA'));
    assignin('base', newSTAT, evalin('base', 'STAT'));
%     save (fullSavePath, newDATA);
    save (fullSavePath, newDATA,newSTAT);
end
