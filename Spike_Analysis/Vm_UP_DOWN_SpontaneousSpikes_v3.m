clc; clear; close all
addpath(genpath('X:\OptoLab_v4.1\function'))

%%
%SAVE OUTPUT
savePath = 'X:\1 TIDIS Lab\Tiago\Campelo_Manuscript\Fig4_PFC_Diazepam\RawData\';
saveFile = 'DZP_SpontSpikes.mat';
condition ='DZP_';
fullSavePath = fullfile(savePath, saveFile);
savefile = true;

splitWake = false;

%LOAD FILES
%BasicDataAnalysis struct with information on Quiet/Motion ID
Condition = 'DZP_BasicAnalysis_10sec';
tmp = load (['X:\1 TIDIS Lab\Tiago\Campelo_Manuscript\Fig4_PFC_Diazepam\RawData\', Condition]);
Basic_DATA = tmp.(sprintf('%sDATA', condition)); 

%Vmstage file containing EEG and corresponding Vm traces
tmp = dir(['X:\1 TIDIS Lab\Tiago\Ephys_PROCESSED\In_vivo\Awake_Sleep\PFC_Processed\Pharmacology\Diazepam\*\*\*_vmStage.mat']);
tmp([tmp.isdir]) = [];
if isempty(tmp)
    fprintf(2,'No files found for: %s\n',files.vmStage);
    return
end
Files = selectionList(fullfile({tmp.folder},{tmp.name}));
if isempty(Files)
    fprintf(2,'No file selected\n');
    return
end

%Define main stages to be analyzed
stages = {'Wake','NREM','REM'};
clear tmp

%% Properties for spike detection
%SPIKE DETECTION - Based on spikeRemoval from T. Rusterholz (23 Jun 2021)
par.threshold = -20; %[mV]
par.fun = @(x,fs)x-movmedian(x,ceil(20/1000*fs/2)*2+1); %odd nbins
par.margin = [3,5]; %[ms]
%interpolation method (using interp1)
par.interpMethod = 'pchip';
%ARTIFACTS
% minimal difference between peaks (artifacts if smaller)
art.minDiff = 0.001; %[s], set to zero if not wanted
%% Load Vm data; filter Vm data;
DATA = [];

%Initiate file loop (across all the recordings)
for fil = 1:numel(Files)
    %Initiate structs to load Vm data, BasicAnalysis data and cell ID
    File_Vm = Files{fil};
    [rPath,rFile] = fileparts(File_Vm);
    tmp_Vm = load(File_Vm);
    label = strrep (rFile, '_vmStage', '');
    fieldName = genvarname(label);
    fsDAT = tmp_Vm.fs;
%     numericPart = regexp(rFile, '[\do]+', 'match');
    numericPart = strrep(rFile, '_vmStage', '');
    IDcell = ['x', numericPart];
    
    %Load traces ID for Quiet and Motion Wake
    if isfield (Basic_DATA.(IDcell), 'Wake')
        if splitWake
            Quiet_ID = Basic_DATA.(IDcell).Wake.Quiet_ID;
            Motion_ID = Basic_DATA.(IDcell).Wake.Motion_ID;
        else
            A = Basic_DATA.(IDcell).Wake.Quiet_ID;
            B = Basic_DATA.(IDcell).Wake.Motion_ID;
            Quiet_ID = [A; B];
            Motion_ID = Quiet_ID;
        end
    end
    
    %Initiate stages loop
    for sta = 1:numel (stages)
        stage = stages{sta};
        %Load the traces ID from Basic Analysis and corresponding to each
        %stage
        if sta == 1 && isfield (Basic_DATA.(IDcell), 'Wake')
            if splitWake
                Trace_ID = [Quiet_ID; Motion_ID]; %Merge both for Wake: splot later
            else
                Trace_ID = Quiet_ID;
            end
            %Initiate structs for Quiet_Wake
            UP_states_bout.Quiet_Wake = [];
            DOWN_states_bout.Quiet_Wake = [];
            %Initiate structs for Motion_Wake
            UP_states_bout.Motion_Wake = [];
            DOWN_states_bout.Motion_Wake = [];
            
        elseif sta >= 2 && isfield (Basic_DATA.(IDcell), stage)
            Trace_ID = Basic_DATA.(IDcell).(stage).VM.trace; 
            UP_states_bout.(stage) = [];
            DOWN_states_bout.(stage) = [];
            
        end %condition stage
        
        %Loop across all bouts
        for bou = 1: numel (tmp_Vm.(stage))
            %But only analyse if the bout was considered in basic analysis
            %(traces_ID filtering)
            
            if any (ismember (Trace_ID, bou))
                %if is member, load the trace/bout for further analysis
                bout = ['b', num2str(bou)];
                peaks = [];
                spont_spikes = [];
                dat = tmp_Vm.(stage)(:);
                dat = dat{bou,1};
                Vm = dat(~isnan(dat));
                Vm_ts = (1:numel(Vm))/fsDAT;
                Vm_ts = Vm_ts(:);
                
                %Vm to use as min (Vm) for spike threshold
                Spk_thresholding = min (Vm);
                
                %Vm Filtering
                fs = tmp_Vm.fs;
                Vm_filt = movmean(Vm - min(Vm), fs/10);
                Vm_filt = Vm_filt(:); %force to column vector
                t = (1:numel(Vm_filt))/fs;
                N = numel(Vm_filt);
                
                %FIND PEAKS
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
                
                if false %False: no plot; Also: nanmean removes spikes!
                    %PLOT PEAKS
                    %close all; clc
                    figure('WindowState','maximized')
                    t = (1:numel(Vm_filt))/fs;
                    plot(t, Vm_filt,'linewidth',2); hold on
                    %up/down
                    labels = Labels(:,1);
                    colors = {'b','r'};
                    for k = 1:size(Labels)
                        label = labels{k};
                        color = colors{k};
                        ind   = strcmpi({peaks.label},label);
                        locs  = [peaks(ind).loc];
                        plot(t(locs), Vm_filt(locs),':','color',color,...
                            'marker','o','linewidth',2)
                        %threshold
                        thres = peaks(find(ind,1,'first')).thres;
                        h = yline(thres,':',sprintf('Threshold %s',label));
                        h.Color = color;
                    end
                    %duration
                    for k = 1:numel(peaks)
                        peak  = peaks(k);
                        color = colors{strcmpi(labels,peak.label)};
                        ind   = [peak.indSTA, peak.indEND];
                        y     = [1,1]*peak.thres;
                        plot(t(ind), y, 'color',color, 'linewidth',2)
                        text(mean(t(ind)),y(1),...
                            sprintf('%.2g s',peak.duration),...
                            'horizontalalignment','center',...
                            'verticalalignment','bottom')
                    end
                end %Plot loop
                
                %% Calculate UP states spikes
                cnt = 0;
                for UPz = 1:size ({peaks.peak},2)
                    if strcmpi (peaks(UPz).label, 'UP')
                        cnt = cnt+1;
                        UP_id = ['UP', num2str(cnt)];
                        if peaks(UPz).indSTA > 0.1*fsDAT && peaks(UPz).indEND + 0.1*fsDAT < numel (Vm)
                            A = peaks(UPz).indSTA - 0.1*fsDAT;
                            B = peaks(UPz).indEND + 0.1*fsDAT;
                        else
                            A = peaks(UPz).indSTA;
                            B = peaks(UPz).indEND;
                        end
                        
                        Vm_spk = Vm (A:B);
                        Vm_spk_ts = Vm_ts (A:B);
                        data_UPstate = Vm_spk;
                        
                        if ~any(Vm_spk >= par.threshold)
                            spont_spikes(UPz).noSPK = 0;
                            spont_spikes(UPz).rate = 0;
                            spont_spikes(UPz).amp = 0;
                            spont_spikes(UPz).threshold = 0;
                            spont_spikes(UPz).isi = 0;
                            
                        else
                            indSTA_spk = find(...
                                [Vm_spk;-inf]>=par.threshold & ...
                                [-inf;Vm_spk]<par.threshold);
                            indEND_spk = find(...
                                [Vm_spk;-inf]<par.threshold  & ...
                                [-inf;Vm_spk]>=par.threshold)-1;
                            
                            noSPK = numel(indSTA_spk);
                            indPKS = NaN(noSPK,1);
                            t_PKS = NaN(noSPK,1);
                            
                            % Spike Analysis
                            for spk = 1:noSPK
                                %spike peak
                                ind = indSTA_spk(spk):indEND_spk(spk);
                                [~,indP] = max(Vm_spk(ind));
                                indPKS(spk) = indP+ind(1)-1;
                                t_PKS(spk) = Vm_spk_ts (indPKS(spk));
                                %amplitude
                                amp(spk) = Vm_spk(indPKS(spk))- min (Vm_spk);
                            end
                            
                            if indPKS (1) > 3
                                threshold = Vm_spk (indPKS(1)-2) - Spk_thresholding;
                            else
                                threshold = Vm_spk (indPKS(1)) - Spk_thresholding;
                            end
                            
                            %interspike interval (msec.)
                            isi  = (diff(indPKS)/fs)* 1000;
                            %Spiking rate
                            rate = noSPK/(numel (Vm_spk)/fs);
                        
                            if rate > 450
                                noSPK = NaN;
                                rate = NaN;
                                amp = NaN;
                                threshold = NaN;
                                isi = NaN;
                            else
                                %AppendSpike data into a new spike struct
                                DATA.(IDcell).(stage).data_UPstate.(bout).(UP_id) = data_UPstate; %Append ALL UPstates withSpikes
                            end %clear data to cut spikes with ISI > duration of UP state
                        
                            %filter data to remove weirdo values
                            if noSPK > 50
                                noSPK = NaN;
                            elseif rate > 20
                                rate = NaN;
                            end
                            
                            spont_spikes(UPz).noSPK = noSPK;
                            spont_spikes(UPz).rate = rate;
                            spont_spikes(UPz).amp = nanmean(amp);
                            spont_spikes(UPz).threshold = nanmean (threshold);
                            spont_spikes(UPz).isi = nanmean(isi);
                        end %condition: UP state w/ Spike or No-Spike
                    end % Condition: Just UP states
                end %Loop UP states
                                
                %UP states / DOWN states ID selection
                bout_peaks.label = {peaks.label}';
                label_DOWN = [];
                label_UP = [];
                
                for lab = 1:size (bout_peaks.label,1)
                    label_ID = bout_peaks.label{lab,1};
                    if strcmpi (label_ID, 'Down')
                        label_DOWN = [label_DOWN; lab];
                        
                    elseif strcmpi (label_ID, 'UP')
                        label_UP = [label_UP; lab]; 
                    end %condition label_UP/DOWN
                end %lab loop
                
                %Divide DATA into UP and DOWN states
                DATA_subfields = {'peak', 'duration','noSPK', 'rate','amp','threshold','isi'};
                
                for subf = 1: size (DATA_subfields,2)
                    subfieldz = DATA_subfields{1, subf};
                    
                    if subf < 3
                        data = [peaks.(subfieldz)]';
                        bouts.peaks.Down.(subfieldz) = data(label_DOWN);
                        bouts.peaks.Up.(subfieldz) = data(label_UP);
                        
                    elseif subf > 2
                        results = [];
                        % Loop to remove empty results
                        for i = 1:numel(label_UP)
                            idx = label_UP(i);
                            data = spont_spikes(idx).(subfieldz);
                            data_UP = data(:);
                            results = [results; data_UP];
                        end %label UP
                        bouts.spikes.Up.(subfieldz) = results;
                    end %condition subfield
                end %subfield loop
                
                peak_state = unique(bout_peaks.label');
                
                for peakl = 1: size (peak_state,2)
                    peak_type = peak_state{1, peakl};
                    
                    if sta == 1
                        if splitWake
                            if any(ismember (Quiet_ID, bou))
                                DATA.(IDcell).Quiet_Wake.(bout).(peak_type).Peaks = bouts.peaks.(peak_type);
                                DATA.(IDcell).Quiet_Wake.(bout).Up.Spikes = bouts.spikes.Up;
                            elseif any (ismember (Motion_ID, bou))
                                DATA.(IDcell).Motion_Wake.(bout).(peak_type).Peaks = bouts.peaks.(peak_type);
                                DATA.(IDcell).Motion_Wake.(bout).Up.Spikes = bouts.spikes.Up;
                            end %condition Trace_ID-bou
                        else
                            DATA.(IDcell).Wake.(bout).(peak_type).Peaks = bouts.peaks.(peak_type);
                            DATA.(IDcell).Wake.(bout).Up.Spikes = bouts.spikes.Up;
                        end
                    
                    elseif sta >= 2
                        if any (ismember (Trace_ID, bou))
                            DATA.(IDcell).(stage).(bout).(peak_type).Peaks = bouts.peaks.(peak_type);
                            DATA.(IDcell).(stage).(bout).Up.Spikes = bouts.spikes.Up;
                        end %condition Trace_ID-bout
                    end %condition stages
                end %condition peak_state
            end %condition UP or DOWN state
        end %condition Trace_ID - Bout
    end %bouts loop
end %stages loop
%% STAT File
STAT = [];

cell_list = fieldnames (DATA);
if splitWake
    stages2 = {'Quiet_Wake', 'Motion_Wake', 'NREM', 'REM'};
else
    stages2 = {'Wake', 'NREM', 'REM'};
end %splitWake loop

%STAT loop for DATA
for sta = 1: size (stages2, 2)
    stage = stages2{1,sta};
    
    for p = 1:size (peak_state,2)
        state_ID = peak_state{1,p};
        
        if p == 1
            mainfields = {'Peaks'};
        else
            mainfields = {'Spikes', 'Peaks'};
        end %condition subfields
        
        for cond = 1: size(mainfields,2)
            conditionz = mainfields{1,cond};
            
            if strcmpi (conditionz, 'Peaks')
                subfields = {'peak', 'duration'};
            elseif strcmpi (conditionz, 'Spikes')
                subfields = {'noSPK', 'rate','amp','threshold','isi'};
            end % condition subfie
            
            for sub = 1: size (subfields,2)
                subconditionz = subfields{1,sub};
                all_bouts_data = [];
                average_cell = [];
                
                JS_all_bouts = [];
                JS_cell_average = [];
                
                for fil = 1: size (cell_list, 1)
                    IDcell = cell_list {fil, 1};
                    
                    if isfield (DATA.(IDcell),stage)
                        bout_list = fieldnames(DATA.(IDcell).(stage));
                        bout_list(strcmpi(bout_list,'data_UPstate')) = [];
                        cell_data = [];
                        
                        
                        if strcmpi(stage,'NREM') && strcmpi(state_ID,'up') ...
                                && strcmpi(subconditionz,'noSPK')
                            aaa=0;
                        end
                        
                        for bou = 1: size (bout_list, 1)
                            bout = bout_list {bou,1};
                            if isfield (DATA.(IDcell).(stage).(bout), (state_ID))
                                data = DATA.(IDcell).(stage).(bout).(state_ID).(conditionz)...
                                    .(subconditionz);
                                cell_data = [cell_data; data];
                            end %isfield loop
                        end %bout_loop
                        
                        %Concatenate all bouts data from all the cells
                        all_bouts_data = [all_bouts_data; cell_data];
                        average_cell = [average_cell; nanmean(cell_data)];

                        %STATs only on UP states with Spikes
                        UPflags = find (cell_data > 0);
                        if numel (UPflags) > 0 && strcmp(conditionz, "Spikes")
                            cell_dataJS = cell_data (UPflags);
                            JS_all_bouts = [JS_all_bouts; cell_dataJS];
                            JS_cell_average = [JS_cell_average; (nanmean (cell_dataJS))];
                        end %condition UPflags
                    end %Isfield IDcell - stage condition
                end %fil loop
                STAT.(stage).(state_ID).(conditionz).All_Bouts.(subconditionz) = all_bouts_data;
                STAT.(stage).(state_ID).(conditionz).Cell_averages.(subconditionz) = average_cell;
                if strcmp(conditionz, "Spikes")
                    STAT.(stage).(state_ID).Spikes.All_Bouts_OnlySpikes.(subconditionz) = JS_all_bouts;
                    STAT.(stage).(state_ID).Spikes.Cell_averages_OnlySpikes.(subconditionz) = JS_cell_average;
                end
            end %subfields loop
        end %condition loop
    end %peak state loop
end %stages loop
%% Plot histograms - UP state duration
figure (1)
if splitWake
    stages_plot = {'Quiet_Wake','Motion_Wake', 'NREM','REM'};
else
    stages_plot = {'Wake', 'NREM','REM'};
end %stages loop

sgtitle('Distribution of duration of UP states (sec)');
for sta = 1:numel(stages_plot)
    stageplt = stages_plot{sta};
    
    subplot (2,2,sta)
    histogram (STAT.(stageplt).Up.Peaks.All_Bouts.duration);
    title (sprintf(stageplt));
end

%% Plot histograms - Spikes Distriubtion
figure (2)

sgtitle('Distribution of spikes per UP states');
for sta = 1:numel(stages_plot)
    stageplt = stages_plot{sta};
    
    subplot (2,2,sta)
    histogram (STAT.(stageplt).Up.Spikes.All_Bouts.noSPK);
    title (sprintf(stageplt));
end

%% SAVE THE DATA & STAT STRUCT
if savefile
    newDATA = [condition 'DATA'];
    assignin('base', newDATA, evalin('base', 'DATA'));
    newSTAT = [condition 'STAT'];
    assignin('base', newSTAT, evalin('base', 'STAT'));
    save (fullSavePath, newDATA,newSTAT);
end %save only if true