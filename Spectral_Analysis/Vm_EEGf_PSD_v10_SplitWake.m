% Load Vm and EEG files into a filelist
close all;  clear; clc; addpath(genpath('X:\OptoLab_v4.1\function'));

%SAVE OUTPUT
savePath = 'X:\1 TIDIS Lab\Tiago\Campelo_Manuscript\Fig1_BasicAnalysis_PFC\Raw_Data';
saveFile = 'PFC_Analysis_PSD_VmEEG_PLOTTING.mat';
savecondition ='PFC_';
fullSavePath = fullfile(savePath, saveFile);
savefile = false;

%LOAD FILES
%BasicDataAnalysis struct with information on the traces with Quiet/Motion ID
Condition = 'PFC_BasicAnalysis_5sec';
tmp = load (['X:\1 TIDIS Lab\Tiago\Campelo_Manuscript\Fig1_BasicAnalysis_PFC\Raw_Data\', Condition]);
Basic_STAT = tmp.(sprintf('%sSTAT', savecondition));
Basic_DATA = tmp.(sprintf('%sDATA', savecondition)); 
folder_ID ='5sec_window'; %String to remove from load file directory to load EEG and Hypnogram
window = 5; %BasicAnalysis: 5 seconds windows of Vm

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

%% Define stages to be analyed; Filtering properties; Frequency bands to be analyzed in FFT.

% Split the Wake into Quiet/Motion_Wake
splitWake = false;

stages = {... {Stages, Filter settings, RGB combination}
    'Quiet_Wake', [13  114 116]'/255; ...
    'Motion_Wake', [57  73 112]'/255; ...
    'NREM', [234 174  51]'/255;...
    'REM', [218  93  41]'/255; ...
    };
%Frequencies bands where EEG & Vm are quantified
Bands = {... {frequency band, field label}
    [0.5 2] , 'slowWave';...
    [2 4]  , 'delta';...
    [6 9]  , 'theta';...
    [8 12] , 'alpha';...
    [10 16] , 'spindle';...
    [30 60], 'gamma30_60';...
    };

%Compute the number of bands
noBND = size(Bands,1);

%FFT EEGcut & VM settings
frqFFT = 0:0.1:60;
wd_fun = @(n)hamming(n)*2;

%Initiate DATA vector to store all the data per recorded cell
DATA = [];
%% Load Files, Filter EEG/Vm traces and subdivide EEG/Vm traces between Quiet and (Motion)Wake
%Define stages to the corresponding EEG stages
stages2 = {'Wake','NREM','REM'};
substages = {'Quiet_Wake','Motion_Wake'};

for fil = 1:numel(Files)
    %Load Vm files
    File_Vm = Files{fil};
    [rPath,rFile] = fileparts(File_Vm);
    tmp = load(File_Vm);
    fsVm = tmp.fs;
    Vm = tmp;
    %Define where EEG/Hypnogram files are located
    if contains(fileparts(File_Vm), folder_ID)
        path = strrep(rPath, '5sec_window', '');
    else
        path = rPath;
    end
    %Load EEGcut of files
    if true
        tmp = load(fullfile (path, ['EEG_filt.mat']));
        EEGcut = tmp.EEGf(:,1); %EEG1: 1; EEG2: 2  !!! CORRECT ME PLEASE !!!
        fsEEGcut = tmp.EEG_fs;
    else
        tmp = load(fullfile (path, ['EEG2.mat']));
        EEGcut = tmp.resampled_data_mV(:); %EEG1: 1; EEG2: 2  !!! CORRECT ME PLEASE !!!
        fsEEGcut = tmp.SampRate;
    end
    EEGcut_ts = (1:numel(EEGcut))/fsEEGcut;
    window_EEG = window*fsEEGcut;
    %Load Hypnogram file
    tmp = load(fullfile (path, ['Hypnogram.mat']));
    %Resample Hypnogram to match fsEEGcut
    Hypnogram = repmat(tmp.Hypnogram(:)',fsEEGcut,1);
    Hypnogram_cut = Hypnogram(:);
    Hypnogram_cut(end+1:numel(EEGcut_ts)) = NaN;
    
    %Define Cell IDs as variable: fieldName
    label = strrep (rFile, '_vmStage', '');
    Cell_ID = genvarname(label);
    
    if false
        hf = figure(1);
        ha = NaN(2,1);
        ha(1) = subplot(211);
        plot(EEGcut_ts, Hypnogram_cut,'k','linewidth',2); hold on; ylim([0.5,3.5])
        ha(2) = subplot(212);
        plot(EEGcut_ts, EEGcut,'k'); hold on
        linkaxes(ha,'x')
        xlim([1,EEGcut_ts(end)])
    end
    
    %% PSD analysis for EEG cut across states - Quiet & Motion wake computed from Vm traces data
    for sta = 1:size (stages2, 2)
        stage = stages2{sta};
        
        %Check the indices of bouts in the corresponding Vm analysis
        if isfield (Basic_DATA.(Cell_ID), stage)
            if sta == 1
                for sub = 1:numel (substages)
                    substagez = substages{sub};
                    if isfield (Basic_DATA.(Cell_ID),(substagez)) %%Alter here for new version!!!!
                        ID_traces.(substagez) = Basic_DATA.(Cell_ID).(substagez).VM.trace; %Alter here for new version
                    else
                        ID_traces.(substagez) = [];
                    end %condition notempty DATA-substage-Trace
                end %substagez - WAKE loop
                %join to use to split traces to wake state in the next for loop
                ID_traces.Wake = [ID_traces.Quiet_Wake; ID_traces.Motion_Wake];
            else
                if isfield (Basic_DATA.(Cell_ID),(stage))
                    ID_traces.(stage) = Basic_DATA.(Cell_ID).(stage).VM.trace;
                end %isempty condition DATA-trace
            end
        end %isfield DATA-stage
        
        %Find transitions for stage in the Hypnogram_cut
        indSTA = find([Hypnogram_cut;NaN]==sta & [NaN;Hypnogram_cut]~=sta);
        indEND = find([Hypnogram_cut;NaN]~=sta & [NaN;Hypnogram_cut]==sta)-1;
        %Crop the calculated indexes to match the bouts for the Vm (based
        %on the variable "window" (5 or 10 seconds, depending on the input)
        indSTA_wind = [];
        indEND_wind = [];
        for ind = 1:size (indSTA,1)
            STA = indSTA(ind);
            END = indEND(ind);
            bouts_nr = floor ((indEND (ind) - indSTA (ind))/window_EEG);
            STA_window = NaN (bouts_nr, 1);
            END_window = NaN (bouts_nr, 1);
            %Use a if loop to crop the first STA/END star to window size
            for row = 1:bouts_nr
                STA2 = STA;
                END2 = STA + window_EEG-1;
                STA_window(row) = STA2;
                END_window(row) = END2;
                %Move to the next window
                STA = STA + window_EEG;
            end %while loop
            %Concatenate indexes for the EEGcut window
            indSTA_wind = [indSTA_wind; STA_window];
            indEND_wind = [indEND_wind; END_window];
        end
        %Append data to a EEG_indexes struct
        EEG_indexes.(stage).STA = indSTA_wind;
        EEG_indexes.(stage).END = indEND_wind;
        
        if false %Plot hypnogram and state transitions for quality check control
            stage = 'NREM';
            indSTA = EEG_indexes.(stage).STA;
            indEND = EEG_indexes.(stage).END;
            close all;
            figure (2); ha = NaN(1,2);
            ha(1) = subplot (2,1,1);
            plot (EEGcut_ts, Hypnogram_cut,'marker','.'); hold on
            plot (EEGcut_ts(indSTA),Hypnogram_cut(indSTA), 'go');
            plot (EEGcut_ts(indEND),Hypnogram_cut(indEND), 'ro');
            ylim([0.5,3.5])
            ha(2) = subplot (2,1,2);
            plot (EEGcut_ts, EEGcut);
            linkaxes(ha, 'x')
            xlim([0, EEGcut_ts(end)]);
        end
        
        %Initiate structs to store FFT_EEG
        if isfield (Basic_DATA.(Cell_ID), stage)
            if sta == 1
                if splitWake
                    structSTA = {'Quiet_Wake','Motion_Wake'};
                else
                    structSTA = {'Wake'};
                end %splitWake loop
            else
                 structSTA = {stage};
            end
        else
            structSTA = [];
        end
                
        for s = 1:size (structSTA,2)
            substagez = structSTA{s};
            FFT_EEG.(substagez) = [];
            for bnd = 1:noBND
                [band, lab] = Bands{bnd,:};
                Bpower.EEG.(substagez).(lab) = [];
            end %band loop
        end %condition to load FFT_EEG structs
        %Loop across all the bouts and sort data for each state
        bouts_nr = size (EEG_indexes.(stage).STA,1);
        first_indice = true; %to be used as a flag in the next for loop
        
        for bou = 1:bouts_nr
%             if sta==1 &&  ~ismember (bou, ID_traces.Wake)
%                 continue
%             end
            STA =  EEG_indexes.(stage).STA(bou);
            END =  EEG_indexes.(stage).END(bou);
            bout = ['b',num2str(bou)];
            % Crop EEGcut to Star and End of the selected bout
            EEG_select = EEGcut (STA:END);
            assert(all(Hypnogram_cut(STA:END)==sta), 'grrrhhh');
            %Remove NaNs from EEGcut
            EEG_bout = EEG_select(~isnan(EEG_select));
            if numel (EEG_bout) > window_EEG/2 %Analye EEGcut that contains, at least, 1/2 window of Vm analysis
                %Select EEGcut bouts that are bigger than 2 seconds
                [mpxx_EEG,f_pxx] = pwelch(EEG_bout,[],fsEEGcut*0,frqFFT,fsEEGcut);
                ind = f_pxx>0 & f_pxx<fsEEGcut/2; %two-sided spectra if f is specified
                mpxx_EEG(ind) = 2*mpxx_EEG(ind);
%                 ind = [1 2];
%                 mpxx_EEG(ind) = 0;
            else
                mpxx_EEG = NaN (numel(frqFFT), 1);
            end %condition EEG_bout not empty
            if splitWake
                if sta == 1 && ismember (bou, ID_traces.Quiet_Wake)
                    stagez = 'Quiet_Wake';
                elseif sta == 1 && ismember (bou, ID_traces.Motion_Wake)
                    stagez = 'Motion_Wake';
                else 
                    stagez = stage;
                end %stagez condition
            else
                if sta == 1
                    stagez = 'Wake';
                else
                    stagez = stage;
                end
            end %splitWake loop
            %Bands loop
            if all(~isnan(mpxx_EEG)) && isfield (Basic_DATA.(Cell_ID),(stage))
                for bnd = 1:noBND
                    [band, lab] = Bands{bnd,:};
                    bp_data = bandpower(mpxx_EEG,frqFFT,band,'psd');
                    DATA.(Cell_ID).(stagez).FFT_EEG.Bpower.(lab).(bout) = bp_data;
                    Bpower.EEG.(stagez).(lab) = [Bpower.EEG.(stagez).(lab); bp_data];
                end %band loop
            end
            if isfield (Basic_DATA.(Cell_ID),(stage))
                DATA.(Cell_ID).(stagez).EEG_bout.(bout) = EEG_bout(:);
                DATA.(Cell_ID).(stagez).FFT_EEG.(bout) = mpxx_EEG(:);
                FFT_EEG.(stagez) = [FFT_EEG.(stagez), mpxx_EEG(:)];
                %Store EEG indexes for later analysis
                DATA.(Cell_ID).(stagez).EEG_indexes.indSTA.(bout) = STA;
                DATA.(Cell_ID).(stagez).EEG_indexes.indEND.(bout) = END;
                if first_indice
                    Firstind.(stagez) = STA;
                    first_indice = false;
                end
            end
            if false
                X = EEGcut_ts([STA,STA,END,END]);
                Z = ones(1,4)*-1;
                set(hf,'currentaxes',ha(1))
                ind = strcmpi(stages(:,1),stage);
                if ~any(ind)
                    ind = contains(stages(:,1),stage);
                end
                col = stages{ind,2};
                for k = 1:numel(ha)
                    set(hf,'currentaxes',ha(k))
                    Y = [min(ylim), max(ylim), max(ylim), min(ylim)];
                    h = fill(X,Y,Z,'edgecolor','none','facecolor',col);
                    uistack(h,'bottom')
                end
                drawnow
            end
            
        end %bout loop
    end %stages2 loop
    %Store average FFT across states
    
    if splitWake
        newstages = {'Quiet_Wake','Motion-Wake','NREM','REM'};
    else
        newstages = {'Wake','NREM','REM'};
    end
    
    for sta = 1:size (newstages,2)
        stage = newstages{1,sta};
        if isfield (DATA.(Cell_ID), stage) && isfield (Basic_DATA.(Cell_ID),(stage))    
            DATA.(Cell_ID).(stage).FFT_EEG.average = nanmean (FFT_EEG.(stage),2);
            [data_maxX, ind_maxX] = max (DATA.(Cell_ID).(stage).FFT_EEG.average(1:(find(f_pxx == 20))));
            Crop_f_pxx = f_pxx(6:(find(f_pxx == 20)));
            DATA.(Cell_ID).(stage).FFT_EEG.MaxFrequency = Crop_f_pxx(ind_maxX);
            DATA.(Cell_ID).(stage).FFT_EEG.MaxPower = data_maxX;
            %Band loop for average
            for bnd = 1:noBND
                [band, lab] = Bands{bnd,:};
                bp_average = Bpower.EEG.(stage).(lab);
                DATA.(Cell_ID).(stage).Bpower_EEG.(lab).average = nanmean (bp_average,1);
            end %band loop
        end %isfield loop
    end %stages loop
    
    %% PSD analysis for the Vm trace
    %Loop across sleep/wake stages
    for sta = 1:size (stages2,2)
        stage = stages2{1, sta};
        %Initiate structs to store FFT_Vm
        if isfield (Basic_DATA.(Cell_ID), stage)
            if sta == 1
                if splitWake
                    structSTA = {'Quiet_Wake','Motion_Wake'};
                else
                    structSTA = {'Wake'};
                end %splitWake loop
            else
                structSTA = {stage};
            end %stage loop
        else
            structSTA = [];
        end %isfield Basic_DATA-stage
                
        for s = 1:size (structSTA,2)
            substagez = structSTA{s};
            FFT_Vm.(substagez) = [];
            for bnd = 1:noBND
                [band, lab] = Bands{bnd,:};
                Bpower.Vm.(substagez).(lab) = [];
            end %band loop
        end %condition to load FFT_EEG structs
        
        if ~isempty (Vm.(stage)) && isfield (Basic_DATA.(Cell_ID),(stage))
            %Load Vm data, nr bouts and traces ID to the corresponding stage
            dataVm = Vm.(stage);
            dataVm = dataVm(:);
            nr_bouts = size (dataVm,1);
            traces = ID_traces.(stage);
            %loop across bout but only analyze the traces used for basic Vm
            %analysis
            
            for bou = 1:nr_bouts            
                if any (ismember (traces, bou))
                    bout = ['b',num2str(bou)];
                    Vm_bout = dataVm{bou} - min (dataVm{bou});
                    [mpxx_Vm,f_pxx] = pwelch(Vm_bout,wd_fun(numel(Vm_bout)),2*fsVm,frqFFT,fsVm);
%                     [mpxx_Vm,f_pxx] = pwelch(Vm_bout, [],[],frqFFT,fsVm);
%                     ind = f_pxx>0 & f_pxx<fsVm/2; %two-sided spectra if f is specified
%                     mpxx_Vm(ind) = 2*mpxx_Vm(ind);
%                     ind = [1 2];
%                     mpxx_Vm(ind) = 0;

                    clear stagez
                    if splitWake
                        if sta == 1 && ismember (bou, ID_traces.Quiet_Wake)
                            stagez = 'Quiet_Wake';
                        elseif sta == 1 && ismember (bou, ID_traces.Motion_Wake)
                            stagez = 'Motion_Wake';
                        else
                            stagez = stage;
                        end %stagez condition
                    else
                        if sta == 1
                            stagez = 'Wake';
                        else
                            stagez = stage;
                        end
                    end %splitWake loop
                    
                    %Store Vm bout and FFT
                    if all(~isnan(mpxx_Vm)) && isfield (Basic_DATA.(Cell_ID),(stage))
                        DATA.(Cell_ID).(stagez).Vm_bout.(bout) = Vm_bout;
                        DATA.(Cell_ID).(stagez).FFT_Vm.(bout) = mpxx_Vm;
                        %Store FFT_vm to be averaged after the bout loop
                        FFT_Vm.(stagez) = [FFT_Vm.(stagez), mpxx_Vm(:)];
                        %Bands loop
                        for bnd = 1:noBND
                            [band, lab] = Bands{bnd,:};
                            bp_data = bandpower(mpxx_Vm,frqFFT,band,'psd');
                            DATA.(Cell_ID).(stagez).FFT_Vm.Bpower.(lab).(bout) = bp_data;
                            Bpower.Vm.(stagez).(lab) = [Bpower.Vm.(stagez).(lab); bp_data];
                        end %band loop
                    end %notempty mpxx_VM and Basic_DATA-stage
                end %ismember traces_ID-bout
            end %bouts loop
        end %condition not empty Vm-stage
    end %stage loop
    
    if splitWake
        newstages = {'Quiet_Wake','Motion-Wake','NREM','REM'};
    else
        newstages = {'Wake','NREM','REM'};
    end
    for sta = 1:size (newstages,2)
        stage = newstages{1,sta};
        if isfield (DATA.(Cell_ID), stage)
            DATA.(Cell_ID).(stage).FFT_Vm.average = nanmean (FFT_Vm.(stage),2);
            [data_maxX, ind_maxX] = max (DATA.(Cell_ID).(stage).FFT_Vm.average(6:end));
            Crop_f_pxx = f_pxx(6:end);
            DATA.(Cell_ID).(stage).FFT_Vm.MaxFrequency = Crop_f_pxx(ind_maxX);
            DATA.(Cell_ID).(stage).FFT_Vm.MaxPower = data_maxX;
            %Band loop for average
            for bnd = 1:noBND
                [band, lab] = Bands{bnd,:};
                bp_average = Bpower.Vm.(stage).(lab);
                DATA.(Cell_ID).(stage).Bpower_Vm.(lab).average = nanmean (bp_average,1);
            end %band loop
        end %isfield loop
    end %stages loop
end %files loop
%% Plot PSD of EEG cut and Vm traces for all the cells across states -> with cumulative Theta and Delta history
if true
    fieldNames = fieldnames (DATA);
    ind1 = find (f_pxx == 0);
    ind2 = find (f_pxx == 20);
    
    figure (3)
    sgtitle('PFC: PSD across all the recorded cells');
    ha = [];
    
    if splitWake
        newstages = {'Quiet_Wake','Motion-Wake','NREM','REM'};
    else
        newstages = {'Wake','NREM','REM'};
    end
    
    for sta = 1:size (newstages,2)
        stage = newstages{1,sta};
        ha(end+1) = subplot (5,2,sta);
        for fil = 1:size (fieldNames,1)
            Cell_ID = fieldNames{fil};
            if isfield (DATA.(Cell_ID), (stage))
                if ~isempty (DATA.(Cell_ID).(stage).FFT_EEG.average)
                    data = DATA.(Cell_ID).(stage).FFT_EEG.average;
                    semilogy (f_pxx(ind1:ind2), data(ind1:ind2)); hold on;
                end %not emtpy condition
            end %isfield loop
        end
        title (strrep(stage,'_','\_'));
        xlabel ('Freq. (Hz)')
        ylabel ('PSD (EEG)');
        
        subplot (4,2, sta+4)
        for fil = 1:size (fieldNames,1)
            Cell_ID = fieldNames{fil};
            if isfield (DATA.(Cell_ID), (stage))
                if ~isempty (DATA.(Cell_ID).(stage).FFT_Vm.average)
                    data = DATA.(Cell_ID).(stage).FFT_Vm.average;
                    semilogy (f_pxx(ind1:ind2), data(ind1:ind2)); hold on;
                end %not emtpy condition
            end %isfield loop
        end
        title (strrep(stage,'_','\_'));
        xlabel ('Freq. (Hz)')
        ylabel ('PSD (Vm)');
        
    end    
    
end %condition true/false
linkaxes(ha,'xy')

%% STAT File with averages for all the data
STAT = [];

MainFields = {'FFT_EEG', 'FFT_Vm','Bpower_EEG','Bpower_Vm'};
SubFields = {'MaxFrequency', 'MaxPower','Max_Frequency'};
fieldCells = fieldnames(DATA);

if splitWake
    newstages = {'Quiet_Wake','Motion_Wake','NREM','REM'};
else
    newstages = {'Wake','NREM','REM'};
end

%Pre-allocate STAT struct
for m = 1: length (MainFields)
    mainfield = MainFields{m};
    for sta = 1:size (newstages,2)
        stage = newstages{1,sta};
        if m == 1 || m == 2
            STAT.(stage).(mainfield) = NaN(length(f_pxx), length(fieldCells)-1);
            subfield = SubFields{1};
            STAT.(stage).(subfield).(mainfield) = NaN(1, length(fieldCells)-1);
            subfield = SubFields{2};
            STAT.(stage).(subfield).(mainfield) = NaN(1, length(fieldCells)-1);
       
        elseif m == 3 || m == 4
            noBND = size(Bands,1);
            for bnd = 1:noBND
                [band, lab] = Bands{bnd,:};
                STAT.(stage).(mainfield).(lab) = NaN(1, length(fieldCells)-1);
            end %band loop
        end
    end %stage loop
end %mainfields loop

for sta = 1:size (newstages,2)
    stage = newstages{1,sta};
    
    for m = 1: length (MainFields)
        mainfield = MainFields{m};
        
        if m == 1 || m == 2
            data_average = NaN (length(f_pxx), length (fieldCells));
            data_average2 = NaN (1, length (fieldCells));
            data_average3 = NaN (1, length (fieldCells));
            for cel = 1: length (fieldCells)
                cell_ID = fieldCells{cel};
                if isfield (DATA.(cell_ID), stage)
                    if ~isempty(DATA.(cell_ID).(stage).(mainfield).average)
                        data = DATA.(cell_ID).(stage).(mainfield).average;
                        data_average(:,cel) = data;
                        subfield = SubFields{1};
                        data2 = DATA.(cell_ID).(stage).(mainfield).(subfield);
                        data_average2(:,cel) = data2;
                        subfield = SubFields{2};
                        data3 = DATA.(cell_ID).(stage).(mainfield).(subfield);
                        data_average3(:,cel) = data3;
                    end %notempty loop
                end %condition loop
            end %cell loop
            STAT.(stage).(mainfield) = data_average;
            subfield = SubFields{1};
            STAT.(stage).(subfield).(mainfield) = data_average2;
            subfield = SubFields{2};
            STAT.(stage).(subfield).(mainfield) = data_average3;
            
        elseif m == 3 || m == 4
            noBND = size(Bands,1);
            for bnd = 1:noBND
                [band, lab] = Bands{bnd,:};
                data_average = NaN (1, length (fieldCells));
                for cel = 1: length (fieldCells)
                    cell_ID = fieldCells{cel};
                    if isfield (DATA.(cell_ID), stage)
                        if DATA.(cell_ID).(stage).(mainfield).(lab).average
                            data = DATA.(cell_ID).(stage).(mainfield).(lab).average;
                            data_average(:,cel) = data;
                        end %not empty condition
                    end %condition loop
                end %cell loop
                STAT.(stage).(mainfield).(lab) = data_average;
            end %band loop
        end % if-elseif condition loop
    end %mainfield loop
end % stage loop

STAT.X_Axis.FFT = f_pxx;
%% Plot PSD grand averages
Plot.yEEG = [];
Plot.yVm = [];
groupColors = [];

if true
    x = f_pxx;
    freqRange = [0.5 20];
    figure(4)
    
    if splitWake
        newstages = {... {Stages, Filter settings, RGB combination}
            'Quiet_Wake', [13  114 116]'/255; ...
            'Motion_Wake', [57  73 112]'/255; ...
            'NREM', [234 174  51]'/255;...
            'REM', [218  93  41]'/255; ...
            };
    else
        newstages = {... {Stages, Filter settings, RGB combination}
            'Wake', [57  73 112]'/255; ...
            'NREM', [234 174  51]'/255;...
            'REM', [218  93  41]'/255; ...
            };
    end
    
    for sta = 1:size(newstages,1)
        stage = newstages{sta, 1};
        %Plot PSD EEG
        clear mea err std x_t
        mea = nanmean(STAT.(stage).FFT_EEG,2);
        std = nanstd(STAT.(stage).FFT_EEG,[],2);
        n = size(STAT.(stage).FFT_EEG,2);
        err = std/sqrt(n);
        mea = mea(find(x==freqRange(1)):find(x==freqRange(2)));
        err = err(find(x==freqRange(1)):find(x==freqRange(2)));
        x_t = freqRange(1):0.1:freqRange(2);
        
        subplot (2,2,1)
        plot (x_t, 10*log(mea), 'Color', newstages{sta,2});
        hold on
        curve1 = 10*log(mea + err);
        curve2 = 10*log(mea - err);
        x2 = [x_t, flip(x_t,2)];
        inBetween = [curve1; flip(curve2,1)];
        fill(x2, inBetween',newstages{sta,2}','facealpha',.2,'EdgeColor','none');
        minY1 = min(mea-err);
        maxY1 = max(mea+err);
        title ('PSD EEG');
        xlabel ('Freq. (Hz)')
        ylabel ('PSD');
        
        %Plot PSD Vm
        clear mea err std x_t
        mea = nanmean(STAT.(stage).FFT_Vm,2);
        std = nanstd(STAT.(stage).FFT_Vm,[],2);
        n = size(STAT.(stage).FFT_Vm,2);
        err = std/sqrt(n);
        mea = mea(find(x==freqRange(1)):find(x==freqRange(2)));
        err = err(find(x==freqRange(1)):find(x==freqRange(2)));
        x_t = freqRange(1):0.1:freqRange(2);
        
        subplot (2,2,3)
        plot (x_t, 10*log(mea), 'Color', newstages{sta,2});
        hold on
        curve1 = 10*log(mea + err);
        curve2 = 10*log(mea - err);
        x2 = [x_t, flip(x_t,2)];
        inBetween = [curve1; flip(curve2,1)];
        fill(x2, inBetween',newstages{sta,2}','facealpha',.2,'EdgeColor','none');
        minY1 = min(mea-err);
        maxY1 = max(mea+err);
        title ('PSD Vm');
        xlabel ('Freq. (Hz)')
        ylabel ('PSD');
        
        %Boxplot data FFT_EEG peak & FFT_Vm peak
        ydata = STAT.(stage).MaxFrequency.FFT_EEG(:);
        Plot.yEEG(:,sta) = ydata;
        ydata = STAT.(stage).MaxFrequency.FFT_Vm(:);
        Plot.yVm(:,sta) = ydata;
        groupColors = [groupColors; newstages{sta,2}'];
        
    end %stages loop
    if splitWake
        groupLabels = {'Quiet', 'Motion', 'NREMs', 'REMs'};
    else
        groupLabels = {'Wake', 'NREMs', 'REMs'};
    end %splitwake loop
    subplot (2,2,2)
    boxplot (Plot.yEEG, 'Labels', groupLabels, 'Colors', groupColors)
    title ('PSD EEG');
    ylabel ('Peak Frequency (Hz)');
    
    subplot (2,2,4)
    boxplot (Plot.yVm, 'Labels', groupLabels, 'Colors', groupColors)
    title ('PSD Vm');
    ylabel ('Peak Frequency (Hz)');  
end %true condition

%% Plot FFT/Vm EEG frequencies
Plot.y = [];

conditions = {'Bpower_EEG','Bpower_Vm'};

for c = 1:size(conditions,2)
    condition = conditions{1,c};
    figure (4+c)
    
    for bnd = 1:noBND
        [band, lab] = Bands{bnd,:};
        Plot.y = [];
        subplot (2,3,bnd)
        sgtitle (strrep(condition,'_','\_'));
        
        if splitWake
            newstages = {'Quiet_Wake','Motion_Wake','NREM','REM'};
        else
            newstages = {'Wake','NREM','REM'};
        end
        
        for sta = 1:size(newstages,2)
            stage = newstages{1, sta};
            data = STAT.(stage).(condition).(lab)(:);
            Plot.y(:,sta) = data;
        end %stages loop
        boxplot (Plot.y, 'Labels', groupLabels, 'Colors', groupColors)
        Freqz = sprintf('Freq: %.1fHz - %.1f Hz', band(1), band(2));
        title (Freqz);
        ylabel ('Peak Frequency (Hz)'); 
    end %band loop
end %conditions loop

%% SAVE THE DATA & STAT STRUCT
if savefile
    newDATA = [savecondition 'DATA'];
    newSTAT = [savecondition 'STAT'];
    assignin('base', newDATA, evalin('base', 'DATA'));
    assignin('base', newSTAT, evalin('base', 'STAT'));
%     save (fullSavePath, newDATA);
    save (fullSavePath, newDATA,newSTAT);
end