clear%% Defined directory for all "_vmStage_basicAna.mat" files
close all; clc; clear;
addpath(genpath('X:\OptoLab_v4.1\function'));

tmp = dir(['X:\1 TIDIS Lab\Tiago\Ephys_PROCESSED\In_vivo\Awake_Sleep\PFC_Processed\Pharmacology\NMDAR_Blockers\**\*_basicAna_10sec.mat']);

%Analysis stages
stages = {'Wake','NREM','REM'};
stages2 = {... {Stages, Filter settings, RGB combination}
    'Quiet_Wake', [13  114 116]'/255; ...
    'Motion_Wake', [57  73 112]'/255; ...
    'NREM', [234 174  51]'/255;...
    'REM', [218  93  41]'/255; ...
    };

%Saving options true or false
zave ='true';

%Alter value of max cumsum according to the window of Vm analysis
% cumsum_size = 25000; %5 sec window
cumsum_size = 50000; %10 sec window

%SAVE OUTPUT
savePath = 'X:\1 TIDIS Lab\Tiago\Campelo_Manuscript\FigX_PFC_DAP5\RawData';
saveFile = "DAP5_BasicAnalysis_10sec.mat";
condition = 'DAP5_';
fullSavePath = fullfile(savePath, saveFile);

tmp([tmp.isdir]) = [];

if isempty(tmp)
    fprintf(2,'No files found for: %s\n',files.vmStage);
    return
end

Files_BasicAna  = selectionList(fullfile({tmp.folder},{tmp.name}));
if isempty(Files_BasicAna)
    fprintf(2,'No file selected\n');
    return
end
cell_IDs = {tmp.name};
rPaths = cellfun(@(x)fileparts(x),Files_BasicAna,'Uniformoutput',false);
Files_Vmstage = cellfun(@(x)ls(fullfile(x,'*_vmStage_10sec.mat')),rPaths,'Uniformoutput',false);
Files_Vmstage = fullfile(rPaths,Files_Vmstage);

%% Extract, normalize the "basic analysis" data
mainfield = {'VM', 'pupil','motion','EMG'};
subfield_VM = {'trace','mean','std','prctile15','prctile95_subtracted','cumsumY'};
subfield_Video = {'mean','std'};
subfield_EMG = {'totalPower'};

for fil = 1:numel(Files_Vmstage)
    BasicAna = Files_BasicAna{fil};
    Vmstage = Files_Vmstage{fil};
    tmp = load(BasicAna);
    tmp_Vmstage = load(Vmstage);
    
    % identify cellID as fieldName
    ID = cell_IDs{1,fil};
    numberStr = erase(ID, '_basicAna_10sec.mat');
    fieldName = ['f' numberStr];

    %Define an initial filtering to remove traces that might have passed
    %visual inspection
    for sta = 1:numel(stages)
        stage = stages{sta};
        if ~isempty (tmp.(stage))
            %Determine Zscores of the cumsum to avoid bouts that are outliers
            main = 1;
            field = mainfield{1};
            sub = 6;
            subfield = subfield_VM{sub};
            data_cumsum = [];
            
            for tra = 1:numel ({tmp.(stage).VM(:).trace})
                data = cell2mat({tmp.(stage).(field)(tra).(subfield)}');
                if numel (data) >= cumsum_size
                    data = data(cumsum_size);
                    data_cumsum = [data_cumsum; data];
                else
                    data = NaN;
                    data_cumsum = [data_cumsum; data];
                end %filter out traces with lower cumsum
            end % loop across traces
            mean = nanmean (data_cumsum);
            std = nanstd (data_cumsum);
            Zscores_cumsum = (data_cumsum - mean) / std;
            %Compute the indexes to be analyzed
            Zscores_indexes.(stage) = find (Zscores_cumsum < 2 & Zscores_cumsum > -2);
        end %not empty tmp-stage condition
    end %stages loop
    
    %Start to loop across all the data to extract the DATA
    for sta = 1:numel(stages)
        stage = stages{sta};
        if isfield (Zscores_indexes, stage)
            select_indexes = Zscores_indexes.(stage);
        end %condition isfield Zscores-stage
        %Initiate data structs
        Zscores_cumsum = [];
        Zscores_motion = [];
        Zscores_EMG = [];
        
        % Visual inspection of the NREM and REM traces
        if sta == 2 || sta == 3
            if ~isempty (tmp_Vmstage.(stage))
                %Create figure and determine the number of subplots
                close all; figure('Name','Visual inspection of traces: NREM/REM','windowstate','maximized'); drawnow;
                %Make sure data is organized in columns
                tmp_stage = tmp_Vmstage.(stage)(:);               
                n   = min([numel(select_indexes),50]);                    
                fac = 1;
                
                while n/fac> 20 %normal: 10
                    fac = fac+1;
                end
                sub1 = floor(n/fac);
                sub2 = ceil(n/sub1);
                cut_off = 75; %normal: 50
                
                %loop across Quiet_ID traces
                for tra = 1:n
                    if numel(select_indexes) <= cut_off   
                        trace = select_indexes(tra);
                    elseif numel(select_indexes) >= cut_off
                        traz = numel(select_indexes) - cut_off + tra;
                        trace = select_indexes(traz);
                    end % trace loop to select the last ones
                    Vm_trace = tmp_Vmstage.(stage){trace};
                    subplot (sub1, sub2, tra);
                    plot (((1:numel(Vm_trace))/tmp_Vmstage.fs), Vm_trace);
                    xlabel('Time (sec.)');
                    ylabel('Vm (mV)');
                    titleStr = sprintf('Trace #%d', trace);
                    title(titleStr);
                end %Quiet_ID loop
                
                % Popup input dialog to get EEGmin, EEGmax, EMGmin, EMGmax, CumSuMmin, and CumSuMmax
                prompt = {'Select traces to keep (use (,) to filter traces)'};
                dlgtitle = 'Traces to keep';
                dims = [1 50]; % Adjusted dimensions for more inputs
                definput = {'1'}; % Default inputs for each field
                answer = inputdlg(prompt, dlgtitle, dims, definput);
                
                if ~isempty(answer)
                    string = strsplit (answer{1,1}, ',');
                    select_indexes = str2double(string);
                    select_indexes = select_indexes(:);
                end
            end %is not empty codition
        end %traces plot NREM/REM
        
        for main = 1:numel(mainfield)
            field = mainfield{main};
            data = [];
            
            if ~isempty (tmp.(stage))
                if any(~isnan(select_indexes))
                    if main == 1
                        for sub = 1:numel (subfield_VM)
                            subfield = subfield_VM{sub};
                            if sub < 6
                                data = cell2mat({tmp.(stage).(field).(subfield)}');
                                data = data (select_indexes);
                                DATA.(fieldName).(stage).(field).(subfield) = data;
                            elseif sub == 6
                                data_cumsum = [];
                                cumsum_trace = [];
                                traces = DATA.(fieldName).(stage).VM.trace;
                                for tra = 1:size ({tmp.(stage).VM(:).trace},2)
                                    trace_ID = tmp.(stage).VM(tra).trace;
                                    if ismember (trace_ID, traces)
                                        if ismember (trace_ID, select_indexes)
                                            data = cell2mat({tmp.(stage).(field)(tra).(subfield)}');
                                            cumsum_trace(:,tra) =  data(:);
                                            data = data(cumsum_size);
                                            data_cumsum = [data_cumsum; data];
                                        end %ismember trace_ID - select indexes
                                    end
                                end %traces loop  
                                DATA.(fieldName).(stage).(field).(subfield) = data_cumsum;
                                if sta > 1
                                    DATA.(fieldName).(stage).(field).averageCumsum = nanmean(cumsum_trace,2);
                                else
                                    DATA.(fieldName).(stage).(field).averageCumsum = cumsum_trace;
                                end %condition stage - cumsumtrace
                                Zscores_cumsum = [Zscores_cumsum; data_cumsum];
                            end %subfield filtering loop
                        end %subfieldVM loop
                        
                    elseif main == 2 || main == 3
                        if isfield (tmp.(stage),'pupil') || isfield (tmp.(stage),'motion')
                            if ~isempty (tmp.(stage).(field))
                                for sub = 1:numel (subfield_Video)
                                    subfield = subfield_Video{sub};
                                    data = cell2mat({tmp.(stage).(field).(subfield)}');
                                    data = data (select_indexes);
                                    DATA.(fieldName).(stage).(field).(subfield) = data;
                                    if main == 3 & sub == 2
                                        Zscores_motion = [Zscores_motion; data];
                                    end %concatenate Zscores of video motions
                                end % %subfield_Video
                            end %not empty pupil/motion
                        end %isfield condition to check if video was processed
                    elseif main == 4
                        for sub = 1:numel (subfield_EMG)
                            subfield = subfield_EMG{sub};
                            data = cell2mat({tmp.(stage).(field).(subfield)}');
                            data = data (select_indexes);
                            DATA.(fieldName).(stage).(field).(subfield) = data;
                            Zscores_EMG = [Zscores_EMG; data];
                        end %loop subfield_EMG
                    end %condition subfield for mainfield
                end %notnan condition
            end %stage not empty condition
        end %mainfield loop
        Zscores.(stage).cumsum = zscore(Zscores_cumsum);
        Zscores.(stage).motion = zscore(Zscores_motion);
        Zscores.(stage).EMG = zscore (Zscores_EMG);
        
        %Define thresholding for each individual recorded cell
        % Plot Z-scores as a scatter plot
        if sta == 1
            if ~isempty(Zscores.(stage).motion)
                figure()
                scatter(Zscores.('Wake').motion, Zscores.(stage).cumsum);
                title ('ZscoreMotion vs ZscoreCumsum: Quiet vs Wake motion');
                xlabel('ZscoreMotion');
                ylabel('ZscoreCumsum');
            end %compute only when motion Zscore contains data.
            
            figure()
            scatter(Zscores.(stage).EMG, Zscores.(stage).cumsum);
            title ('ZscoreMotion vs ZscoreCumsum: Quiet vs Wake motion');
            xlabel('ZscoreEMG');
            ylabel('Zscorecumsum');
            
            % Popup input dialog to get EEGmin, EEGmax, EMGmin, EMGmax, CumSuMmin, and CumSuMmax
            prompt = {'Enter MotionQuiet:', 'Enter Motion_Motion:', 'Enter EMGQuiet:', 'Enter EMGMotion:', 'Enter CumSumQuiet:', 'Enter CumSumMotion:'};
            dlgtitle = 'Input for EEG and EMG Ranges';
            dims = [1 50]; % Adjusted dimensions for more inputs
            definput = {'0.1', '0.1', '0.2', '0.2', '0', '0'}; % Default inputs for each field
            answer = inputdlg(prompt, dlgtitle, dims, definput);
            
            % Check if user provided inputs
            if ~isempty(answer)
                % Convert input strings to numbers
                MotionQuiet.(stage) = str2double(answer{1});
                Motion_Motion.(stage) = str2double(answer{2});
                EMGQuiet.(stage) = str2double(answer{3});
                EMGMotion.(stage) = str2double(answer{4});
                CumSumQuiet.(stage) = str2double(answer{5});
                CumSumMotion.(stage) = str2double(answer{6});
            else
                disp('User cancelled the input dialog.');
            end %loop if answer was provided
        end %sta
    end %stage loop
    
    %Split Wake Traces into Quiet & Motion Wake
    Quiet_ID = [];
    Motion_ID = [];
    
    stage = 'Wake';
    if isfield (DATA.(fieldName),'Wake')
        
        traces_Wake = DATA.(fieldName).Wake.VM.trace;
        
        for tra = 1:numel(traces_Wake)
            id = traces_Wake(tra);
            Zcumsum = Zscores.Wake.cumsum(tra);
            Zemg = Zscores.Wake.EMG(tra);
            
            if isempty (Zscores.Wake.motion)
                Zmotion = NaN;
            else
                Zmotion = Zscores.Wake.motion(tra);
            end
            
            if ~isnan(Zmotion)
                if Zcumsum < CumSumQuiet.(stage) & Zemg < EMGQuiet.(stage) & Zmotion < MotionQuiet.(stage)
                    Quiet_ID = [Quiet_ID; id];
                    
                elseif Zcumsum > CumSumMotion.(stage) & Zemg > EMGMotion.(stage) & Zmotion > Motion_Motion.(stage)
                    Motion_ID = [Motion_ID; id];
                end
            else
                if Zcumsum < CumSumQuiet.(stage) && Zemg < EMGQuiet.(stage)
                    Quiet_ID = [Quiet_ID; id];
                    
                elseif Zcumsum > CumSumMotion.(stage) & Zemg > EMGMotion.(stage)
                    Motion_ID = [Motion_ID; id];
                end
            end
        end
        %% Visual inspection of the Quiet/Wake traces
        %Create figure and determine the number of subplots
        close all; figure('Name','Visual inspection of Quiet Wake Traces','windowstate','maximized'); drawnow;
        %Make sure data is organized in columns
        tmp_Vmstage.Wake = tmp_Vmstage.Wake(:);
        n   = min([numel(Quiet_ID),50]);
        fac = 1;
        
        while n/fac>10
            fac =fac+1;
        end
        sub1 = floor(n/fac);
        sub2 = ceil(n/sub1);

        %loop across Quiet_ID traces
        for tra = 1:n
            trace = Quiet_ID(tra);
            Vm_trace = tmp_Vmstage.Wake{trace,1};
            subplot (sub1, sub2, tra);
            plot (((1:numel(Vm_trace))/tmp_Vmstage.fs), Vm_trace);
            xlabel('Time (sec.)');
            ylabel('Vm (mV)');
            title ('Visual inspection of traces: Quiet Wake');
            titleStr = sprintf('Trace #%d', trace);
            title(titleStr);
        end %Quiet_ID loop
        
        % Popup input dialog to get EEGmin, EEGmax, EMGmin, EMGmax, CumSuMmin, and CumSuMmax
        prompt = {'Select traces to keep (use (,) to filter traces)'};
        dlgtitle = 'Traces to keep';
        dims = [1 50]; % Adjusted dimensions for more inputs
        definput = {'1'}; % Default inputs for each field
        answer = inputdlg(prompt, dlgtitle, dims, definput);
        
         if ~isempty(answer)
             string = strsplit (answer{1,1}, ',');
             Quiet_ID = str2double(string);
             Quiet_ID = Quiet_ID(:);
         end
         
         %Plot Motion_ID traces
         pause (1)
         close all; figure('Name','Visual inspection of Motion Wake Traces','windowstate','maximized'); drawnow;
         n   = min([numel(Motion_ID),50]);
         fac = 1;
         
         while n/fac>10
             fac =fac+1;
         end
         sub1 = floor(n/fac);
         sub2 = ceil(n/sub1);
         
         %loop across Quiet_ID traces
         for tra = 1: n
             trace = Motion_ID(tra);
             Vm_trace = tmp_Vmstage.Wake{trace,1};
             subplot (sub1, sub2, tra);
             plot (((1:numel(Vm_trace))/tmp_Vmstage.fs), Vm_trace);
             title ('Visual inspection of traces: Motion Wake');
             xlabel('Time (sec.)');
             ylabel('Vm (mV)');
             titleStr = sprintf('Trace #%d', trace);
             title(titleStr);
         end %Quiet_ID loop
         
         % Popup input dialog to get EEGmin, EEGmax, EMGmin, EMGmax, CumSuMmin, and CumSuMmax
         prompt = {'Select traces to keep (use (,) to filter traces)'};
         dlgtitle = 'Traces to keep';
         dims = [1 50]; % Adjusted dimensions for more inputs
         definput = {'1'}; % Default inputs for each field
         answer = inputdlg(prompt, dlgtitle, dims, definput);
         
         if ~isempty(answer)
             string = strsplit (answer{1,1}, ',');
             Motion_ID = str2double(string);
             Motion_ID = Motion_ID(:);
         end
         close all;
         
         %% Assigned Waketraces to DATA struct
         subfield_VM = {'trace','mean','std','prctile15','prctile95_subtracted','cumsumY','averageCumsum'};
         for main = 1:numel (mainfield)
             field = mainfield{main};
             DATA.(fieldName).Wake.Motion_ID = Motion_ID;
             DATA.(fieldName).Wake.Quiet_ID = Quiet_ID;
             
             if main == 1
                 for sub = 1:numel (subfield_VM)
                     subfield = subfield_VM{sub};
                     
                     if ~isempty (Quiet_ID) && all(~isnan (Quiet_ID))
                         data_Quiet = [];
                         data = [];
                         for q = 1:numel(Quiet_ID)
                             quiet_trace = Quiet_ID(q);
                             trace_id = find (DATA.(fieldName).Wake.VM.trace == quiet_trace);
                             if sub < 7
                                 data = DATA.(fieldName).Wake.(field).(subfield)(trace_id);
                                 data_Quiet = [data_Quiet; data];
                             else
                                 data = DATA.(fieldName).Wake.(field).(subfield)(:,trace_id);
                                 data_Quiet(:,q) = data;
                             end %use only traces that are not empty
                         end %subfield loop
                     end %loop Quiet_ID
                     DATA.(fieldName).Quiet_Wake.(field).(subfield) = nanmean(data_Quiet,2);
                 end %condition - notempy Quiet_ID
                 
                 if ~isempty (Motion_ID) && all(~isnan (Motion_ID))
                     data_Wake = [];
                     data = [];
                     for w = 1:numel(Motion_ID)
                         motion_trace = Motion_ID(w);
                         trace_id = find (DATA.(fieldName).Wake.VM.trace == motion_trace);
                         if sub < 7
                             data = DATA.(fieldName).Wake.(field).(subfield)(trace_id);
                             data_Wake = [data_Wake; data];
                         else
                             data = DATA.(fieldName).Wake.(field).(subfield)(:,trace_id);
                             data_Wake(:,q) = data;
                         end %subfield loop
                     end %loop Quiet_ID
                     
                     DATA.(fieldName).Motion_Wake.(field).(subfield) = nanmean (data_Wake,2);
                 end %condition not empty Motion_ID
                 
             elseif main == 2 || main == 3
                 if isfield (DATA.(fieldName).Wake,'pupil') || isfield (DATA.(fieldName).Wake,'motion')

                     for sub = 1:numel (subfield_Video)
                         subfield = subfield_Video{sub};
                         
                         if ~isempty (Quiet_ID) && all(~isnan (Quiet_ID))
                             data_Quiet = [];
                             data = [];
                             
                             for q = 1:numel(Quiet_ID)
                                 quiet_trace = Quiet_ID(q);
                                 trace_id = find (DATA.(fieldName).Wake.VM.trace == quiet_trace);
                                 data = DATA.(fieldName).Wake.(field).(subfield)(trace_id);
                                 data_Quiet = [data_Quiet; data];
                             end %loop Quiet_ID
                             DATA.(fieldName).Quiet_Wake.(field).(subfield) =  nanmean(data_Quiet,2);
                         end %condition - notempy Quiet_ID
                         
                         if ~isempty (Motion_ID) && all(~isnan (Motion_ID))
                             data_Wake = [];
                             data = [];
                             
                             for w = 1:numel(Motion_ID)
                                 motion_trace = Motion_ID(w);
                                 trace_id = find (DATA.(fieldName).Wake.VM.trace == motion_trace);
                                 data = DATA.(fieldName).Wake.(field).(subfield)(trace_id);
                                 data_Wake = [data_Wake; data];
                             end %loop Motion_ID
                             DATA.(fieldName).Motion_Wake.(field).(subfield) = nanmean (data_Wake,2);
                         end %condition not empty Motion_ID
                     end %subfield_video loop
                 end %isfield DATA-Pupil & DATA-motion
                 
                 elseif main == 4
                     for sub = 1:numel (subfield_EMG)
                         subfield = subfield_EMG{sub};
                         
                         if  ~isempty (Quiet_ID) && all(~isnan (Quiet_ID))
                             data = [];
                             data_Quiet = [];
                             
                             for q = 1:numel(Quiet_ID)
                                 quiet_trace = Quiet_ID(q);
                                 trace_id = find (DATA.(fieldName).Wake.VM.trace == quiet_trace);
                                 data = DATA.(fieldName).Wake.(field).(subfield)(trace_id);
                                 data_Quiet = [data_Quiet; data];
                             end %loop Quiet_ID
                             DATA.(fieldName).Quiet_Wake.(field).(subfield) = data_Quiet;
                         end %condition - notempy Quiet_ID
                         
                         if ~isempty (Motion_ID) && all(~isnan (Motion_ID))
                             data = [];
                             data_Wake = [];
                             
                             for w = 1:numel(Motion_ID)
                                 motion_trace = Motion_ID(w);
                                 trace_id = find (DATA.(fieldName).Wake.VM.trace == motion_trace);
                                 data = DATA.(fieldName).Wake.(field).(subfield)(trace_id);
                                 data_Wake = [data_Wake; data];
                             end %loop Motion_ID
                             DATA.(fieldName).Motion_Wake.(field).(subfield) = data_Wake;
                         end %condition not empty Motion_ID
                     end %subfield loop
             end %mainfields condition
         end %mainfield loop
    end %isfield DATA
end %file loop

%% Create STAT file
STAT = [];
% subfield_VM = {'mean','std','prctile15','prctile95_subtracted','averageCumsum','cumsumY'}; %remove traces from average struct
subfield_VM = {'mean','std','prctile15','prctile95_subtracted','cumsumY'}; %remove traces from average struct
cell_list = fieldnames (DATA);

for sta = 1: size (stages2, 1)
    stage = stages2 {sta};
    
    for main = 1:numel(mainfield)
        mainfiel = mainfield{main};
        
        if main == 1
            subfieldz = subfield_VM;
        elseif main == 2 || main == 3
            subfieldz = subfield_Video;
        else
            subfieldz = subfield_EMG;
        end %main-to-subfield loop
        
        for sub = 1:numel (subfieldz)
            subfield = subfieldz{sub};
            data_cells = [];
            
            for fil = 1:numel(cell_list)
                fieldName = cell_list {fil};
                
                if isfield (DATA.(fieldName),(stage))
                    if isfield (DATA.(fieldName).(stage),(mainfiel))
                        if isfield (DATA.(fieldName).(stage).(mainfiel),(subfield))
                            if main == 1 && sub == 6 % UPDATE HERE IF ALTERING SUBFIELD VM
                                data = DATA.(fieldName).(stage).(mainfiel).(subfield);
                                data_cells(:,fil) = data;
                            else
                                data = nanmean (DATA.(fieldName).(stage).(mainfiel).(subfield));
                                data_cells = [data_cells; data];
                            end %subfield average loop
                        end %isfield DATA-subfield
                    end %isfield DATA-mainfield
                end %isfield DATA-stage
            end %cell/fill loop
            
            %Fix the issue of initiating the struct file with an empty cell
            %(where data is = 0)
            notZeroColumns = ~all(data_cells == 0, 1);
            data_cells = data_cells(:, notZeroColumns);
            STAT.(stage).(mainfiel).(subfield) = data_cells;
        end %subfield loop
    end %mainfield loop
end %stages loop
%% Plot the results for all the conditions
goforit = true;
groupColors  = [];
groupLabels = {'Quiet', 'Motion', 'NREMs', 'REMs'};

if false
    figure (1);
    for sub = 1:size(subfield_VM,2)
        substagez = subfield_VM{1, sub};
        subplot (2,4, sub);
        Plot  = [];
        
        %Loop data across all stages
        for sta = 1:size (stages2, 1)
            stage = stages2{sta,1};
            %Boxplot data FFT_EEG peak & FFT_Vm peak
            ydata = STAT.(stage).VM.(substagez);
            Plot(:,sta) = ydata;
            
            if goforit
                groupColors = [groupColors, stages2{sta,2}'];
            end
        end %stage loop
        goforit = false;
        
        boxplot (Plot, 'Labels', groupLabels)
        
    end %subfield loop
end %true/false loop

%% SAVE THE DATA & STAT STRUCT
if zave
    newDATA = [condition 'DATA'];
    newSTAT = [condition 'STAT'];
    assignin('base', newDATA, evalin('base', 'DATA'));
    assignin('base', newSTAT, evalin('base', 'STAT'));
    save (fullSavePath, newDATA,newSTAT);
end