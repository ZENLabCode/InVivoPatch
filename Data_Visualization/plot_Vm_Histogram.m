% Load Vm and EEG files into a filelist
close all;  clear; clc
addpath(genpath('X:\OptoLab_v4.1\function'));

%Vmstage file containing EEG and corresponding Vm traces
tmp = dir(['X:\1 TIDIS Lab\Tiago\Ephys_PROCESSED\In_vivo\Awake_Sleep\PFC_Processed\Pharmacology\NMDAR_Blockers\*\*\*_vmStage_basicAna_Movie.mat']);
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

%LOAD FILES
%Load Vm file
File_Vm = Files{1};
[rPath,rFile] = fileparts(File_Vm);

%Stages definition
stages = {... {Stages, Filter settings, RGB combination}
    'Quiet_Wake', [13  114 116]'/255; ...
    'Motion_Wake', [57  73 112]'/255; ...
    'NREM', [234 174  51]'/255;...
    'REM', [218  93  41]'/255; ...
    };

%BasicDataAnalysis struct with information on the traces with Quiet/Motion ID
Condition = 'PFC_DAP5_BasicAnalysis_5sec';
savecondition ='DAP5_';
tmp = load (['X:\1 TIDIS Lab\Tiago\Campelo_Manuscript\Fig3_PFC_DAP5\RawData\', Condition]);
Basic_DATA = tmp.(sprintf('%sDATA', savecondition)); 
folder_ID ='5sec_window'; %String to remove from load file directory to load EEG and Hypnogram
window = 5; %BasicAnalysis: 5 seconds windows of Vm

%% Load Vm and sort data to average hystograms for the same cell

%Resample Vm for plotting:
fsVM = 100;
Vm_ts = (1:(5*fsVM))/fsVM;
for fil = 1:numel(Files)
    %Load Vm files
    File_Vm = Files{fil};
    tmp = load(File_Vm);
    [rPath,rFile] = fileparts(File_Vm);
    substringToRemove = "_basicAna_Movie";
    Spikestrign = strrep(rFile, substringToRemove, '');
    tmp2 = load (strcat(rPath, '\', Spikestrign));
    %Define Cell IDs as variable: fieldName
    label = strrep (rFile, '_vmStage_basicAna_Movie', '');
    Cell_ID = ['x',label];
    
    %initiate figure
    figure (fil)
    
    %loop across stages
    for sta = 1:size (stages,1)
        stage = stages{sta,1};
        
        if sta < 3
            filter = 'Wake';
        else
            filter = stage;
        end %use to filter existing traces on the next loop
        
        if ~isempty(tmp.(filter))
            
            if sta == 1 && isfield (Basic_DATA.(Cell_ID), filter)
                VM = tmp.Wake.VM;
                Vmspikes = tmp2.Wake;
                Vmspikes = Vmspikes(:);
                TraceID.Quiet_Wake = Basic_DATA.(Cell_ID).Wake.Quiet_ID;
                TraceID.Motion_Wake = Basic_DATA.(Cell_ID).Wake.Motion_ID;
                
            elseif sta > 2  && isfield (Basic_DATA.(Cell_ID), filter)
                VM = tmp.(stage).VM;
                Vmspikes = tmp2.(stage);
                Vmspikes = Vmspikes(:);
                TraceID.(stage) = Basic_DATA.(Cell_ID).(stage).VM.trace(:);
            end %condition stage
            
            %Save HistX struct
            HistX = VM(1).histX;
            
            %Determine which trace to plot according to the state
            if sta == 1 && isfield (Basic_DATA.(Cell_ID), filter)
                tracez = TraceID.Quiet_Wake;
            elseif sta == 2 && isfield (Basic_DATA.(Cell_ID), filter)
                tracez = TraceID.Motion_Wake;
            elseif sta > 2 && isfield (Basic_DATA.(Cell_ID), filter)
                tracez = TraceID.(stage);
            end
            
            %Initiate Average struct for histogram of distributions
            Average_HistY = [];
            
            if isfield (Basic_DATA.(Cell_ID), filter)
                for tra = 1:size(tracez,1)
                    trace = tracez (tra);
                    
                    if ~isnan (trace)
                        if tra < 4
                            subplot (4,4,(tra - 1) * 4 + sta)
                            n = numel(Vmspikes{trace,1});
                            VM_noSAM = floor(n/tmp2.fs*fsVM);
                            VM_data = interp1((1:n)/tmp2.fs, Vmspikes{trace,1},(1:VM_noSAM)/fsVM);
                            plot (Vm_ts, VM_data, 'Color', stages{sta,2});
                        end
                        
                        HistY = VM(trace).histY;
                        Average_HistY = [Average_HistY, HistY];
                    end
                    
                    %         Average_HistY = Average_HistY(:,1);
                    if ~isempty(Average_HistY)
                        Average_HistY = nanmean (Average_HistY(:,1:1),2);
                        
                        IndPlot = find (HistX < -40);
                        
                        subplot(4, 4, 3 * 4 + sta)
                        b = bar (HistX(IndPlot), Average_HistY(IndPlot), 'FaceColor', stages{sta,2}, 'FaceAlpha', 0.3); hold on;
                        b.BarWidth = 1; % Adjust the bar width to remove space between bars
                        b.EdgeColor = 'none'; % Remove the black lines around the bars
                        plot(b.XData, smooth(b.YData,4), 'Color', stages{sta,2}, 'linewidth', 1)
                        xlabel('Bins (mV)');
                        ylim([0 20]);
                        ylabel('Counts (%)');
                    end %is not emtpy Average_HistY
                end %traces loop
            end
        end %notempty data-filter
    end %stages loop
end %file loop

