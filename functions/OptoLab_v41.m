% OptoLab 4.1
%------------------------------------------------------------------------
% A MATLAB toolbox for optogenetics/sleep studies
%
%   Developed by
%      Mojtaba Bandarabadi, PhD
%      Email: m.bandarabadi@gmail.com
%      Zentrum f?r experimentelle Neurologie (ZEN)
%      Inselspital, Bern, 2015-2017
%
%   Updated on by
%       Thomas Rusterholz, 2018
%       Update from OptoLab 4
%
%------------------------------------------------------------------------

%% ---------------------------- Initiation ---------------------------- %%
OptoDir = com.mathworks.mlservices.MLEditorServices.getEditorApplication;
OptoDir = char(OptoDir.getActiveEditor.getStorageLocation.getFile);
OptoDir = fileparts(OptoDir);
addpath(genpath(OptoDir))
%TR: re-initiation each section with init_OptoLab();


%% ------------------------- Channel selection ------------------------ %%
%% Channel selection
clear; close all; clc
dform = '*.dat';  % '*.ncs','*.dat','*.mat'
chanSel = 'sel';  % 'all','sel','list'
switch chanSel
    case 'all'
        chanNames = dir(dform);
        chanNames = natsort({chanNames.name})
    case 'sel'
        chanNames = cellstr(uigetfile(dform,'Select','multiselect','on'))
    case 'list'
        tmp = sprintf('amp-A-%03i.dat;',[0:5,7,26:27]);
        chanNames = regexp(tmp(1:end-1),';','split');
end
clear chanSel dform
%% -------------------------------------------------------------------- %%


%% -------------------------- Converting ------------------------------ %%
%% Neuralynx *.ncs to *.mat
init_OptoLab(); clc;
rPath='C:\Dexmedetomidine\N17 Dexmedetomidine\'; %N11, 12 & 14:17
tmp=dir(fullfile(rPath,'*.ncs')); cd(rPath);
chanNames={tmp.name};
SampRate = 1000;   % sampling rate
portNumber = [];   % stim port number, or empty ([])
f_Nlx2Mat(chanNames,SampRate,portNumber)

%% Intan *.dat to *.mat
init_OptoLab(); clc;
% cd 'T:\Projects\test_Data\from_datFiles'
% cd 'Z:\Thomas\testData_OptoLab5\Intan_Carlos\'
% cd 'C:\Projects\testData\INTAN'
%EEG and EMG channels, e.g. {'amp-A-017','amp-A-018'}
idxEEG = {}; %renaming to EEG (or EEG1, EEG2, ...)
idxEMG = {... {EMG channels, REF channels}; %re-referencing, filt & save
    'amp-A-008','amp-A-010';... PS: no re-reference if no 2nd column2!
    'amp-A-016';...
    };

idxEMG = {};
%sampling rate in [Hz] (resampling if not empty)
SampRate = 1000; %PS: recorded SampRate must be dividable by this
%read Intan rhd-file (TR, new script, returns data in a structure)
dataIntan = read_Intan_RHD2000_file_zenlab();
%export data
res = f_Intan2Mat(dataIntan,SampRate,idxEEG,idxEMG); %TR, new script

% %save channel info
% xlswrite('channel_info.xlsx',res)
% xls_cellFit('channel_info.xlsx');

%% Matlab *.mat to *.bin for SlipAnalysis
init_OptoLab(); clc;
% chanNamesEEG = {'EEG1','EEG2'}; % EEG channels (e.g. {'EEG1','EEG2',...})
% chanNamesEMG = {'EMG'};   % EMG channels (e.g. {'EMG1','EMG2',...})
tmp=sprintf('amp-C-%03i.dat;',1:31);
chanNamesEEG=regexp(tmp(1:end-1),';','split');
chanNamesEMG={};
%TR - new function (no downsampling)
f_Mat2Bin_noDS(chanNamesEEG,chanNamesEMG,true)
clear

%% Binary *.bin to *.mat (A-M systems recordings)
% THIS SECTION NOW WORKS WITH FULL FILE PATH
init_OptoLab(); clc;
clc; clear; close all
... cd C:\Projects\testData\AMSystems
    ... cd C:\Projects\testData\fromMary
    ... cd C:\Projects\test_Data\fromCaro2
    ... cd 'Z:\Thomas\forCaro\stoke file for spectrum analysis'
    ... cd 'Z:\2 AstroSleep_Group\Irene\EEG_EMG\Sleep Anlysis\CMT stroke\Sleep Characterization\SHAM\CMT15\BL\D'
    ... cd 'Z:\1 TIDIS Lab\Yudong\LH project\EEG_LH\Scored\Vgat\Vgat LH 1-Scored\VgatLH1_Session_1'
    
%RECORDING INFO (do manually as some exp-files have wrong/missing info)
%get_binFileInfo() or get_binFileInfo(file), exp- or ini-file)
rec=get_binFileInfo(fullfile(pwd,'*.exp'));
check_binFile0(rec)%to check & correct manually
%fList=matlab.codetools.requiredFilesAndProducts('check_binFile0')
% select file
rec.binFile   = ls('*.bin'); %one file only!!!
%rec.binFile   = 'Z:\1 TIDIS Lab\Yudong\LH project\EEG_LH\Scored\POA\POA2\POA2_DZPsession_2';
rec.binPath   = pwd; %path of bin file
%tmp = dir('Z:\1 TIDIS Lab\Yudong\LH project\EEG_LH\Scored\POA\POA2\POA2_DZPsession_2\POA2_DZPsession_2_2022-12-15_14-47-07.bin');
%files = fullfile({tmp.folder},{tmp.name});
%files = 'Z:\1 TIDIS Lab\Yudong\LH project\EEG_LH\Scored\POA\POA2\POA2_DZPsession_2_2022-12-15_14-47-07.bin';

%% just to search files only once
tmp = dir('C:\Users\I0344502\Desktop\Paula\**\*.bin');
files = fullfile({tmp.folder},{tmp.name});
file = selectionList(files,[],struct('title','Select ONE file!!!'));
while numel(file)>1
    file = selectionList(files,[],struct('title','Select ONE file!!!'));
end
if numel(file)==0
    fprintf(2,'No File Selected!\n')
    return
end%% CONVERTING
init_OptoLab(); clc;
% path = 'X:\2 AstroSleep_Group\Astro Stroke\MA group\MA GROUP';
%path = 'Z:\1 TIDIS Lab\Yudong\LH project\EEG_LH\Scored\POA\POA2';
%addpath('Z:\OptoLab_v4.1\function\convert')
%cd(path)
%tmp = dir(fullfile(path,'**','*.bin'));
% tmp = dir('E:\OXT40_Pene_S1.bin');
% namelist = {tmp.name};
% pathslist = {tmp.folder};
% rPaths = fullfile(pathslist,namelist);
% rPaths=selectionList(rPaths);
rPaths = file;
for k = 1:numel(rPaths)
    [rec.binPath,rec.binFile,tmp] = fileparts(rPaths{k});
    rec.binFile = [rec.binFile,tmp];
    rec.precision = 'uint16'; %saved precision ('int16' or 'uint16')
    rec.SampRate  = 512; %recorded sampling rate [Hz]
    rec.channels  = {'EEG1','EEG2', 'EMG'}; %ALL channels, this order
%     rec.channels  = {'EEG1','EEG2','EEG3','EEG4','EMG'}; %ALL channels, this order
    %rec.channels  = {'EEG','EMG'}; %ALL channels, this order
    rec.chanSTIM  = {}; %{'Stim'}; %'Stim', stim channels (out of rec.channels)
    rec.acquisitionRangeMax = 10;
    rec.offset    = 5;
    rec.minStimGap = 0.02*0.5; %treats stims as one stimulus when gap is < this, in seconds
    %calibration function (find correct one, or ask someone and make test-plot)
    rec.bin2data=@(x) (x/2^15*rec.acquisitionRangeMax - rec.offset) *1000;
    %  *1000 for [V] --> [mV]
    % if offset==0, then probably -->
    % rec.bin2data=@(x) (x/2^16*rec.acquisitionRangeMax) *1000;
    % plot and check!!
    
    %CONVERT
    SampRate = rec.SampRate/512*1000; %wished sampling rate (resamples)
    testPlot = true; %true or false, to make test-plot (do at least for 1st file)
    f_Bin2Mat(rec,SampRate,testPlot);
%     f_Bin2Mat_BRIDGE(rec,SampRate,testPlot);
end

%Test plot
% test_plotConverted(rec.channels);

%% EEGLAB, *.mat to *.set & *.fdt (for Armand's sleep scoring tool)
init_OptoLab(); clc;
% cd 'C:\Projects\testData\INTAN\data_inDAT'
% cd 'C:\Projects\testData\INTAN\data_inRHD'
% cd 'C:\Projects\testData\INTAN\data_inRHD2'

chanNames={'EEG1','EEG2','EMG'};
chanLabels={}; %channel labels, e.b. {'EEG','EMG'} (auto if empty)
%   optional: a montage file string, e.g. chanLabels='montage.xyz'
savename='EEGLABdata'; %savename (exports *.set and *.fdt)
f_mat2eeglab(chanNames,savename,chanLabels)

%%%TO TEST (SCORING)
%% res = f_sleepScoring_old('EEGLABdata_old',1);
Stages={... stages to score {num, label} num as number to score keyboard
    1,'Wake';...
    2,'REM';...
    3,'NREM';...      ,'N'   ,'NREM';...
    };
rFile='EEGLABdata_old';
res = f_sleepScoring(rFile,Stages);
%save
hypnogram=res.hypnogram;
SampRate=res.srate; %for having sample rate of hypnogram
save(fullfile(fileparts(rFile),'Hypnogram.mat'),'hypnogram','SampRate')
%test
% a=load('Hypnogram.mat')

%% save('Hypnogram.mat','-struct','res')


%% -------------------------------------------------------------------- %%




%% ------------------------ Hypnogram functions ----------------------- %%
%% Correct Hypnogram and convert to mat format
init_OptoLab(); clc;
f_hypnoCorrect()   % f_hypnoCorrect('4stages') for 4 stages scoring

%% Sleep quantification (percentage, number, duration, and distance) and
%  transposing hypnogram
init_OptoLab(); clc;
winLen = 3600;     % window size in sec to calculate segmented statistics
resHypno = f_hypnoStatistics(winLen);
clear winLen
%% Scoring (csc-eeg-tools from Armand Mensen)
init_OptoLab(); clc;
% cd 'C:\Projects\testData\INTAN\data_inRHD2'
file='EEGLABdata.set';
res = f_sleepScoring(file);
save('Hypnogram.mat','-struct','res')




%% ------------ Average raw traces / spectrograms over stim ----------- %%
%% Average raw traces over stim
init_OptoLab(); clc;
clc; close all; %PS: hynogram epochs must be 1s
cd 'T:\Projects\test_Data\fromCaro3';
chanNames={'EEG'}; %,'EMG.mat','EOG.mat'}


numPulse = 1;      % number of stim pulses during each trial
margin = [10,10];  % time before and after each trial [s]
dsFactor = 10;     % downsampling factor
stageSel = 'probability'; % select stage based on 'probability' or 'start'
%NEW SCRIPT
resMeanRaw = f_rawStim(chanNames,numPulse,margin,stageSel,dsFactor);
save('resMeanRaw.mat','resMeanRaw')

%%%TEST
% tic
resMeanRawO = f_rawStim_old20180712(chanNames,numPulse,margin,stageSel,dsFactor);
% xlim/1000
% toc
% isequal(resMeanRaw,resMeanRawO)


%% Average spectrograms over stim pulses
init_OptoLab(); clc;
clc; close all; %PS: hynogram epochs must be 1s
% cd 'C:\Projects\test_Data\fromCaro1';
cd 'T:\Projects\test_Data\fromCaro3';
chanNames={'EOG','EEG','EMG'}; %,'EMG.mat','EOG.mat'}
chanFilter={[0,1],NaN,NaN}; %bandpass, no filter if NaN
specNames=chanNames(2); %calculate/plot spectra (out of chanNames)
numPulse = 1;      % number of stim pulses during each trial
margin = [10,10];  % time before and after in sec
minfreq = 0.5;     % min frequency for plot
maxfreq = 20;      % max frequency for plot
yscale = 'log';    % Y axis scale; 'log' or 'lin'
idStim = 1;        % index of laser
stageSel = 'probability'; % select stage based on 'probability' or 'start'
%NEW SCRIPT ,SUB_SCRIPT f_spectrogram.m
resSpec=f_spectrogramStim_new(chanNames,specNames,chanFilter,numPulse,...
    margin,minfreq,maxfreq,yscale,idStim,stageSel);

% resSpec = f_spectrogramStim_old20180712(chanNames,numPulse,margin,minfreq,maxfreq,...
%     yscale,idStim,stageSel);


%% New script (remove somewhere else)
init_OptoLab(); clc;
cd 'Z:\Thomas\forCaro\stoke file for spectrum analysis'
clc; close all;
opt.rPath=pwd; %read path
opt.chanNames={'EEG2'}; %
opt.Stages={1,'Wake';2,'NREM';3,'REM'}; %{number, label}
opt.SampRateHyp=1; %Sampling rate of hypnogram (psd per stage)

% f_psd_newForCaro(opt)
% f_psd_newForCaro2(opt)
% f_psd_newForCaro3(opt)


%%


%% Average spectrograms over two laser probes
init_OptoLab(); clc;
numPulse = 5;      % number of stim pulses during each trial
margin = [10,10];  % time before and after in sec
minfreq = 0.5;     % min frequency for plot
maxfreq = 20;      % max frequency for plot
yscale = 'log';    % Y axis scale; 'log' or 'lin'
resSpec2 = f_spectrogramStim2(chanNames,numPulse,margin,minfreq,...
    maxfreq,yscale);
%% -------------------------------------------------------------------- %%


%% -------------------------- Power spectrum -------------------------- %%
%% Stage specific power spectral density (PSD)
init_OptoLab(); clc;
% cd('Z:\Thomas\testData\Caro\27.08.2018D')
% chanNames={'amp-C-020.dat','amp-C-016.dat'};
cd C:\Projects\testData\CarlosCoherence
chanNames={'amp-D-008','amp-D-012'};
clc; close all;
opt.SampRateHyp=1000; %Sampling rate of saved hypnogram (for upsampling)!!!
opt.timeIndex=[]; %all if empty, interactive selection if NaN
[resPSD,timeIndex]=f_psd(chanNames,opt);
export_PSD(resPSD,'',timeIndex)

%%% same for coherence (Carlos)

%% PSD of SSFO experiments
init_OptoLab(); clc;
stages = '3stages'; % '3stages' (W,NR,R), '4stages' (aW,qW,NR,R)
resPSdSSFO = f_psdSSFO(chanNames,stages);

%% PSD before/during/after stim
init_OptoLab(); clc;
margin = 10;        % margin before/after stim in sec
mrgPl = 1;          % merge too close pulses (0 or 1)
resPSdStim = f_psdStim(chanNames,margin,mrgPl);

%% Power in a specific sleep stage over time
init_OptoLab(); clc;
clc;
winLen = 300;       % window size in sec
stepLen = 150;      % moving step in sec
stage = 2;          % sleep stage index
mrgPl = 0;          % merge too close pulses (0 or 1)
%TR, NEW SCRIPT
resPSdTime = f_psdTime(chanNames,winLen,stepLen,stage,mrgPl);
%TR, NEW SCRIPT
f_psdTime_plot(resPSdTime)

%% Stage specific PSD before, while and after stimuli train
cd C:\Projects\testData\CarlosCoherence

init_OptoLab(); clc;
chanNames={'amp-D-008','amp-D-012'};
opt.margin=5; %in [s], before and after stimulus train!
opt.stimVariable='stimTimes'; %e.g 'stimTimes', 'stimeTimes2' ,...
opt.SampRateHyp=1000; %Sampling rate of saved hypnogram!!!
opt.minTrainGap=1; %minimal gap [s] between stim trains (to find trains)
opt.Stages={1,'Wake';2,'NREM';3,'REM'}; %{number, label}
opt.noise.freq=50:50:500; %for 50 or 60 Hz and harmonics removal
opt.noise.df=0.5; %freq +- df will be replaced by following data
%psd
res=f_psdByStim(chanNames,opt);
%plot (2nd property structure, 3rd function handle (e.g. @(x)smooth(x,20))
props.axis={'ylim',[10^-3,10^3]};
plot_psdByStim(res,props,@(x)smooth(x,20));


%% -------------------------------------------------------------------- %%


%% ------------------- Coherence ---------------------------------------%%
%% before and after stimuli train start
cd C:\Projects\testData\CarlosCoherence

init_OptoLab(); clc;
chanNames={'amp-D-008','amp-D-012'};
opt.margin=5; %in [s], before and after stimulus train!
opt.stimVariable='stimTimes'; %e.g 'stimTimes', 'stimeTimes2' ,...
opt.f=0.5:0.5:40; %frequency output vector
opt.SampRateHyp=1000; %Sampling rate of saved hypnogram!!!
opt.minTrainGap=1; %minimal gap [s] between stim trains (to find trains)
opt.Stages={1,'Wake';2,'NREM';3,'REM'}; %{number, label}
%coherence
res = f_coherenceByStim(chanNames,opt);
%plot
plot_coherenceByStim(res);

%% Stage specific coherence
cd C:\Projects\testData\CarlosCoherence

init_OptoLab(); clc;
chanNames={'amp-D-008','amp-D-012'};
opt.f=0.5:0.5:40; %frequency output vector
opt.SampRateHyp=1000; %Sampling rate of saved hypnogram!!!
opt.Stages={1,'Wake';2,'NREM';3,'REM'}; %{number, label}
%coherence
res = f_coherence(chanNames,opt);
%plot
s.name='coherence'; %savename without extension (no saving if empty)
s.print={'-depsc','-r300'}; %variablxes for command print
plot_coherence(res,s);



%% ------------------------- Spike analysis --------------------------- %%
cd 'T:\Projects\test_Data\spike_Ivan\BE52slash17_2\baseline'; %spikes
chanNames={'amp-A-021.dat'} %%% tmp TR

%% ------------ Spike/LFP per pulse trail -----------------------------%%

%% plot spikes per pulse trail
init_OptoLab(); clc;
cd('Z:\Thomas\testData\Mary\Test'); close all;
chanNums={0,1}; %channel number (as cell)
chanNums={0}; %channel number (as cell)

data.chanNames=cellfun(@(x) sprintf('amp-A-%03i',x),...
    chanNums,'uniformoutput',false); %channel files
data.scoringLength=1; %scoring length hypnogram [s] (for upsampling)
data.margin=[10,10]; %time before/after each stimulus pulse train [s]
data.minPulseTrainInterval=1.1; %in [s] (to separate pulse trains)
% use slightly larger than max within trains, to compensate inaccuracy

%plot pulse trains +- margin, return indieces pulse trains (cell)
res = f_spikesPerPulseTrain(data); %res is cell
sname='pulseTrainIndex.xlsx';
xlswrite(sname,res); %save as xls-file

%%% POWER, %f_powerPerPulseTrain(chanNames,res of filename)
res2=f_powerPerPulseTrain(data.chanNames,sname);

%%% SPIKES COUNT, f_spikesCountPerPulseTrains(chanNames,res of filename)
res3=f_spikesCountPerPulseTrains(data.chanNames,sname);
sname2='pulseTrain_spikeCount.xlsx';
delete(sname2)
for cha=1:numel(data.chanNames)
    channel=data.chanNames{cha};
    xlswrite(sname2,res3(:,:,cha),channel)
end
xls_deleteSheets(sname2,data.chanNames,-1)
xls_cellFit(sname2)


%% Detect and sort single channel spikes
init_OptoLab(); clc;
chanNames={'amp-A-000.dat'} %including extension!!!
timeVect = [];
%saves spike_FILENAME.mat
f_spikeDetection(chanNames,timeVect) % extract spikes
%apply data to spike_FILENAME.mat, saves data_wc_amp-A-021.run
f_spikeSorting(chanNames)            % sort detected spikes


%% Raster plot of units
init_OptoLab(); clc;
idxPlot = [1 1000]; % start and end in sec
winLen = 500;       % window size for spike rate in ms
resSpikeRate = f_spikeRaster(chanNames,idxPlot,winLen); %NEW SCRIPT

%% Average spike rate over stimulations
init_OptoLab(); clc;
numPulse = 5;       % number of stim pulses during each trial
binLen = 200;       % bin length for averaging spikes in ms
margin = 2000;      % margin from stim in ms
stageSel = 'start'; % select stage based on 'probability' or 'start'
%TR NEW SCRIPT (current data has no stimTimes, find other to test)
resSpikeStim = f_spikeStim(chanNames,margin,numPulse,binLen,stageSel);

%% Average spikes over transitions (requires scoring)
init_OptoLab(); clc;
margin = 5000;      % margin from transition in ms
binLen = 50;       % bin length for averaging spikes in ms
transition = 'N2W'; % transition type:'W2N','N2R','N2W','R2W'
%TR, should work
resSpikeTrn = f_spikeTransition(chanNames,margin,transition,binLen);

%% Statistics of spikes over sleep/wake cycles (requires scoring)
init_OptoLab(); clc;
%TR, should work
SampRate=1000;
timeIndex=[1,30*60*SampRate]; %e.g. [1,30*60*SampRate] 1st 30 min of rec
resSpikeStages = f_spikeStage(chanNames,timeIndex);
%% -------------------------------------------------------------------- %%



%% ------------------------- Spindle analysis ------------------------- %%
%% Detection for mice
%% cd 'C:\Projects\test_Data\fromCaro2';
init_OptoLab(); clc;
cd C:\Projects\testData\fromMary
chanNames={'EEG'}; %TR
idxNREM = 2;  % index of NREM in Hypnogram
for cha=1:numel(chanNames)
    channel=chanNames{cha};
    f_spindleDetectionMice(channel,'',sprintf('spindle_%s',channel),...
        idxNREM)
end
%% Plot detected spindles
init_OptoLab(); clc;
% histograms:'hist'; plot all: 'all'; plot one: 'one'
chanNames = {'EEG1'};
f_spindlePlot(chanNames,'all');

%% Spindle statistics (rate, REM transition)
init_OptoLab(); clc;
clc; close all
resSpinStat = f_spindleStatistics(chanNames);

%% Raster plot of spindles
init_OptoLab(); clc;
freqBand = [];  % freq. band to calculate power, or leave empty []
resSpinRaster = f_spindleRaster(chanNames,freqBand);

%% Average of spectrograms, slow waves and spindle envelop
init_OptoLab(); clc;
clc; close all; chanNames={'EEG1','EEG2','EEG3','EEG4'} %temp TR
stage = 2;       % sleep stage index
margin = [2,2];  % margin ([before after]) in sec
%TR, NEW SCRIPT
% resSpinSpectSWS = f_spindleSpectrogramSWS(chanNames,stage,margin);
resSpinSpectSWS = f_spindleSpectrogramSWS_2(chanNames,stage,margin);
%% -------------------------------------------------------------------- %%

%% ------------------------- Ripple analysis ------------------------- %%
%% Ripple Detection (NOT YET FINISHED!!!)
init_OptoLab(); clc;
rPaths={'T:\Projects\test_Data\spike_Ivan\BE52slash17_2\baseline'};
chanNames={'amp-A-021'};

% rPaths={'Z:\Thomas\RIPPLETEST'};
% chanNames={'amp-A-004'}; %!!!!! neu, dat-file, nicht konvertiert

% rPaths=pwd; %one path only (current path)
% chanNames={'amp-A-021'};

%ripple detection settings
opt.fpass=[140,250]; %passpand [Hz]
opt.minRIPdur=10; %min ripple duration [ms]
opt.minGAPdur=0; %min gap duration [ms] (joins ripples with smaller gap!)
opt.funThreshold=@(x) mean(x)+4*std(x); %threshold function
f_rippleDetection(rPaths,chanNames,opt,[]); %saves ripples_CHANNAME.mat
%f_rippleDetection_2(chanNames,fpass,fileHYP) --> not better

%test plot
filename=fullfile(rPaths{1},sprintf('ripples_%s.mat',chanNames{1}));
test_plotRipples(filename)

%% events analyses per pulse train or ripples, spindles and up/down states
init_OptoLab(); clc; clear;
% cd 'T:\Projects\test_Data\spike_Ivan\BE52slash17_2\baseline'; %temp
% chanNames={'amp-A-021'};
% cd 'C:\Projects\testData\CarlosCoherence'
% chanNames={'amp-D-008'};
% cd  'Z:\Thomas\testData\Mary\Test'
% chanNames={'amp-A-00'};
cd 'Z:\Ivan\files_and_data\preliminary_data\BE52slash17_2\baseline'
chanNames={'amp-A-021'};

cd 'Z:\Carlos\Reuniens Project\Opto\ReuOpt_Batch1\ReuOpt_3_e_181213_151934'
chanNames={'amp-A-008'};


%input arguments as structure
opt.readPath=''; %current path if empty
opt.chanNames=chanNames;
opt.margin=2.5; %time before/after each stimulus pulse train [s]
opt.minTrainGap=1; %minimal gap between trains [s], for to separate trains
%Must not be exact, just > maximal stimuli gap within
%trains and <= minimal train gap (smaller is even
%better as timing is mostly not exact)
%E.g. if gaps between stimuli within trains is ~0.1s
%and train gaps ~10s, minTrainGap=1 works fine!
opt.stimVariable='stimTimes'; %e.g 'stimTimes', 'stimeTimes2' ,...
opt.testPlot=false; %plot to check pulse trains


% opt.SampRateHyp=1000; %Sampling rate of saved hypnogram!!!
% opt.Stages={1,'Wake';2,'NREM';3,'REM'}; %{number, label}
% opt.noise.freq=50:50:500; %for 50 or 60 Hz and harmonics removal
% opt.noise.df=0.5; %freq +- df will be replaced by following data


%
res = f_eventsPerPulseTrain(opt);

%%


%% ------------------------ Up/Down analysis -------------------------- %%
%% Up/Down state detection
init_OptoLab(); clc;
idxNREM = 2;     % index of NREM sleep
freqBand = [];   % frequency band to calculate power
resUpDown = f_upDownDetection(chanNames,idxNREM,freqBand);

%% Plot detected Up/Down states
init_OptoLab(); clc;
idxNREM = 2;     % index of NREM sleep
f_upDownPlot(chanNames,idxNREM)

%% Average spike rate of up/down states. The first channel in the list is
%  taken as the reference
init_OptoLab(); clc;
binLen = 10;         % bin length for averaging spikes in ms
margin = 500;        % margin from start in ms
state = 'downStart'; % transition: 'upStart','upEnd','downStart','downEnd'
idStim = []';     % inside:'in1','in2'; outside:'out1','out2'; all: []
resSpikeUpDown = f_upDownSpike(chanNames,margin,binLen,state,idStim);
%% -------------------------------------------------------------------- %%



%% ----------- Estimate coherency between a set of channels ----------- %%
% Coherency between channels in a frequency band before/during/after stim
init_OptoLab(); clc;
margin = 10;         % margin before/after stim in sec
freqBand = [0.5 4];  % frequency of interest to measure coherency
mrgPl = 1;           % merge to close pulses (0 or 1)
resSyncStim = f_syncStim(chanNames,margin,freqBand,mrgPl);
%% -------------------------------------------------------------------- %%
