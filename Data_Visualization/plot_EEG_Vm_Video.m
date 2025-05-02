% plot_representation
%-------------------------------------------------------------------------
% Plots axes:
%   - Hypnogram
%   - VM data (downsampled)
%   - EMG
% Highlights selected stage bouts (they have margins)

%% ADD:
%  - spectra (color coded)
%    e.g. 5 s windows
%  - eye tracking

clc; clear; %close all
addpath(genpath(fullfile('X:\OptoLab_v4.1\function')))

%PARAMETERS
%----------
%FILES
%defined VM- & vidoe -file
id    = '21o21004';
ino   = 'IN0';
rPath = 'X:\1 TIDIS Lab\Tiago\Ephys_PROCESSED\In_vivo\Awake_Sleep\RC_Processed\20211022_TC_P038_Day1\';
fileVID = 'X:\1 TIDIS Lab\Tiago\Pupil_Tracking_RAW\In_vivo_Patch\20211021_TC_P038_Day1\Image 2021_10_21_124630,884_processed_proc.mat';
%other files
subPath = sprintf('%s_%s',id,ino);
Files = {...{label, full filename}; when changing labels, change code
    ... PS: when chaning label, must also be changed in code
    'VM'   , fullfile(rPath,sprintf('%s_%s.mat',id,ino));...
    'Bouts', fullfile(rPath,subPath,sprintf('%s_vmStage.mat',id));...
    'HYP'  , fullfile(rPath,subPath,'Hypnogram.mat');...
    'EMG'  , fullfile(rPath,subPath,'EMG.mat');...
    'EEG'  , fullfile(rPath,subPath,'EEG2.mat');...
    'Video', fileVID;...
    };

%OPTIONS
%stages (case sensitive, fieldname in *_vmStage.mat)
Stages = {... {stages, color}
    'NREM', [0.9290 0.6940  0.1250];...
    'REM' , [0.8500 0.3250  0.0980];...
    'Wake', [0      0.4470  0.7410];...
    };
%VM sample rate (fs)
% - original is 100000 Hz, which is huge (interpolates data)
%   PS: EMG has 1000, bout index 5000, hypnogram will be upsampled to EMG
% - fsVM must be <= original
fsVM = 100; %[Hz], must be set, normally 5000
fs   = 45; %all the rest (if not empty)
%Filter function
% Dependent on data, freq-band and sampling rate, @(data,fBand,fs)
% E.g. - filt.fband.EMG = [10,50];
%        filt.fun = @(data,fs)passband_fourier(data,filt.fband.EMG,fs)
%        Passband filter 10-50 Hz
filt.fun  = @(data,fBand,fs)passband_fourier(data,fBand,fs);
filt.fband.EMG = [10 250]; %[Hz] (no filter if empty) & [10 80]
filt.fband.EEG = [1 40]; %[Hz] (no filter if empty) [.9 40]
filt.fband.VM  = [0,fsVM/3]; %upper should be < fsVM/2 (only if fsVM < 100000)
%spectra
spec.flim = [0,30]; %foa all, set [0,inf];
spec.cmap = jet; %colormap

%PLOT PROPERTIES
props.fig = {};
%axes (uses default axes size, might be changed by figure position)
props.axi = {};
props.fig_createAxes = {[70,50,50],[45,45,45],'pixel'}; %dx, dy & unit
%plot
props.plotHYP = {'color',ones(1,3)*.3};
props.plotEEG = {'color',ones(1,3)*.3};
props.plotEMG = {'color',ones(1,3)*.3};
props.plotVM  = {'color',ones(1,3)*.3};
props.plotSTA = {}; %colors see above
%for plotting stages
plotSTA.VM  = true;
plotSTA.EEG = true;
plotSTA.EMG = true;
plotSTA.HYP = false;

%MAIN SCRIPT
%------------
scriptName = mfilename;
fprintf('%s\n%s\n',scriptName,repmat('-',size(scriptName)))

%INIT
%all files must exist
ind = find(cellfun(@(x)exist(x,'file')~=2,Files(:,2)));
if ~isempty(ind)
    fprintf(2,'Missing Files:\n')
    for k = 1:numel(ind)
        %[~,rPath,rExt] = fileparts(Files{ind(k),2});
        fprintf(2,'  %s - %s\n',Files{ind(k),:});
    end
    return
end
%number of ...
noSTA = size(Stages,1);

%% LOAD DATA
fprintf('Load Data\n')
warning('off','MATLAB:dispatcher:UnresolvedFunctionHandle') %annoying
for fil = 1:size(Files,1)
    [label,file] = Files{fil,:};
    [~,rPath,rExt] = fileparts(file);
    fprintf('  - %s (%s)\n',[rPath,rExt],label)
    D.(label) = load(file);
    %resample
    if ~(isempty(fs) || isnan(fs))
        switch lower(label)
            case {'eeg','emg'}
                D.(label).resampled_data_mV = resample(...
                    D.(label).resampled_data_mV, fs, D.(label).SampRate);
                D.(label).SampRate = fs;
                fprintf('    resampled to %g Hz\n',fs)
            case 'bouts'
                fields = fieldnames(D.(label).info.indBouts);
                for k = 1:numel(fields)
                    field = fields{k};
                    D.(label).info.indBouts.(field) = round(...
                        D.(label).info.indBouts.(field)/D.(label).fs*fs);
                end
                D.(label).fs = fs;
                fprintf('    resampled to %g Hz\n',fs)
            case {'vm','hyp','video'}
                %nothin to do, e.g. hypnogram upsampled later
            otherwise
                disp(D.(label))
                error('Upsampling not implemented')
        end
    end
end
%fprintf('YES\n')
%return%%%

%% PREPARA DATA
fprintf(' Prepare Data\n');
%bouts
Bouts.fs = D.Bouts.fs;
%for uncutted, Bouts.indEEG = D.Bouts.info.indEEG;
Bouts.indEEG = 1:numel(D.EEG.resampled_data_mV)/D.EEG.SampRate*D.Bouts.fs;
Bouts.indBOU = D.Bouts.info.indBouts;
%VM
VM.fs    = fsVM;
VM.data  = D.VM.resampled_data_mV(:);
VM.noSAM = numel(VM.data);
if D.VM.SampRate>VM.fs
    %check
    if VM.fs/2<filt.fband.VM(2)
        error('fsVM/2 is smaller than filt.fband.VM(2) (%g < %g)',...
            VM.fs/2,filt.fband.VM(2))
    end
    %filter
    fprintf('  - VM passband filter, %g - %g Hz\n',...
        min(filt.fband.VM(:)),max(filt.fband.VM(:)));
    VM.data = filt.fun(VM.data,filt.fband.VM,VM.fs);
    %interpolation
    fprintf('  - VM interpolation, %i --> %i Hz\n',...
        D.VM.SampRate,VM.fs);
    n = VM.noSAM;
    VM.noSAM = floor(n/D.VM.SampRate*VM.fs);
    VM.data = interp1((1:n)/D.VM.SampRate,VM.data,(1:VM.noSAM)/VM.fs);
end
%EMG
EMG.fs    = D.EMG.SampRate;
EMG.data  = D.EMG.resampled_data_mV(:);
EMG.noSAM = numel(EMG.data);
if ~isempty(filt.fband.EMG)
    fprintf('  - EMG passband filter, %g - %g Hz\n',...
        min(filt.fband.EMG(:)),max(filt.fband.EMG(:)));
    EMG.data = filt.fun(EMG.data,filt.fband.EMG,EMG.fs);
end
%EEG
EEG.fs    = D.EEG.SampRate;
EEG.raw   = D.EEG.resampled_data_mV(:);
EEG.data  = EEG.raw;
EEG.noSAM = numel(EEG.data);
if ~isempty(filt.fband.EEG)
    fprintf('  - EEG passband filter, %g - %g Hz\n',...
        min(filt.fband.EEG(:)),max(filt.fband.EEG(:)));
    EEG.data = filt.fun(EEG.data,filt.fband.EEG,EEG.fs);
end
%Hypnogram (same fs as tight above)
if EEG.fs~=EMG.fs
    error('uups')
end
HYP.fs    = EEG.fs;
HYP.noSAM = EEG.noSAM;
tmp = repmat(D.HYP.Hypnogram(:)',EEG.fs,1);
HYP.data  = tmp(:);
HYP.data(numel(HYP.data)+1:HYP.noSAM) = NaN;

% fac = floor(HYP.noSAM/numel(HYP.data));
% if fac<1
%     error('Downsampling Hypnogram is not implemented')
% elseif fac>0
%     fprintf('  - Hypnogram upsampling by factor %i\n',fac);
%     HYP.data = repmat(HYP.data',fac,1);
%     HYP.data = HYP.data(:);
% end
% HYP.data(numel(HYP.data)+1:HYP.noSAM) = NaN;

%% VIDEO
fname = Files{strcmpi(Files(:,1),'video'),2};
[rPath,rFile,rExt] = fileparts(fname);
fprintf('  - %s\n',[rFile,rExt])
tmp = load(fname,'pupil','filenames','motSVD_1');
%movie (needed ???)
fname = squeeze([tmp.filenames])';
[~,rFile,rExt] = fileparts(fname);
if exist(fname,'file')~=2 %if was moved to another folder
    fname = fullfile(rPath,[rFile,rExt]);
end
fprintf('  - %s\n',[rFile,rExt])
%append
file = squeeze([D.Video.filenames])';
if exist(file,'file')==2
    VID.movie  = VideoReader(file);
    VID.fs     = VID.movie.FrameRate;
else
    VID.fs     = 10;
    warning('Video file not found. Frame Rate manually set to %g',VID.fs)
end
VID.pupil  = D.Video.pupil{1};
VID.pupil  = double(VID.pupil.area_smooth(:));
VID.pupil  = VID.pupil/max(VID.pupil);
VID.motion = double(mean(abs(D.Video.motSVD_1),2));
VID.t      = (1:numel(VID.pupil))/VID.fs;


%% CALCULATE EEG SPECTRA
N  = 4*EEG.fs; %N = 5 for 5-s window4
steps = N/2; %hope it's an integer ;-), else return error anyway
ind = round(Bouts.indEEG/Bouts.fs*EEG.fs);
tmp = ind(1)-steps:steps:ind(end)-steps; %index time points
IND = bsxfun(@plus,repmat(tmp,N,1),(1:N)');
ind = IND>0 & IND<=EEG.noSAM;
DAT = NaN(size(IND));
DAT(:,2:end-1) = EEG.raw(IND(:,2:end-1));
ind = IND(IND(:,1)>0);
DAT(ind,1) = EEG.raw(ind);
ind = IND(IND(:,end)>=EEG.noSAM);
DAT(ind,end) = EEG.raw(ind);
% IND(IND<1|IND>EEG.noSAM) = NaN;
% DAT = EEG.raw(IND);
%init
[rows,cols] = size(DAT);
SPEC.fs  = EEG.fs/steps;
SPEC.t   = (0:cols-1)/SPEC.fs;
[~,SPEC.f] = pwelch(randn(rows,1),hann(rows),0,[],EEG.fs);
SPEC.PXX = NaN(numel(SPEC.f),cols);
%spectra
indN = any(isnan(DAT),1);
for k = find(indN)
    tmp = DAT(:,k);
    idx = ~isnan(tmp);
    SPEC.PXX(:,k) = pwelch(tmp(idx),hann(sum(idx)),0,SPEC.f,EEG.fs);
end
SPEC.PXX(:,~indN) = pwelch(DAT(:,~indN),hann(rows),0,SPEC.f,EEG.fs);

%% FIGURE/AXES
% clc; close all

[dx,dy] = props.fig_createAxes{1:2};
len = 880; hig = 80; %axis length, hight
noAXI = 7;
pos = [10, 10, [1,0,1]*dx(:)+len, [1,(noAXI-1),1]*dy(:)+noAXI*hig];
hf = figure('position',pos,props.fig{:});
movegui(hf,'center'); drawnow
ha = fig_createAxes(hf,[noAXI,1],props.fig_createAxes{:});
set(ha,'unit','normalized'); drawnow
axi = 0; %axes indes
%bout index
indBOU = max([1, round(min(Bouts.indEEG)/Bouts.fs*EMG.fs)]) : ...
    min([EMG.noSAM, round(max(Bouts.indEEG)/Bouts.fs*EMG.fs)]);
t = ((1:EMG.noSAM)-indBOU(1)+1)/EMG.fs; %[s], EEG, EMG, 111
tMax = t(indBOU(end));

%PLOT Hypnogram
axi = axi+1;
set(hf,'CurrentAxes',ha(axi));
plot(t(indBOU),HYP.data(indBOU),props.plotHYP{:},'tag','HYP'); hold on;
set(gca,'ylim',[0.7,3.2],'ytick',1:3,...
    'yticklabel',{'Wake','NREM','REM'})
title({strrep(id,'_','\_'),'Hypnogram'})
set(gca,'tag','HYP');
% 
%PLOT VM
axi = axi+1;
set(hf,'CurrentAxes',ha(axi));
tVM = (1:VM.noSAM)/VM.fs; %[s]
tMax = max([tVM(end),tMax]);
plot(tVM,VM.data,props.plotVM{:},'tag','VM'); hold on
title('VM')
ylabel('[mV]')
set(gca,'tag','VM');

%PLOT EMG
axi = axi+1;
set(hf,'CurrentAxes',ha(axi));
y = EMG.data(indBOU);
plot(t(indBOU),y,props.plotEMG{:},'tag','EMG'); hold on;
if isempty(filt.fband.EMG)
    title('EMG')
else
    title(sprintf('EMG (passband %g - %g Hz)',...
        min(filt.fband.EMG(:)),max(filt.fband.EMG(:))))
end
% xlabel('Time [s]')
ylabel('[mV]')
set(gca,'tag','EMG','ylim',[-1,1]*max(abs(y))*1.05)

%PLOT EEG
axi = axi+1;
set(hf,'CurrentAxes',ha(axi));
y = EEG.data(indBOU);
plot(t(indBOU),y,props.plotEEG{:},'tag','EEG'); hold on;
if isempty(filt.fband.EEG)
    title('EEG')
else
    title(sprintf('EEG (passband %g - %g Hz)',...
        min(filt.fband.EEG(:)),max(filt.fband.EEG(:))))
end
% xlabel('Time [s]')
ylabel('[mV]')
set(gca,'tag','EEG','ylim',[-1,1]*max(abs(y))*1.05)

%PLOT SPECTRA
axi = axi+1;
set(hf,'CurrentAxes',ha(axi));
indF = SPEC.f>=spec.flim(1) & SPEC.f<=spec.flim(2);
pcolor(SPEC.t,SPEC.f(indF),10*log10(SPEC.PXX(indF,:)))
caxis([-50,-30]) %spec limit for colors
shading interp %flat looks equal to imagesc
%imagesc(SPEC.t,SPEC.f(indF),10*log10(SPEC.PXX(indF,:)))
axis xy;
% xlabel('Time [s]')
ylabel('Frequency [Hz]')
title('EEG Power Spectra [dB]')
colormap(spec.cmap);
set(gca,'tag','SPEC')

%PLOT MOTION
axi = axi+1;
set(hf,'CurrentAxes',ha(axi));
plot(VID.t,VID.motion,'tag','motion')
% xlabel('Time [s]')
title('Motion')
set(gca,'tag','motion')

%PLOT PUPIL
axi = axi+1;
set(hf,'CurrentAxes',ha(axi));
plot(VID.t,VID.pupil,'tag','pupil')
xlabel('Time [s]')
title('Pupil')
set(gca,'tag','pupil')
tMax = max([VID.t(end),tMax]);

%% PLOT STAGE BOUTS
labels = fieldnames(plotSTA);
noLAB  = numel(labels);
for lab = 1:noLAB
    label = labels{lab};
    if ~plotSTA.(label)
        continue
    end
    set(hf,'CurrentAxes',findobj(ha,'type','axes','tag',label));
    tmp = findobj(gca,'type','line','tag',label);
    t     = tmp.XData;
    data  = tmp.YData;
    fs    = 1/mean(diff(t));
    noSAM = numel(t);
    
    %for each stages
    hp = NaN(1,noSTA);
    for sta = 1:noSTA
        [stage,color] = Stages{sta,:};
        %get indices
        IND = Bouts.indBOU.(stage)-Bouts.indEEG(1)+1;
        IND = round(IND/Bouts.fs*fs);
        IND(IND<1) = 1;
        IND(IND>noSAM) = noSAM;
        noK = size(IND,1);
        if noK==0
            hp(sta) = plot(NaN,NaN,'color',color); %for legend
        end
        for k = 1:noK
            ind = IND(k,1):IND(k,2);
            dat = data(ind);
            hp(sta) = plot(t(ind),data(ind),'color',color);
        end
     
    end
    % legend(hp,Stages(:,1)')
end

%SETTINGS
linkaxes(ha,'x')
set(ha,'xlim',[0,tMax])
