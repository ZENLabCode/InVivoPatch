% plot_representation
%-------------------------------------------------------------------------
% Plots axes:
%   - Hypnogram
%   - VM data (downsampled)
%   - EMG
% Highlights selected stage bouts (they have margins)
clc; clear; close all
addpath(genpath('x:\OptoLab_v4.1\function'))


%PARAMETERS
%----------
%FILES VM (to select, other data files read from same folder)
% - auto searches for files using dir command
% - opens list to select from files found
filesVM = 'X:\1 TIDIS Lab\Tiago\Ephys_PROCESSED\In_vivo\Awake_Sleep\dPFC_Benzodiazepines_Processed\20230428_TC_P113_Day1\Evoked_Spikes\**\*_IN0.mat';

%OPTIONS
%stages (case sensitive, fieldname in *_vmStage.mat)
stages = {'NREM','REM','Wake'}; %out of {'NREM','REM','Wake'}
% colors = num2cell(lines(numel(stages)),2); %as cell
colors = {[ 0.9290    0.6940    0.1250],[0.8500    0.3250    0.0980],[ 0    0.4470    0.7410]};

%VM sample rate (fs)
% - original is 100000 Hz, which is huge (interpolates data)
%   PS: EMG has 1000, bout index 5000, hypnogram will be upsampled to EMG
% - fsVM must be <= original
fsVM = 5000; %[Hz]
%Filter function
% Dependent on data, freq-band and sampling rate, @(data,fBand,fs)
% E.g. - fBandEMG = [10,50]; 
%        filtFun = @(data,fs)passband_fourier(data,fBandEMG,fs)
%        Passband filter 10-50 Hz
fBandEMG = [10 250]; %[Hz] (no filter if empty) & [10 80]
fBandEEG = [1 40]; %[Hz] (no filter if empty) [.9 40]
fBandVM  = [0,fsVM/3]; %upper should be < fsVM/2 (only if fsVM < 100000)
filtFun  = @(data,fBand,fs)passband_fourier(data,fBand,fs);

%PLOT PROPERTIES
props.fig = {};
%axes (uses default axes size, might be changed by figure position)
props.axi = {};
props.fig_createAxes = {[70,50,50],[70,70,70],'pixel'}; %dx, dy & unit
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

%FILES
tmp = dir(filesVM);
if isempty(tmp)
    fprintf(2,'No filesVM found for: %s\n',filesVM)
    return
end
fnames = fullfile({tmp.folder},{tmp.name});
%remove files that does not yet have newly exported info field
warning off
ind = true(size(fnames));
for k = 1:numel(fnames)
    tmp = load(fnames{k},'info');
    ind(k) = isfield(tmp,'info');
end
fnames(~ind) = [];
if isempty(fnames)
    fprintf(2,'All VM files have no bout indices saved yet\n')
    fprintf(2,['Re-run ''extract_VM_hypnogram.m'' ',...
        '(AND ''spikeRemoval.m'')\n'])
    return
end
warning on
%select files
fnames = selectionList(fnames);
noPAT = numel(fnames); %vm-files
if noPAT==0
    fprintf(2,'No file selected\n')
    return
end

%porps must not be empty
tmp = fieldnames(props);
for k = 1:numel(tmp)
    if isempty(props.(tmp{k}))
        props.(tmp{k}) = {'visible','on'};
    end
end
%number of ...
noSTA = numel(stages);

%% PATH LOOP
nnPAT = numel(num2str(noPAT));
indent = blanks(2*nnPAT+2);
for pat = 1:noPAT
    [rPath,rFile] = fileparts(fnames{pat});
    id = regexp(rFile,'_','split'); id  = id{1};
    fprintf('%0*i/%i: %s\n',nnPAT,pat,noPAT,rPath)
    
    %FILES (specific order!)
    files = {...
        sprintf('%s.mat',rFile);...
        sprintf('%s_vmStage.mat',id);...
        'Hypnogram.mat';...
        'EMG.mat';....
        'EEG2.mat';...
        };
    %check
    ind = find(cellfun(@(x)exist(x,'file')~=2,fullfile(rPath,files)));
    if ~isempty(ind)
        fprintf(2,'%s Missing Files:\n',indent)
        for k = 1:numel(ind)
            fprintf(2,'%s   %s\n',indent,tmp{ind(k)})
        end
        continue
    end
    
    %% LOAD DATA
    fprintf('%s Load Data\n',indent);
    clear D
    for k = 1:numel(files)
        file = files{k};
        fprintf('%s  - %s\n',indent,file);
        D.(sprintf('f%i',k)) = load(fullfile(rPath,file));
    end
    
    %% PREPARE DATA (car of file order)
    fprintf('%s Prepare Data\n',indent);
    %VM
    VM = D.f1.resampled_data_mV(:);
    if D.f1.SampRate<fsVM
        error('Not for downsampling!!!')
    elseif D.f1.SampRate>fsVM
        %check
        if fsVM/2<fBandVM(2)
            error('fsVM/2 is smaller than fBandVM(2) (%g < %g)',...
                fsVM/2,fBandVM(2))
        end
        %filter
        fprintf('%s  - VM passband filter, %g - %g Hz\n',indent,...
            min(fBandVM(:)),max(fBandVM(:)));
        VM = filtFun(VM,fBandVM,fsVM);
        %interpolation
        fprintf('%s  - VM interpolation, %i --> %i Hz\n',...
            indent,D.f1.SampRate,fsVM);
        n1 = numel(VM); n2 = floor(n1/D.f1.SampRate*fsVM);
        VM = interp1((1:n1)/D.f1.SampRate,VM,(1:n2)/fsVM);
    end
    %BOUTS
    fsBOU = D.f2.fs;
    BOU.indEEG = D.f2.info.indEEG;
    BOU.indBOU = D.f2.info.indBouts;
    %EMG
    fsEMG = D.f4.SampRate;
    if isempty(fBandEMG)
        EMG = D.f4.resampled_data_mV(:);
    else
        fprintf('%s  - EMG passband filter, %g - %g Hz\n',indent,...
            min(fBandEMG(:)),max(fBandEMG(:)));
        EMG = filtFun(D.f4.resampled_data_mV(:),fBandEMG,fsEMG);
    end
    %EEG
    fsEEG = D.f4.SampRate;
    if isempty(fBandEEG)
        EEG = D.f5.resampled_data_mV(:);
    else
        fprintf('%s  - EEG passband filter, %g - %g Hz\n',indent,...
            min(fBandEEG(:)),max(fBandEEG(:)));
        EEG = filtFun(D.f4.resampled_data_mV(:),fBandEEG,fsEEG);
    end    
    %Hypnogram
    Hypnogram = D.f3.Hypnogram(:);
    fac = floor(numel(EMG)/numel(Hypnogram));
    if fac<1
        error('Downsampling Hypnogram is not implemented')
    elseif fac>0
        fprintf('%s  - Hypnogram upsampling by factor %i\n',indent,fac);
        Hypnogram = repmat(Hypnogram',fac,1);
        Hypnogram = Hypnogram(:);
    end
    Hypnogram(numel(Hypnogram)+1:numel(EMG)) = NaN;
    
    %% FIGURE/AXES
    close all; clc
    
    [dx,dy] = props.fig_createAxes{1:2};
    N = 4; %axes
    lenY = N*150+[1,N-1,1]*dy(:);
    hf = figure('position',[300 400 1000 lenY],props.fig{:});
    movegui(hf,'center'); drawnow
    ha = fig_createAxes(hf,[N,1],props.fig_createAxes{:});
    set(ha,'unit','normalized'); drawnow
    
    %PLOT Hypnogram & EMG
    noSAM = numel(EMG);
    ind = round(BOU.indEEG/fsBOU*fsEMG); %hypnogram & emg    
    ind = max([1,ind(1)]):min([noSAM,ind(2)]); %due to rounding
    t = ((1:noSAM)-ind(1)+1)/fsEMG; %[s]
    tMax = t(ind(2));
    %Hypnogram
    set(hf,'CurrentAxes',ha(1));
    plot(t(ind),Hypnogram(ind),props.plotHYP{:}); hold on;
    set(gca,'ylim',[0.7,3.2],'ytick',1:3,...
        'yticklabel',{'Wake','NREM','REM'})
    title({id,'Hypnogram'})
    %EMG
    set(hf,'CurrentAxes',ha(4));
    plot(t(ind),EMG(ind),props.plotEMG{:}); hold on;
    if isempty(fBandEMG)
        title('EMG')
    else
        title(sprintf('EMG (passband %g - %g Hz)',...
            min(fBandEMG(:)),max(fBandEMG(:))))
    end
    xlabel('Time [s]')
    ylabel('[mV]')
    %EEG
    set(hf,'CurrentAxes',ha(3));
    plot(t(ind),EEG(ind),props.plotEEG{:}); hold on;
    if isempty(fBandEEG)
        title('EEG')
    else
        title(sprintf('EEG (passband %g - %g Hz)',...
            min(fBandEEG(:)),max(fBandEEG(:))))
    end
    xlabel('Time [s]')
    ylabel('[mV]')
    
    %Stage bouts
    hp1 = NaN(noSTA,1);
    hp2 = NaN(noSTA,1);
    hp3 = NaN(noSTA,1);
    for sta = 1:noSTA
        stage = stages{sta};
        color = colors{sta};
        IND = round(BOU.indBOU.(stage)/fsBOU*fsEMG);
        IND(IND<1) = 1; IND(IND>noSAM) = noSAM;
        noK = size(IND,1);
        for k = 1:noK
            ind = IND(k,1):IND(k,2);
            if plotSTA.HYP
                set(hf,'CurrentAxes',ha(1));
                hp1(sta) = plot(t(ind),Hypnogram(ind),'color',color);
            end
            if plotSTA.EMG
                set(hf,'CurrentAxes',ha(4));
                hp2(sta) = plot(t(ind),EMG(ind),'color',color);
            end
            if plotSTA.EEG
                set(hf,'CurrentAxes',ha(3));
                hp3(sta) = plot(t(ind),EEG(ind),'color',color);
            end
        end
        if noK==0
            set(hf,'CurrentAxes',ha(1));
            hp1(sta) = plot(NaN,NaN,'color',color);
            set(hf,'CurrentAxes',ha(4));
            hp2(sta) = plot(NaN,NaN,'color',color);
            set(hf,'CurrentAxes',ha(3));
            hp3(sta) = plot(NaN,NaN,'color',color);
        end
    end
    if plotSTA.HYP
        set(hf,'CurrentAxes',ha(1));
        legend(hp1,stages)
    end
    if plotSTA.EMG
        set(hf,'CurrentAxes',ha(4));
        legend(hp2,stages)
    end
    if plotSTA.EEG
        set(hf,'CurrentAxes',ha(3));
        legend(hp3,stages)
    end
    
    
    %PLOT VM
    set(hf,'CurrentAxes',ha(2));
    noSAM = numel(VM);
    t = (1:noSAM)/fsVM; %[s]
    tMax = max([t(end),tMax]);
    plot(t,VM,props.plotVM{:}); hold on
    title('VM')
    ylabel('[mV]')
    %Stage bouts
    if plotSTA.VM
        hp = NaN(noSTA,1);
        for sta = 1:noSTA
            stage = stages{sta};
            color = colors{sta};
            IND = BOU.indBOU.(stage)-BOU.indEEG(1)+1;
            IND = round(IND/fsBOU*fsVM);
            IND(IND<1) = 1; IND(IND>noSAM) = noSAM;
            for k = 1:size(IND,1)
                ind = IND(k,1):IND(k,2);
                 hp(sta) = plot(t(ind),VM(ind),'color',color);
            end
            if isempty(IND)
                hp(sta) = plot(NaN,NaN,'color',color);
            end
        end
        legend(hp,stages)
    end
    
    %SETTINGS
    linkaxes(ha,'x')
    set(ha,'xlim',[0,tMax])
    zoom xon
    
end %path loop
