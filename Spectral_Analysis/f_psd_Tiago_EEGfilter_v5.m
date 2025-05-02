function [psd, EEG] = f_psd_Tiago_EEGfilter(channels)

% Bands and stages to be analyzed

Bands = {... {frequency band, field label}
    [0.5 1.5] , 'slowWave';...
    [2 4]  , 'delta';...
    [5 10]  , 'theta';...
    [8 12] , 'alpha';...
    [10 16] , 'spindle';...
    [30 60], 'gamma30_60';...
    };

Stages = {... {number, name}
    1,'Wake';...
    2,'NREM';...
    3,'REM';...
    };

%number of ...
noSTA = size(Stages,1);
noBND = size(Bands,1);
noCHA = numel(channels);

%FILENAMES
files = cell(size(channels));
for cha = 1:noCHA
    [rPath,rFile] = fileparts(channels{cha});
    if isempty(rPath)
        rPath = pwd;
    end
    files{cha} = fullfile(rPath,[rFile,'.mat']);
end
fileHyp = fullfile(rPath,'Hypnogram.mat');

%check
if exist(fileHyp,'file')~=2
    warning('Hypnogram.mat not found (no PSD return)')
    psd = NaN; EEG = []; return
end
ind = cellfun(@(x)exist(x,'file')~=2,files);
if all(ind)
    warning('No EEG file exist (no PSD return)')
    psd = NaN; EEG = []; return
elseif any(ind)    
    fprintf('[\bChannels not exist (removed)]\b\n')
    fprintf(' [\b- %s]\b\n',channels{ind})
    channels(ind) = [];
    files(ind)     = [];
    noCHA = numel(channels);
end

%LOAD DATA
%EEG channels
for cha = 1:noCHA
    tmp = load(files{cha});
    tmp.SampRate = double(tmp.SampRate);
    if cha==1
        fs    = tmp.SampRate;
        noSAM = numel(tmp.resampled_data_mV);
        EEG   = NaN(noSAM,noCHA);
    elseif fs~=tmp.SampRate || noSAM~=numel(tmp.resampled_data_mV)
        error('ggrrhhh')
    end
    EEG(:,cha) = tmp.resampled_data_mV(:);   
end

%Hypnogram
tmp = load(fileHyp);
% if  numel(tmp.Hypnogram)*fs-noSAM>fs || (noSAM-numel(tmp.Hypnogram)*fs)>60*fs
%     % Assumption: hypnogram epoch length = 1-s (respectively fs of 1 Hz)
%     % - scoring should not exceed eeg data (maybe last epoch not a full one)
%     % - if more than a minute not scored, not sure if epoch lengt of 1-s
%     %   is correct (one minute is quite conservative)
%     % --> re-check and comment error if you are sure it's correct
%     error('Check upsampling')
% end

Hypnogram = repmat(tmp.Hypnogram(:)',fs,1); %up-sampling
Hypnogram = Hypnogram(:);
Hypnogram(end+1:noSAM,1) = NaN;

%FILTER EEG
EEG0 = EEG; %keep for test plot
if true %maybe to switch off, at least for comparison
    IND = false(noSAM,noCHA);
    for cha = 1:noCHA
        eeg = EEG(:,cha);
        eegA = abs(eeg);
        difA = [0;diff(eegA)];
        eegM = movmedian(eegA,fs*10)+3*nanstd(eegA);
        %find artifacts
        ind = find(eegA>eegM);
        ind([false;diff(ind)<=2*fs]) = [];
        for k = 1:numel(ind)
            ind0 = ind(k);
            ind0 = min([ind0,find(difA(1:ind0)<0,1,'last')]);
            ind0 = ind0+(0:2*fs-1);
            eegA(ind0) = NaN;
            IND(ind0,cha) = true;
        end
    end
    
    %interppolation & filter
    indFilt = any(IND,2);
    x = 1:noSAM;
    warning('off','signal:internal:filteringfcns:ForcedHighpassDesign')
    for cha = 1:noCHA
        eeg = EEG(:,cha);
        eeg(indFilt) = interp1(x(~indFilt),eeg(~indFilt),x(indFilt),...
            'linear',nanmean(eeg));
        eeg(isnan(eeg)) = 0;
        %new filter
        %eeg = f_eegfilt(eeg(:)', fs, 0.5, fs/2);
        eeg = bandpass(eeg,[0.5, fs/2],fs);
        % eeg(ind) =  passband_fourier(eeg(ind),[0.5,fs/2],fs);
        %eeg = passband_fourier(eeg,[0.5,fs/2],fs);
        EEG(:,cha) = eeg(:);
    end
    EEG(indFilt,:) = NaN;
    warning('on','signal:internal:filteringfcns:ForcedHighpassDesign')
end

%TEST PLOT
if false
    close all; clc
    hf = figure; hf.WindowState = 'maximized'; drawnow
    ha = NaN(noCHA,1);
    for cha = 1:noCHA
        ha(cha) = subplot(noCHA,1,cha);
        plot(EEG0(:,cha)); hold on
        plot(EEG(:,cha));
        title(channels{cha})
    end
    linkaxes(ha,'xy')
    zoom on
    return
end

%REMOVE BANDS BEYOND FS/2
ind = true(noBND,1);
for bnd = 1:noBND
    ind(bnd) = max(Bands{bnd,1})>fs/2;
end
Bands(ind,:) = [];
noBND = size(Bands,1);

%CHANNEL LOOP
for cha = 1:noCHA
    channel = channels{cha};
    eeg = EEG(:,cha);
    
    %STAGES LOOP
    for sta = 1:noSTA
        stage = Stages{sta,1};
           
        %PWELCH
        tmp = eeg(Hypnogram==stage);
        tmp(isnan(tmp))=[];
        
        if ~isempty (tmp)
            %Define analyis windows (wd) in seconds*fs
            wd = 2*fs;
            wd_fft = hamming(wd);
            wd_number = round (length(tmp)/wd);
            startIndex = 1;
            window = 1;
            pxx_window = cell(1,  wd_number);
            
            while startIndex + wd -1 <= length (tmp)
                EEG_window = tmp (startIndex : startIndex + wd -1);
                
                if ~any (isnan(EEG_window))
                    try
                        [pxx,f] = pwelch(EEG_window,wd_fft,[],[],fs);
                        f(1)   = []; %remove offset (not wanted)
                        pxx(1) = [];
                        startIndex = startIndex + wd;
                        window = window+1;
                    catch
                        f   = 0.5:0.5:fs/2;
                        pxx = NaN(size(f));
                    end %try-catch loop
                    
                else
                    startIndex = startIndex + wd;
                end %NonNaN analysis
                
                %REMOVE SOCKET NOISE
                ind0  = -1:1; %replacments index
                noise = 50:50:fs/2; %50Hz & harmonics, for europe
                noise(noise+max(ind0)>max(f)) = [];
                
                for k = 1:numel(noise)
                    [~,ind] = min(abs(f-noise(k)));
                    ind1 = ind+ind0;         %index to replace
                    ind2 = ind1+numel(ind0); %signal after
                    
                    if ind2(end)>numel(f)
                        ind2 = ind1-numel(ind0); %signal before
                    end
                    pxx(ind1) = pxx(ind2);
                end %noise loop
                
                pxx_window{window} = pxx;
                
            end %while loop
            
            %APPEND RESULTS
            psd(cha).chanName = channel;
            psd(cha).EEG_fs = fs;
            psd(cha).indFilt  = indFilt;
            psd(cha).freq = f;
            pxx_window = cell2mat (pxx_window);
            pxx = nanmean (pxx_window, 2);
            psd(cha).psd(:,sta) = pxx;
            psd(cha).stages = Stages(:,2)';
            %band power
            for bnd = 1:noBND
                [band, lab] = Bands{bnd,:};
                if all(isnan(pxx))
                    psd(cha).(lab)(:,sta) = NaN;
                else
                    psd(cha).(lab)(:,sta) = bandpower(pxx,f,band,'psd');
                end
            end %band power loop
        else
            continue
        end %condition "contains stage"
    end %stage loop
end %channel loop