clc; clear; close all;
addpath('X:\1 TIDIS Lab\Tiago\MatLab_Scripts\In_Vivo_Patch\Spectral_Analysis')
addpath(genpath('X:\OptoLab_v4.1\function'))

%% PARAMETERS
%----------
%MAIN PATH
% Opens browser to select path, starting in this one
mPath = 'X:\1 TIDIS Lab\Tiago\Ephys_PROCESSED\In_vivo\Awake_Sleep\PFC_Processed\';

channels = {'EEG1','EEG2'};
%options for function f_psd_Tiago.m
opt.SampRateHyp = 1;

%MAIN SCRIPT
%-----------
scriptName = mfilename;
scriptPath = fileparts(which(scriptName));
addpath(scriptPath)

%GET DATA PATH
rPath = uigetdir(mPath,'Select Paths');
if ~ischar(rPath)
    fprintf(2,'No path selected!\n')
    return
end
cd(rPath)

[psd,EEGf] = f_psd_Tiago_EEGfilter_v5(channels);
%% Plot EEG non-filtered & filtered
if false
    noCHA = numel(channels);
    hf = figure('WindowState','maximized'); drawnow
    ha = NaN(noCHA,1);
    for cha = 1:numel(channels)
        channel = channels{cha};
        ha(cha) = subplot(noCHA,1,cha);
        tmp = load(channel);
        tmp.SampRate = double(tmp.SampRate); %convert to double to avoid errors with classes.
        if cha==1
            t   = (1:numel(tmp.resampled_data_mV))/tmp.SampRate;
        end
        plot(t,tmp.resampled_data_mV); hold on
        plot(t,EEGf(:,cha)); hold on
        %text
        legend({'original','filtered'},'location','northeastoutside')
        title(channel)
        xlabel('Time [s]')
    end
    tmp = cellfun(@(x)strrep(x,'_','\_'),strsplit(rPath,filesep),...
        'uniformoutput',false);
    sgtitle(strjoin(tmp(end-1:end),', '))
    %settings
    linkaxes(ha,'x')
    set(ha,'xlim',[0,t(end)])
    zoom on
end

EEG_fs = psd(1).EEG_fs;
%% SAVE
save('EEG_filt.mat','psd', 'EEGf', 'EEG_fs');