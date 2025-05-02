% spikeRemoval_plot
%-------------------------------------------------------------------------
% Plots results from spikeRemoval. For to check whether everthiny worked
% well.
%
% Thomas Rusterholz, 23 Jun 2021
%------------------------------------------------------------------------
clc; clear; close all
addpath(genpath('X:\OptoLab_v4.1\function'))

%PARAMETERS
%----------
%FILES
% - auto searching files using dir command
%   Removes files where spikeRemoval was not yet run
% - opens selection list to select ONE of files found (a lot of figures)
files = 'X:\1 TIDIS Lab\Tiago\Ephys_PROCESSED\In_vivo\Awake_Sleep\PFC_Processed\*\*\*_vmStage.mat';

%STAGES (case sensitive! variable names in read files)
stages = {'Wake','NREM','REM'};
traces = []; %all if empty, to test some

% stages = {'NREM'};
% traces = 8;

%PLOT PROPERTIES
props.axis = {'ylim',[-80,40]};
%plot
props.plotThreshold = {'color',[.5 .5 .5],'linewidth',1.5,...
    'linestyle',':'};
props.plotVM = {'color','b','linewidth',1.5};
props.plotSpikes = {'color','r','linewidth',1.5,'marker','.',...
    'MarkerSize',8};
props.plotAMP  = {'color','r','linewidth',1.5,'linestyle',':'};
props.plotFWHM = {'color','r','linewidth',1.5,'linestyle',':'};
props.plotArtifacts = {'color','k','linewidth',1.5,'linestyle',':'};
%text
% spike will be numbered. However, multiple spikes may obtain one label
% when they are close in time.
% E.g. label '2-3' if isi of spikes 2 and 3 is smaller than limISI
limISI = 0.1;


%MAIN SCRIPT
%-----------
scriptName = mfilename;
fprintf('%s\n%s\n',scriptName,repmat('-',size(scriptName)))
%non empty props
tmp = fields(props);
for k = 1:numel(tmp)
    props.(tmp{k}) = ['visible','on',props.(tmp{k})];
end

%FIND/SELECT FILES
tmp = dir(files);
if isempty(tmp)
    fprintf(2,'No files found for: %s\n',files)
    return
end
Files = fullfile({tmp.folder},{tmp.name});
%remove files that have no spike removal
% fprintf('Searching files ...\n')
% ind = false(size(Files));
% for k = 1:numel(Files)
%     ind(k) = ~ismember('no_spike',who('-file',Files{k}));
% end
% Files(ind) = [];
%select one file only
opt.title = 'Select ONE file only';
fname = selectionList(Files,[],opt);
while numel(fname)>1
    fname = selectionList(Files,[],opt);
end
if numel(fname)==0
    fprintf(2,'No File Selected!\n')
    return
end
fname = fname{1};

%LOAD DATA
[rPath,rFile] = fileparts(fname);
DATA = load(fname);
fs = DATA.fs;
threshold = DATA.no_spike.parameters.threshold;

%% LOOP STAGES
% clc; close all
for sta = 1:numel(stages)
    stage = stages{sta};
    DataRAW = DATA.(stage);
    DataSPK = DATA.no_spike.(stage);
    if isempty(DataSPK)
        fprintf('[\bNo %s trace!]\b\n',stage)
        continue
    end
    
    %LOOP TRACES
    Traces = 1:numel(DataRAW);
    if ~isempty(traces)
        Traces(~ismember(Traces,traces)) = [];
    end
    noTRA = numel(Traces);
    for tra = Traces
        datR = DataRAW{tra};
        datS = DataSPK(tra);
        t = (1:numel(datR))/fs;
        noSPK = numel(datS.indPKS);
        
        %PLOT
        hf = figure;
        %threshold / data
        plot([0,t(end)],[1,1]*threshold,props.plotThreshold{:});
        hold on
        hp  = plot(t,datS.data,props.plotVM{:});
        legSTR = {'Vm spike free'};
        %spikes
        for spk = 1:noSPK %spikes
            %plot spikes
            ind = datS.indSTA(spk):datS.indEND(spk);
            hp(2) = plot(t(ind),datR(ind),props.plotSpikes{:});
            legSTR{2} = 'Spikes';
            %plot amp & fwhm
            ind = datS.indPKS(spk);
            y   = datR(ind)-[0,datS.amp(spk)];
            plot(t(ind)+[0,0],y,props.plotAMP{:});
            plot(datS.t_fwhm(spk,:),sum(y)/2*[1,1],props.plotFWHM{:});
        end
        %artifacts
        ind = isnan(datS.data);
        if any(ind)
            d = datR;
            d(~ind) = NaN;
            hp(end+1) = plot(t,d,props.plotArtifacts{:});
            legSTR{end+1} = 'Artifacts';
        end
        
        %text
        text(t(end),threshold,'Threshold','rotation',90,...
            'HorizontalAlignment','center','VerticalAlignment','top')        
        
        xlabel('Time [s]')
        ylabel('Vm [mV]')
        str = strrep(rFile,'_','\_');
        title({sprintf('%s, %s Trace %i/%i',str,stage,tra,noTRA),...
            sprintf('Spikes N = %i, Spike Rate %g Hz',noSPK,datS.rate)})
        legend(hp,legSTR)
        %settings
        set(gca,'xlim',[0,t(end)],props.axis{:})
        zoom xon; drawnow
        
        %ADD SPIKE NUKBERS
        ind1 = 1; ind2 = 0;
        while ind1<=noSPK
            ind2 = min([noSPK,...
                find(datS.isi(ind1:end)>limISI,1,'first')+ind1-1]);            
            ind = datS.indPKS([ind1,ind2]);
            str = strjoin(cellfun(@num2str,...
                num2cell(unique([ind1,ind2])),'uniformoutput',false),'-');
            text(sum(t(ind))/2,max(datR(ind)),str,...
                'HorizontalAlignment','center',...
                'VerticalAlignment','bottom')
            ind1 = ind2+1;
        end
    end %loop traces
end %loop stages

return
%% TEST SPECTRA (current data only)
[pxxR,fR] = pwelch(datR,8*fs,6*fs,2*fs,fs);
[pxxS,fS] = pwelch(datS.data,8*fs,6*fs,2*fs,fs);

[pxxR,fR] = pwelch(datR,[],0,[],fs);
[pxxS,fS] = pwelch(datS.data,[],0,[],fs);

%plot data
close all; figure;
subplot(311);
plot(t,datR,'r'); hold on
plot(t,datS.data,'b');
xlabel('Time [s]')
ylabel('Vm [mV]]')
legend('Vm Raw','Vm spike free');
xlim([0,t(end)])
%plot power
subplot(312);
semilogy(fR,pxxR,'r'); hold on
semilogy(fR,pxxS,'b')
legend('Vm Raw','Vm spike free');
xlabel('Frequency [Hz]')
ylabel('Power [mV^2]')
%plot diff power
subplot(313);
dif = pxxR-pxxS;
plot(fR,dif,'linewidth',2);
xlabel('Frequency [Hz]')
ylabel('Diff Power')
yl = [min(dif),max(dif)];
ylim(yl+[-1,1]*diff(yl)*0.1)
zoom xon