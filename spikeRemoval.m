 
 % spikeRemoval
%-------------------------------------------------------------------------
% Removes spikes from Vm data. Detects spikes, estimates start & end points
% to replace with interpolated data. Additionally sets data that seems to
% be artifacts to NaN.
%
% Additionally calculates spike amplitude, fwhm, isi & rate.
% Data will be appended to read file (variable 'no_spike' so far).
% 
%
% Thomas Rusterholz, 23 Jun 2021
%------------------------------------------------------------------------
clc; clear; close all
addpath(genpath('X:\OptoLab_v4.1\function'))

%PARAMETERS
%----------
%FILES
% - auto searching files using dir command
% - opens selection list to select files from files found
files = 'X:\1 TIDIS Lab\Tiago\Ephys_PROCESSED\In_vivo\Awake_Sleep\PFC_Processed\Pharmacology\NMDAR_Blockers\*\*\*_vmStage_10sec.mat';
%files = 'Z:\Tiago\Ephys_PROCESSED\In_vivo\Awake_Sleep\dPFC_Processed\**\*_vmStage.mat';


%SPIKE DETECTION / REMOVAL
% - Detects spike by data crossing the threshold (>= threshold)
% - Replaces data around detected spikes with interpolated data.
%   Start & end indices are estimated, by selecting the closest index
%   before/after threshold crossing out of 3 conditions (if fullfilled)
%     1. Last points before/after threshold crossing with
%        positive/negative slope (Note: spike peak is upwards).
%     2. Using function fun(data, sampling rate).
%        First points before/after threshold crossing < 0.
%        PS: The function is meant to be of form data minus filtered data.
%            The filter should be set to get rid of spikes (low-pass)
%            Using moving median across 20 ms seems to be fine.
%     3. Margin to limit start/end indices distance to threshold crossing.
%        Start index within margin(1) before transition threshold
%        End   index within margin(2) after  transition threshold
par.threshold = -20; %[mV] (-10 seems to be fine for most)
par.fun = @(x,fs)x-movmedian(x,ceil(20/1000*fs/2)*2+1); %odd nbins
par.margin = [3,5]; %[ms]
%interpolation method (using interp1)
par.interpMethod = 'pchip'; %best was 'pchip', mostly fine was 'spline'

%STAGES (case sensitive! variable names in read files)
stages = {'Wake','NREM','REM'};

%ARTEFACTS
% minimal difference between peaks (artifacts if smaller)
art.minDiff = 0.002; %[s], set to zero if not wanted


%MAIN SCRIPT
%-----------
scriptName = mfilename;
fprintf('%s\n%s\n',scriptName,repmat('-',size(scriptName)))
%init info string
info = char({sprintf('Exported by %s, %s',scriptName,date);'';...
    'For each stage and trace (e.g. noSpike.NREM(1).indPKS)';...
    ' - indSTA : start index of spike replacment';...
    ' - indEND : end index of spike replacment';...
    ' - indPKS : index of spike peak';...
    ' - isi    : inter-spike interval [s] (peak to peak)';...
    ' - amp    : spike amplitude (spike peak minus replaced value)';...
    ' - fwhm   : full with at half maximum of spikes (using amp)';...
    ' - t_fwhm : interpolated timing of fwhm, [s]';...
    ' - data   : trace data with replaced spike data by interpolation';...
    ' - rate   : amount of spikes divided by trace duration, [Hz]';...
    ' - threshold: Vm at spike peaks minus minimum data ';...
    '';sprintf(...
    'Spike replacement by interpolation (interp1, method ''%s'')',...
    par.interpMethod);...
    });
indent = blanks(2);

%FIND/SELECT FILES
tmp = dir(files);
if isempty(tmp)
    fprintf(2,'No files found for: %s\n',files)
    return
end
Files = selectionList(fullfile({tmp.folder},{tmp.name}));
noFIL = numel(Files);
if noFIL==0
    fprintf(2,'No file selected\n')
    return
end

%BASIC STRUCTURE INIT
res0.indSTA = [];
res0.indEND = [];
res0.indPKS  = [];
res0.isi     = [];
res0.amp     = [];
res0.fwhm    = [];
res0.t_fwhm  = [0,2];
res0.rate    = 0; %zero for no spike
res0.threshold = [];
res0.data    = [];

%FILE LOOP
nnFIL = numel(num2str(noFIL));
indent = blanks(2*nnFIL+2);
for fil = 1:noFIL
    fname = Files{fil};
    [rPath,rFile] = fileparts(fname);
    fprintf('%*i/%i: %s\n',nnFIL,fil,noFIL,fname);
    
    %% LOAD DATA
    DATA = load(fname);
    fs = DATA.fs;
    margin = round(abs(par.margin)/1000*fs); %negative not allowed
    %init
    clear RES;
    RES.info = info;
    RES.parameters = par;
    
    %LOOP STAGES
    for sta = 1:numel(stages)
        stage = stages{sta};
        DataSTA = DATA.(stage);
        fprintf('%s  - %s\n',indent,stage);   
        
        %LOOP TRACES
        noTRA = numel(DataSTA);
        str = sprintf(', trace 0\n');
        fprintf('\b%s',str)
        for tra = 1:noTRA
            fprintf(repmat('\b',size(str)))
            str = sprintf(', trace %i of %i\n',tra,noTRA);
            fprintf('%s',str)
            data  = DataSTA{tra};
            data  = data(:);
            noSAM = numel(data);
            t     = (1:noSAM)'/fs;
            
            %% FIND SPIKES
            indSTA = find(...
                [data;-inf]>=par.threshold & ...
                [-inf;data]<par.threshold);
            indEND = find(...
                [data;-inf]<par.threshold  & ...
                [-inf;data]>=par.threshold)-1;
            %start/end correction
            noSPK = numel(indSTA);
            indPKS = NaN(noSPK,1);
            dataF = par.fun(data,fs);  %diff to filtered data
            for spk = 1:noSPK %start/end correction
                %spike peak
                ind = indSTA(spk):indEND(spk);
                [~,indP] = max(data(ind));
                indPKS(spk) = indP+ind(1)-1;
                
                %index left/right to thereshold crossing
                indLEF = indSTA(spk)-margin(1):indSTA(spk);
                indRIG = indEND(spk):indEND(spk)+margin(2);
                indLEF(indLEF<1)     = [];
                indRIG(indRIG>noSAM) = [];
                
                %STAR/END INDEX ESTIMATIONS (3 different tests)
                indS = NaN(1,3);
                indE = NaN(1,3);
                %min distance to threshold crossing
                indS(1) = max([indLEF(1)-1,1]);       %>=1
                indE(1) = min([indRIG(end)+1,noSAM]); %<=noSAM
                %positive slope at spike rise
                dat = data(indLEF);
                ind = find(dat(1:end-1)>dat(2:end),1,'last'); %neg slope
                if ~isempty(ind)
                    ind = ind+indLEF(1);
                    if ind<indSTA(spk) %before indSTA
                        indS(2) = ind;
                    end
                end
                %negative slope at spike decay
                dat = data(indRIG);
                ind = find(dat(1:end-1)<dat(2:end),1,'first'); %pos slope
                if ~isempty(ind)
                    ind = ind+indRIG(1)-1;
                    if ind>indEND(spk) %after indEND
                        indE(2) = ind;
                    end
                end
                %fun
                ind1 = find(dataF(indLEF)<0,1,'last')+indLEF(1)-1;
                ind2 = find(dataF(indRIG)<0,1,'first')+indRIG(1)-1;
                if ~isempty(ind1)
                    indS(3) = ind1;
                end
                if ~isempty(ind2)
                    indE(3) = ind2;
                end
                
%                 %%% test plot
%                 close all; hf = figure;
%                 ind = indLEF(1):indRIG(end);
%                 plot(t(ind),data(ind),'b','marker','.'); hold on;
%                 hline(par.threshold)
%                 plot(t(indSTA(spk)),data(indSTA(spk)),'rx','linewidth',2)
%                 plot(t(indEND(spk)),data(indEND(spk)),'rx','linewidth',2)
                
                %estimated start/end
                indSTA(spk) = max(indS);
                indEND(spk) = min(indE);
                
%                 %%%
%                 plot(t(indSTA(spk)),data(indSTA(spk)),'ro','linewidth',2)
%                 plot(t(indEND(spk)),data(indEND(spk)),'ro','linewidth',2)
%                 plot(t(indS(2)),data(indS(2)),'c>','linewidth',2)
%                 plot(t(indE(2)),data(indE(2)),'c<','linewidth',2)
%                 ind = false(size(data));
%                 ind(indSTA(spk):indEND(spk)) = true;
%                 dd = interp1(t(~ind),data(~ind),t(ind),par.interpMethod);
%                 plot(t(ind),dd,'g'); zoom xon
%                 return
            end
            
            %% ARTIFACTS I (no spikes)
            ind = diff(indPKS)<art.minDiff*fs;
            data2 = data;
            if any(ind)
                IND1 = unique([find(ind(:));find(ind(:))+1]);
                IND2 = IND1;
                while ~isempty(IND1)
                    ind = find(ismember(IND1,(0:numel(IND1)-1+IND1(1))));
                    data2(indSTA(ind(1)):indEND(ind(end))) = NaN;
                    IND1(ind) = [];
                end
                indSTA(IND2) = [];
                indEND(IND2) = [];
                indPKS(IND2) = [];
                clear IND1 IND2
                noSPK = numel(indSTA);
            end      
            
            %INIT
            res = res0;
            res.indSTA = indSTA;
            res.indEND = indEND;
            res.indPKS = indPKS;
            res.data   = data2;
            if noSPK==0
                RES.(stage)(tra) = res;
                continue
            end
            
            %SPIKE REMOVAL / ANALYSIS
            res.amp = NaN(noSPK,1);
            T1 = NaN(noSPK,1); T2 = NaN(noSPK,1);
            for spk = 1:noSPK
                %artifacts II (trailing ending spikes)
                ind1 = indSTA(spk):indPKS(spk);
                ind2 = indPKS(spk):indEND(spk);
                if (ind1(1)==1 && all(data(ind1)>par.threshold)) || ...
                        (ind2(end)==noSAM && all(data(ind2)>par.threshold))
                    res.data(indSTA(spk):indEND(spk)) = NaN;
                    continue
                end
                %spike removal
                indS = false(size(data));
                indS(indSTA(spk):indEND(spk)) = true;
                res.data(indS) = interp1(t(~indS),data(~indS),t(indS),...
                    par.interpMethod);
                %amplitude
                amp = data(indPKS(spk))-res.data(indPKS(spk));
                res.amp(spk) = amp;
                
                %time points fwhm
                dat = data(indS);
                tt  = t(indS);
                y   = res.data(indPKS(spk))+amp/2;
                [~,ind] = max(dat);
                ind1 = 1:ind(1);
                ind2 = ind(1):numel(dat);
                if any(dat(ind1)<y) && any(dat(ind2)<y)
                    T1(spk) = interp1(dat(ind1),tt(ind1),y,...
                        par.interpMethod);
                    T2(spk) = interp1(dat(ind2),tt(ind2),y,...
                        par.interpMethod);
                end
            end
            %remove invalid spikes
            ind = isnan(res.amp);
            if any(ind)
                res.amp(ind)    = [];
                res.indSTA(ind) = [];
                res.indEND(ind) = [];
                res.indPKS(ind) = [];
                T1(ind) = []; T2(ind) = [];
                noSPK = numel(res.indPKS);
            end
            %additional calcs
            res.isi  = diff(indPKS)/fs;
            res.rate = noSPK/t(end);
            res.fwhm = T2-T1;
            res.t_fwhm = [T1,T2];
            res.threshold = res.data(res.indPKS)-min(res.data);
            if ~isequal(fields(res),fields(res0))
                error('Check init res0')
            end            
            RES.(stage)(tra) = res;
            
            %PLOT (for testing)
            if false
                hf = figure;
                plot([0,t(end)],[1,1]*par.threshold,':','linewidth',1.5,...
                    'color',zeros(1,3)+0.5); hold on
                hp(1) = plot(t,res.data,'b','linewidth',1.5);
                for spk = 1:noSPK %spikes
                    ind = res.indSTA(spk):res.indEND(spk);
                    hp(2) = plot(t(ind),data(ind),'r','linewidth',1.5,...
                        'marker','.','MarkerSize',8);
                    %amp and fwhm
                    ind = res.indPKS(spk);
                    plot(t(ind)+[0,0],data(ind)-[0,res.amp(spk)],'r:',...
                        'linewidth',1.5,'marker','o','MarkerSize',6);
                    y = res.data(ind)+res.amp(spk)/2;
                    plot(res.t_fwhm(spk,:),[y,y],'r:','linewidth',1.5,...
                        'marker','o','MarkerSize',6);                 
                end
                %text
                text(t(end),par.threshold,'Threshold','rotation',90,...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment','top')
                xlabel('Time [s]')
                ylabel('Vm [mV]')
                title({strrep(rFile,'_','\_'),...
                    sprintf('%s, trace %i/%i',stage,tra,noTRA),...
                    sprintf('Spikes N = %i',noSPK)})
                legend(hp,'Spike free','Spikes')
                %settings
                set(gca,'xlim',[0,t(end)])
                zoom xon; drawnow
            end
        end %loop traces
        if ~isfield(RES,stage)
            RES.(stage) = [];
        end
    end %loop stages   
    
    %APPEND
    no_spike = RES; %other variable name
    save(fname,'no_spike','-append')   
    fprintf('%s Data appended!\n',indent)
end %loop files