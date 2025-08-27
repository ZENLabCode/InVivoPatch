%%  sealTest_expFit.m
%  --------------------------------------------------------------------
%  • Quantify membrane resistance (R_m) for every Vm epoch
%  • Classify each epoch by sleep stage (Wake / NREM / REM)
%  • DATA.(cellID).Wake/NREM/REM  = 1×nBout vectors  (MΩ)
%  • STAT matrices (cells × bouts) per stage
%  • Clean grid plot with per-stage medians
%  • Saves DATA & STAT to  RSP_Spontaneous_Evoked.mat
%  --------------------------------------------------------------------
clc; clear; close all
addpath(genpath('C:\Users\U1306265\Documents\OptoLab_v4.1\function'))

%% USER SETTINGS --------------------------------------------------------
fsDAT    = 5000;   fsHYP = 1;        % target VM rate, hypnogram rate
t0 = 0.255;  t1 = 0.280;            % slice window (s)
Vstep_mV = -5;                      % command step (mV)
tailFrac = 0.10;                    % last 10 % → Isteady
thr      = 150;                     % drop R_m > 150 MΩ

rawStages  = {'Wake','NREM','REM'};
dispStages = {'qWake','NREMs','REMs'};
stageCol   = [0.000 0.447 0.741;
              0.929 0.694 0.125;
              0.850 0.325 0.098];

saveDir  = 'C:\Users\U1306265\Desktop\NewFigures\SealtTest_Spontaneous\dAP5';
saveName = 'dAP5_Spontaneous_Evoked.mat';
%% ---------------------------------------------------------------------

%% FILE SELECTION -------------------------------------------------------
pat = fullfile(saveDir,'**','*_IN0_epochs.mat');
fl  = dir(pat);  fl([fl.isdir]) = [];
if isempty(fl), error('No *_IN0_epochs.mat files found.'); end
Files = selectionList(fullfile({fl.folder},{fl.name}));
if isempty(Files), error('No file selected.'); end
clear fl
%% ---------------------------------------------------------------------

%% MAIN ANALYSIS → DATA -------------------------------------------------
DATA = struct();

for f = 1:numel(Files)
    fileVm          = Files{f};
    [rPath,rFile]   = fileparts(fileVm);
    cellID          = matlab.lang.makeValidName(strrep(rFile,'_IN0_epochs',''));

    S = load(fileVm);                                % epoch file struct

    % ---------- choose / create 5-kHz Vm matrix -----------------------
    if isfield(S,'resampled_data_mV') && S.SampRate == fsDAT
        VM = S.resampled_data_mV;                    % already 5 kHz
    else
        VM = resample(S.Vm_epochs_mV, fsDAT, S.SampRate);
    end
    [nSam,nEp] = size(VM);           epDur = nSam/fsDAT;

    % ---------- load & upsample hypnogram -----------------------------
    H = load(fullfile(rPath,'Hypnogram.mat'),'Hypnogram').Hypnogram(:);
    H = repmat(H', fsDAT/fsHYP, 1);  H = H(:);
    if numel(H) < nSam*nEp, H(end+1 : nSam*nEp) = NaN; end

    % ---------- slice indices -----------------------------------------
    i0 = floor(t0*fsDAT) + 1;
    i1 = floor(t1*fsDAT);

    stageRm.Wake = [];  stageRm.NREM = [];  stageRm.REM = [];

    % ---------- epoch loop --------------------------------------------
    for ep = 1:nEp
        tr = VM(:,ep);  if numel(tr)<i1, continue, end

        seg   = tr(i0:i1);
        nTail = max(5, round(tailFrac*numel(seg)));
        I_ss  = mean(seg(end-nTail+1:end));          % pA
        if isnan(I_ss) || I_ss==0, continue, end

        Rm = (Vstep_mV*1e-3) / (I_ss*1e-12) / 1e6;   % MΩ
        if Rm > thr, continue, end

        % dominant hypnogram code for this epoch
        hIdx   = ((ep-1)*epDur*fsDAT+1) : (ep*epDur*fsDAT);
        codes  = H(hIdx);  codes = codes(~isnan(codes));
        if isempty(codes), continue, end
        stCode = mode(codes);

        switch stCode
            case 1, stageRm.Wake(end+1) = Rm;
            case 2, stageRm.NREM(end+1) = Rm;
            case 3, stageRm.REM(end+1)  = Rm;
        end
    end

    DATA.(cellID).Wake = stageRm.Wake;
    DATA.(cellID).NREM = stageRm.NREM;
    DATA.(cellID).REM  = stageRm.REM;

    fprintf('%s | Wake:%2d  NREM:%2d  REM:%2d\n',...
        cellID, numel(stageRm.Wake), numel(stageRm.NREM), numel(stageRm.REM));
end

%% BUILD STAT (cells × bouts) ------------------------------------------
cellNames = fieldnames(DATA);
STAT      = struct();
stageMaxB = zeros(1,3);

for k = 1:3
    st = rawStages{k};
    stageMaxB(k) = max(cellfun(@(c) isfield(DATA.(c),st)*numel(DATA.(c).(st)), cellNames));
    M = NaN(numel(cellNames), stageMaxB(k));
    for c = 1:numel(cellNames)
        if isfield(DATA.(cellNames{c}),st)
            v = DATA.(cellNames{c}).(st);
            M(c,1:numel(v)) = v;
        end
    end
    STAT.(st) = M;
end
STAT.cellNames = cellNames;

%% GRID PLOT (fixed axes + median lines) -------------------------------
figure('Color','w','Name','Rm per bout (median labelled)');

for c = 1:numel(cellNames)
    vals = []; for k=1:3, vals=[vals, STAT.(rawStages{k})(c,:)]; end
    vals = vals(~isnan(vals)); if isempty(vals), continue, end
    yLim   = [min(vals) max(vals)]; if diff(yLim)==0, yLim=yLim+[-1 1]; end
    yTicks = linspace(yLim(1), yLim(2), 3);

    for k = 1:3
        ax = subplot(numel(cellNames),3,(c-1)*3+k); hold(ax,'on');

        y = STAT.(rawStages{k})(c,:);  g = ~isnan(y);
        plot(ax, find(g), y(g), '-', 'Color',stageCol(k,:), 'LineWidth',1.4);

        if any(g)
            med = median(y(g));
            lightC = stageCol(k,:) + 0.5*(1-stageCol(k,:));
            plot(ax, [0 stageMaxB(k)], [med med],'--','Color',lightC,'LineWidth',2.4);
            text(stageMaxB(k)*0.02, med+0.04*diff(yLim), sprintf('%.0f',med),...
                 'Color',lightC,'FontWeight','bold','Parent',ax);
        end

        xlim(ax,[0 stageMaxB(k)]);   xticks(ax,0:5:stageMaxB(k));
        ylim(ax,yLim);               yticks(ax,yTicks);
        box(ax,'off');

        if k==1
            set(ax,'YTickLabel',num2str(yTicks','%.0f'));
        else
            set(ax,'YTickLabel',[]);
        end
        if c<numel(cellNames), set(ax,'XTickLabel',[]); end
        if c==1, title(ax,dispStages{k}); end
    end
end
sup = axes(gcf,'visible','off');
sup.YLabel.Visible='on';   sup.XLabel.Visible='on';
ylabel(sup,'Membrane resistance  R_m  (M\Omega)');
xlabel(sup,'Bout #');

%% SAVE DATA & STAT -----------------------------------------------------
if ~exist(saveDir,'dir'), mkdir(saveDir); end
save(fullfile(saveDir,saveName),'DATA','STAT','-v7');
fprintf('\nSaved DATA and STAT to %s\n', fullfile(saveDir,saveName));
