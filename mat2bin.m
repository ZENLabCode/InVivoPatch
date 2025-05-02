% mat2bin.m
%-------------------------------------------------------------------------
% Reads mat files and save them as bin-file.
% Made for EEG/EMG data

clc; clear; close all; fclose all;
addpath(genpath('X:\OptoLab_v4.1\function')) %for selectionList.m

%PARAMETERS
%----------
%FILES
% - no extension needed (will be corrected in script anyway)
%    - EEG (or EMG) files must be mat-files
%    - Save extension will be *.bin by default
Files.EEG = {'EEG1','EEG2','EMG'}; %cell, saves channels in this order
Files.bin = 'EEG'; %savename (default *.bin)

%PATHS (path list)
tmp = dir(fullfile('X:\1 TIDIS Lab\Tiago\Ephys_PROCESSED\In_vivo\Awake_Sleep\RC_Processed\*\Evoked_Spikes\*\EEG1.mat'));
rPaths = selectionList(unique({tmp.folder})); %unique paths!

%BIN-FILE DATA INFO
% Settable
%  - precision  : Save precision
%                 So far only made for 'int16' or 'uint16'
%                 (for scoring with SlipAnalysis)
%  - aqRange    : Aquisition Range (range from min to max!!!)
%                 In our lab with EEG data in [mV] we mostly use:
%                   aqRange = 10, for possible data range -5 to 5 mV
%                   (Crops data beyond range, however, EEG data beyond
%                    this range is anyway artifact)
%  - fac        : Factor to loaded data. Mainly for to change data unit.
%                 E.g. loaded data is in microvolts but you want it in
%                      millivolts, then set fac = 1/1000
% Not settable (uses default or available values)
%  - offset     : 0 for int16 or aqRange/2 for uint16
%                 (for EEG, values are above AND below zero)
%  - sampleRate : Uses sampling rate from read data files
%                 (maybe add option for re-sampling somewhen)
% Transformation functions, for physical (EEG) and digital (saved) data
%  - physical to digital data: (x + offset)/aqRange*2^16;
%  - digital to physical data: x/2^16*aqRange - offset
%                              Back-transformation when reading bin-files
opt.precision = 'int16'; %'int16' or 'uint16'
opt.aqRange   = 10; %aquisition range
opt.fac       = 1; %1 for no change


%CHECK OPTIONS (rough check for suspicious data)
% Has no effect to data export!
% Only prints some additional infos to command window (in orange)
% Print infos for the following conditions:
%  - pct : if percentage of data points beyond aquisition range is > pct
%          --> Maybe aqRange is set too small OR maybe wrong unit
%          PS: aquisition range is [-aqRange aqRange]/2
%          E.g.:
%            - pct = 0, print out if any data point is beyond range
%            - pct = 1, print out if 1% of data points are beyond range
%                       (for to limit print out, acceptable artifacts)
%  - lowA: if percentage of data beyond low amplitude range is < pct
%          --> Maybe wrong unit, e.g. may happen if aquision range is given
%              in millivolts but loaded data is in volts.
%          Then, change unit with opt.fac OR aqRange is set too large
%          anyway (reduce for better data resolution or live with it)
chk.pct  = 0.01; %Guess 0.01 is fine for possible artifacts
chk.lowA = opt.aqRange/1000; %.../1000 is probably fine for checking units

%MAIN SCRIPT
%-----------
scriptName = mfilename;
fprintf('%s\n%s\n',scriptName,repmat('-',size(scriptName)))

%INIT
if ~iscell(Files.EEG)
    Files.EEG = {Files.EEG}; %just in case
end
%number of ...
noPAT = numel(rPaths);
noFIL = numel(Files.EEG);
%default extensions
for k = 1:noFIL %remove, default *.mat
    [~,tmp] = fileparts(Files.EEG{k});
    Files.EEG{k} = [tmp,'.mat'];
end
[~,tmp] = fileparts(Files.bin);
Files.bin = [tmp,'.bin']; %default extension
%offset
switch opt.precision
    case 'int16'
        opt.offset = 0;
    case 'uint16'
        opt.offset = opt.aqRange/2;
    otherwise
        error('Precision ''%s'' not implemented',opt.precision)
end

%PATH LOOP
nnPAT  = numel(num2str(noPAT));
indent = blanks(2*nnPAT+2);
for pat = 1:noPAT
     rPath = rPaths{pat};
     fprintf('%*i/%i: %s\n',nnPAT,pat,noPAT,rPath)     
     %check for missing data files
     ind = cellfun(@(x)exist(x,'file')~=2,fullfile(rPath,Files.EEG));
     if any(ind)
         fprintf(2,'%s Missing Files: %s\n',indent,...
             strjoin(Files.EEG(ind),', '))
         continue
     end
     
     %LOAD DATA
     fprintf('%s Load Data (N = %i)\n',indent,noFIL)
     for fil = 1:noFIL
         rFile = Files.EEG{fil};
         fprintf('%s   %s\n',indent,rFile)
         tmp = load(fullfile(rPath,rFile));
         if fil==1
             noSAM = numel(tmp.resampled_data_mV);
             fs    = tmp.SampRate;
             Data  = NaN(noFIL,noSAM);
         elseif ~isequal(fs,tmp.SampRate)
             error('uups')
         end
         Data(fil,:) = tmp.resampled_data_mV*opt.fac;
     end
     
     %CHECK
     limA = opt.aqRange/2;
     limM = chk.lowA/2;
     %data points beyond range
     ind = Data(:)<-limA | Data(:)>limA;
     pct = sum(ind)/numel(ind)*100;
     if pct>chk.pct
         fprintf(['%s [\bData found beyond aquisition range [%g %g]',...
             ']\b\n'],indent,-limA,limA)
         fprintf('%s   [\bData Range: [%g %g]]\b\n',...
             indent,min(Data(:)),max(Data(:)))
         fprintf(['%s   [\bData Points beyond range: N = %i (%g %%)',...
             ']\b\n'],indent,sum(ind),pct)
     end
     %unit check
     ind = Data(:)<-limM | Data(:)>limM;
     pct = sum(ind)/numel(ind)*100;
     if pct<chk.pct
         fprintf(['%s [\bData is very low for aquisition range [%g %g]',...
             ']\b\n'],indent,-limA,limA)
         fprintf('%s   [\bData Range: [%g %g]]\b\n',...
             indent,min(Data(:)),max(Data(:)))
     end
     %test plot
     if false
         figure; t = (1:noSAM)/fs; xl = [0,t(end)];
         hp = NaN(noFIL+2,1);
         hp(1:noFIL) = plot(t,Data'); hold on;
         hp(end-1) = plot(xl,+[limA limA],':r','linewidth',1.5);
         hp(end-1) = plot(xl,-[limA limA],':r','linewidth',1.5);
         hp(end)   = plot(xl,+[limM limM],':k','linewidth',1.5);
         hp(end)   = plot(xl,-[limM limM],':k','linewidth',1.5);
         legend(hp,[Files.EEG(:)','Aquisiton Range','Low Amp'])
         set(gca,'xlim',xl); zoom on
     end
     
     %DATA TRANSFORMATION
     DataBIN = (Data + opt.offset)/opt.aqRange * 2^16;
     DataBIN = cast(DataBIN,opt.precision);
     fprintf('%s Bin-file Data Info (all channels)\n',indent)
     fprintf('%s   Save precision     : %s\n',indent,opt.precision)
     fprintf('%s   Sampling Rate [Hz] : %g\n',indent,fs)
     fprintf('%s   Aquisition Range   : %g\n',indent,opt.aqRange)
     fprintf('%s   Offset             : %g\n',indent,opt.offset)
     
     %SAVE DATA
     fid = fopen(fullfile(rPath,Files.bin),'w');
     fwrite(fid,DataBIN,opt.precision);
     fclose(fid);
     fprintf('%s Saved: %s\n',indent,Files.bin)
end %path loop