% export_abfFiles
%-------------------------------------------------------------------------
% Reads abf-files (tested for version ABF2) and saves data to mat-files.
%
% Thomas Rusterholz, 11 Mai 2021
%------------------------------------------------------------------------

% file was copied from Z:\Tiago\Ephys_RAW\In_vivo\Awake_Sleep
clc; clear; close all; fclose all;
zenDrive = 'X:'; %used in script!!!
addpath(genpath(fullfile(zenDrive,'\OptoLab_v4.1\function'))) %for selectionList.m

%PARAMETERS
%-----------
%FILES
% - char, for to find files with dir command
%   e.g. 'C:\data\**\*.abf'
% - opens selection list for to select files from files found
%files = 'Z:\Tiago\Ephys_RAW\In_vivo\Awake_Sleep_OpenBCI\**\*.abf.lnk';
files = 'X:\1 TIDIS Lab\Tiago\Ephys_PROCESSED\In_vivo\Awake_Sleep\PFC_Processed\Pharmacology\NMDAR_Blockers\20210512_TC_P018_Day1\*.lnk';

%SAVENAME FUNCTION
% - saves mat-file (*.mat, other extensions will be replaced)
% - overwrites files!!!
% - based on read path & file
fun.savename = @(rPath,rFile)fullfile(rPath,[rFile,'.mat']);

%TEST PLOT
% - For huge data, slows down everything
plotData  = true; %plot data, logical
plotChans = []; %channel numbers to plot, plots all if empty


%MAIN SCRIPT
%-----------
scriptName = mfilename;
fprintf('%s\n%s\n',scriptName,repmat('-',size(scriptName)))

%FIND FILES
tmp = dir(files);
tmp([tmp.isdir]) = [];
if numel(tmp)==0
    fprintf(2,'No File Found For: %s\n',files)
    return
end
rnames = fullfile({tmp.folder},{tmp.name});
%select files
rnames = selectionList(rnames);
noFIL = numel(rnames);
if noFIL==0
    fprintf(2,'No File Selected!\n')
    return
end

%FILE LOOP
nnFIL  = numel(num2str(noFIL));
indent = blanks(2*nnFIL+2);
t0 = tic;
for fil = 1:noFIL
    tic
    rname = rnames{fil};
    [rPath,rFile,rExt] = fileparts(rname);
    fprintf('%*i/%i: %s\n',nnFIL,fil,noFIL,rname)
    if strcmpi(rExt,'.lnk')
        aserver = actxserver('WScript.Shell');
        tmp = aserver.CreateShortcut(rname);
        rname = tmp.TargetPath;
        rname(1:2) = zenDrive; %comming from a link
        [~,rFile,rExt] = fileparts(rname);
        aserver.release;
    end
    
    %LOAD DATA
    [Resampled_data_mV,~,hdr] = abfload(rname,'doDispInfo',false);
    C = permute(Resampled_data_mV,[1 3 2]);
    C = reshape(C,[],size(Resampled_data_mV,2),1);
    SampRate = 10^6/hdr.si;
    %start date/time vector
    recStart = datevec(num2str(hdr.uFileStartDate),'yyyymmdd');
    recStart(6) = hdr.uFileStartTimeMS/1000; % ~recTime(1) 
    recStart = datevec(datenum(recStart));
    %check
    if ~strcmpi(hdr.fFileSignature,'ABF2')
        warning(sprintf(['Script was only tested for abf-files of ',...
            'version ABF2.\nVersion %s needs to be tested too!'],...
            hdr.fFileSignature))
    end
    
%     %TEST PLOT
%     if plotData
%         fprintf('%s Plot Data\n',indent)
%         %channels to plot (index)
%         ind = plotChans(plotChans>0 & plotChans<=noCHA);
%         if isempty(ind)
%             ind = 1:noCHA;
%         end
%         %plot
%         hf = figure;
%         t = (1:size(Resampled_data_mV,1))/SampRate;
%         for k = ind
%             plot(t,Resampled_data_mV(:,k)); hold on
%             drawnow
%         end
%         %text
%         title({[rFile,rExt],datestr(recStart,'dd-mmm-yyyy HH:MM:SS')})
%         xlabel('Time [s]');
%         ylabel(sprintf('[%s]',hdr.recChUnits{ind(1)}))
%         legend(hdr.recChNames(ind))
%         %settings
%         set(gca,'xlim',([0,t(end)]))
%         drawnow
%         zoom xon
%     end
    
    %SAVE
    channels = hdr.recChNames;
    if ~iscell(channels)
        channels = {channels};
    end
    channels = cellfun(@(x)strrep(x,' ',''),channels,'uniformoutput',false);
    noCHA = numel(hdr.recChNames);
    for cha = 1:noCHA
        channel = channels{cha};
        [sPath,sFile] = fileparts(fun.savename(rPath,rFile));
        sname = fullfile(sPath,[sFile,'_',channel,'.mat']);
        info.info     = sprintf('Exported by %s.m, %s',scriptName,date);
        info.file     = rname;
        info.hdr      = hdr;
        info.recStart = recStart;
        resampled_data_mV = C(:,cha);
        save(sname,'info','resampled_data_mV','SampRate')
        fprintf('%s Saved: %s\n',indent,sname)
    end
    fprintf('%s ',indent); toc
end %file loop
toc(t0)
%just in case
fclose all;