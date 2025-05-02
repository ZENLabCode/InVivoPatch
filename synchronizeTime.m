% synchronizeTime
%-------------------------------------------------------------------------
% To synchronize time of abf-recordings with TTL. Plots data for to set
% desired time shift to all abf data with TTL.
%   Appends start date/time vector to TTL mat-file (from OpenBCI txt-files)
% and shifted start date/time vector to abf mat-files.
%
% PS: - needs TTL mat-file (exported from OpenBCI txt-files) and exported 
%       ABF data mat-files (also needs abf-files or links to find
%       corresponding filenames).
%     - ABF data will be resampled to same sample rate as TTL
%     - You can use the arrow keys to shift the data in steps of 1/SampRate
%       In x-direction, Alt+arrow for 20 steps and Ctrl+arrow for 100 steps
%       (see fun_keyboard.m)
%
%
% Thomas Rusterholz, 25 Mai 2021
%------------------------------------------------------------------------
clc; clear; close all; fclose all;
addpath(genpath('X:\OptoLab_v4.1\function')) %for selectionList.m

%PARAMETERS
%-----------
%FILES
% - char, for to find OpenBCI txt-files with dir command
%   e.g. files = 'C:\data\**\OpenBCI*.txt'
%   Searches TTL mat-file and ABF-files in same path
% - opens selection list for to select files from files found
files = ['X:\1 TIDIS Lab\Tiago\Ephys_PROCESSED\In_vivo\Awake_Sleep\PFC_Processed\20210623_TC_P024_Day1\OpenBCI-*.txt'];
%to find TTL and ABF data files (without path!!!)
fileTTL = 'TTL.mat'; %must be unique in corresponding path
filesABF = {... takes first that exist
    '*_IN1.mat';...
    '*_IN0.mat';...
    };
% fileABF = '*abf*'; %PS: needs abf-files but reads mat-files!

%OPTIONS
%time margin
% - plots data of abf-files +- margin. Should be larger than expected time
%   delay to TTL. In tested data, a margin of 50 seconds was fine.
margin = 500; %50; %in [s]
%variable name exported to ABF mat-files
%(uses 'startDate' for TTL mat-file by default)
% PS: use always the same name. Will be used in other scrips. So probably
%     do not change or only the first time you use it.
varName = 'startDateShifted';
%reExport
% if true , you can select files from all files found.
% if false, removes files that already contain the variable varName
%           It's for to list only files that were not yet synchronized.
reExport = true;

%test mode will save nothing!
testMode = false;

%MAIN SCRIPT
%-----------
scriptName = mfilename;
fprintf('%s\n%s\n',scriptName,repmat('-',size(scriptName)))
if testMode
    warning('Data will NOT be Saved (testing mode)')
end

%FIND FILES
%search OpenBCI files
tmp = dir(files);
tmp([tmp.isdir]) = [];
if numel(tmp)==0
    fprintf(2,'No File Found For: %s\n',files)
    return
end
files = tmp;
%search corresponding TTL and abf-files
noFIL = numel(files);
Files = cell(noFIL,3); %OpenBCI, TTL and ABF files 
for fil = 1:noFIL
    rPath = files(fil).folder;
    rFile = files(fil).name;
    fileB = fullfile(rPath,rFile);
    
    %find TTL file
    tmp = dir(fullfile(rPath,fileTTL));
    tmp([tmp.isdir]) = [];
    switch numel(tmp)
        case 0
            continue
        case 1
            fileT = fullfile(rPath,tmp.name);
        otherwise
            error('TTL file is not unique for: %s',fullfile(rPath,fileTTL))
    end
    if numel(tmp)==0
        continue %no such file
    end
    
    %find abf-files
    for k = 1:numel(filesABF)
        tmp = dir(fullfile(rPath,filesABF{k}));
        if numel(tmp)~=0
            break
        end
    end
    if numel(tmp)==0
        continue %no such file
    end
    filesA = {tmp.name};
    filesA = cellfun(@(x)[x(1:find(x=='.',1,'first')),'mat'],filesA,...
        'UniformOutput',false); %change extension
    filesA = fullfile(rPath,filesA);
    
    %append
    Files(fil,:) = {fileB,fileT,filesA};
end
clear fileB fileT filesA files %clean up
%only files having all needed files
Files(cellfun(@isempty,Files(:,1)),:) = [];
noFIL = size(Files,1);
if noFIL==0
    fprintf(2,'No Path Found with TTL & ABF-files\n')
    return
end
%remove files that were already performed
if ~reExport
    ind = false(noFIL,1);
    for fil = 1:noFIL
        files = Files{fil,3};
        ind = false(size(files));
        for k = 1:numel(files)
            tmp = matfile(files{k});
            if ismember(varName,fields(tmp))
                ind(fil) = true;
                break
            end
        end
    end
    Files(ind,:) = [];
    noFIL = size(Files,1);
    if noFIL==0
        fprintf(2,'All Files Are Already Synchronized\n')
        return
    end
end
%select files
opt.unique = false; opt.sort = false;
[~,ind] = selectionList(Files(:,1),[],opt);
Files = Files(ind,:);
noFIL = size(Files,1);
if noFIL==0
    fprintf(2,'No File Selected!\n')
    return
end

%FILE LOOP
nnFIL  = numel(num2str(noFIL));
indent = blanks(2*nnFIL+2);
for fil = 1:noFIL
    [fileBCI,fileTTL,filesABF] = Files{fil,:};
    [rPath,rFile,rExt] = fileparts(fileBCI);
    fprintf('%*i/%i: %s\n',nnFIL,fil,noFIL,rPath)
    
    %% GET START DATE/TIME TTL
    fprintf('%s Read %s\n',indent,[rFile,rExt]);
    fid = fopen(fileBCI,'r');
    str = strtrim(fgetl(fid));
    while str(1)=='%' || isempty(str)
        str = strtrim(fgetl(fid));
    end
    fclose(fid);
    str = regexp(str,' ','split');
    str = regexp(str{end},':','split');
    %create datevec
    [~,rFile] = fileparts(fileBCI);
    startDate = datevec(rFile,'OpenBCI-RAW-yyyy-mm-dd_HH-MM-SS');
    startDate(4) = str2double(str{1});
    startDate(5) = str2double(str{2});
    startDate(6) = str2double(str{3});
    
    %% LOAD DATA
    %data TTL
    [~,lab,ext] = fileparts(fileTTL);
    fprintf('%s Load %s\n',indent,[lab,ext]);
    tmp = load(fileTTL);
    ttl = tmp.resampled_data_mV(:);
    if all(isnan(ttl))
        fprintf(2,'%s All TTL is NaN\n',indent)
        continue
    end
    fs  = tmp.SampRate;
    noSAM = numel(ttl);
    t = (1:noSAM)'/fs; %in [s]
    if ~testMode
        save(fileTTL,'startDate','-append')
    end
    %print start/end data
    fprintf('%s   Rec Start: %s\n',indent,...
        datestr(startDate,'dd-mmm-yyyy HH:MM:SS.FFF'))
    fprintf('%s   Rec End  : %s\n',indent,...
        datestr(datenum(startDate+[0,0,0,0,0,t(end)]),...
        'dd-mmm-yyyy HH:MM:SS.FFF'))
    
    %data ABF
    noABF = numel(filesABF);
    VM    = cell(noABF,4); %t, data & label
    omit  = false(1,noABF);
    for k = 1:noABF
        fname = filesABF{k};
        [~,lab,ext] = fileparts(fname);
        fprintf('%s Load %s\n',indent,[lab,ext]);
        tmp = load(fname);
        %resample (that's probably enough)
        datV = resample(tmp.resampled_data_mV(:),fs,tmp.SampRate);
        tV = (1:numel(datV))'/fs + etime(tmp.info.recStart,startDate);
        %print start/end data
        fprintf('%s   Rec Start: %s\n',indent,...
            datestr(tmp.info.recStart,'dd-mmm-yyyy HH:MM:SS.FFF'))
        fprintf('%s   Rec End  : %s\n',indent,...
            datestr(datenum(tmp.info.recStart+[0,0,0,0,0,tV(end)]),...
            'dd-mmm-yyyy HH:MM:SS.FFF'))
        if tV(end)<t(1) || tV(1)>t(end)
            omit(k) = true;
            fprintf(2,'%s   Omit: start/end dates do not overlap with TTL\n',...
                indent)
            continue
        end
        if isfield(tmp,varName)
            dt = etime(tmp.(varName),startDate) - ...
                etime(tmp.info.recStart,startDate);
            VM(k,:) = {tV,datV,lab,dt};
        else
            VM(k,:) = {tV,datV,lab,NaN};
        end
    end
    %correct skipped files
    VM(omit,:)     = [];
    filesABF(omit) = [];
    noABF = numel(filesABF);
    if noABF==0
        fprintf(2,'%s All Omitted\n',indent)
        continue
    end
    
    %% FIGURE
    %clc; close all;
    fprintf('%s Plot:\n',indent)
    hf = figure('WindowState','maximized',...
        'KeyPressFcn',@fun_keyboard,...
        'CloseRequestFcn',@fun_close);
    %axes
    cols = ceil(sqrt(noABF));
    rows = ceil(noABF/cols);
    dx = [90,90,60]; dy = [60,80,60];
    ha1 = fig_createAxes(hf,[rows,cols],dx,dy,'pixel');
    ha2 = fig_createAxes(hf,[rows,cols],dx,dy,'pixel');
    set(ha1,'box','off','nextplot','add',...
        'xaxislocation','top','xticklabel',[],...
        'yaxislocation','right','unit','normalized')
    set(ha2,'box','off','nextplot','add','color','none',...
        'unit','normalized')
    set(ha1(noABF+1:end),'visible','off')
    set(ha1(noABF+1:end),'visible','off')
    
    % PLOT
    g.dt    = 1/fs;
    g.hp    = NaN(noABF,1);
    g.t     = cell(noABF,1);
    g.shift = zeros(noABF,1);
    for k = 1:noABF
        [tV,datV,lab,dt] = VM{k,:};
        
        
        %TTL CROPPED
        [~,ind1] = min(abs(t-tV(1)));
        [~,ind2] = min(abs(t-tV(end)));
        ind1 = max([ind1-margin*fs,1]);
        ind2 = min([ind2+margin*fs,noSAM]);
        datT = ttl(ind1:ind2);
        tT = t(ind1:ind2);
        
        %ESTIMAT TIME SHIFT
        if isnan(dt) %for new estimation
            %point 2
            ind0 = floor(fs/10); %start index to ignore inital signal
            indP = []; df = 1/100; f = 3+df;
            while isempty(indP)
                f = f-df;
                lim  = mean(datV) + f*std(datV);
                indP = find(datV(ind0+1:end)>lim,1,'first')+ind0;
            end
            dif  = abs(diff(datV(1:indP)));
            indPV = []; df = 1; f = 11;
            while isempty(indPV)
                f = f-df;
                lim = dif(end)/f;
                indPV = find(dif<lim,1,'last');
            end
            %point 1
            [tmp,ind0] = min(abs(tT-tV(indPV))); %start point
            ind  = max([1,ind0-50*fs]):ind0+1;
            indP = []; df = 1/100; f = 3+df;
            while isempty(indP) %positive peak
                f = f-df;
                lim  = mean(datT(ind)) + f*std(datT(ind));
                indP = find(datT(ind)>lim,1,'last')+ind(1)-1;
            end
            ind   = indP-3*fs:indP+1;
            indPT = []; df = 1/100; f = 3+df;
            while isempty(indPT)
                f = f-df;
                %                 lim   = mean(datT(ind)) - f*std(datT(ind));l
                lim   = mean(datT(ind)) - f*std(datT(ind));
                indPT = find(datT(ind)<lim,1,'first')+ind(1)-1;
            end
            lim = mean(datT(ind)) - std(datT(ind))/4;
            indPT = find(datT(1:indPT)>lim,1,'last')-1; drawnow
            %time shift
            dt =  tT(indPT)-tV(indPV);
            fprintf('%s   New estimated delay %s: %g s\n',indent,lab,dt)
        else
            fprintf('%s   Previous saved delay %s: %g s\n',indent,lab,dt)
        end
        
        %PLOT
        %ttl
        set(hf,'CurrentAxes',ha1(k))
        plot(tT,datT,'b');
        ylabel('TTL','color','b')
        %vm
        set(hf,'CurrentAxes',ha2(k))
        hp = plot(tV+dt,datV,'r');
        ylabel('Vm','color','r')
        %text
        title(strrep(lab,'_','\_'))
        xlabel({'Time [s]',sprintf('Shift Vm: %.3f',dt)})
        %settings
        linkaxes([ha1(k),ha2(k)],'x')
        xl = [min([tT(1),tV(1)+dt]),max([tT(end),tV(end)+dt])];
        set(gca,'xlim',xl+[-1,1]*0.01*diff(xl))
        
        %append
        g.hp(k)    = hp;
        g.t{k}     = tV;
        g.shift(k) = dt;
    end
    set(ha1(noABF+1:end),'visible','off')
    set(ha2(noABF+1:end),'visible','off')
    sgtitle({regexprep(rPath,{'\','_'},{'\\\\','\\_'}),...
        'Zoom to ROI and deselect zoom',...
        'Click into axis, use arrow keys to shift Vm data',...
        'Close figure when all done'})
    zoom xon
    guidata(hf,g)
    waitfor(hf)
    
    %APPEND DATA
    if ~testMode
        fprintf('%s Append ''%s'' to:\n',indent,varName)
        clear res
        for k = 1:noABF
            %use always 
            fname = filesABF{k};
            [sPath,sFile,sExt] = fileparts(fname);
            tmp = regexp(sFile,'_','split');
            sFile = strjoin([tmp(1),'IN0'],'_');
            fname = fullfile(sPath,[sFile,sExt]);
            %
            tmp = load(fname,'info');
            d = tmp.info.recStart + [0,0,0,0,0,g.shift(k)];
            res.(varName) = datevec(datenum(d));
            save(fname,'-struct','res','-append')
            %text
            [~,rFile,rExt] = fileparts(fname);
            fprintf('%s   %s\n',indent,[rFile,rExt])
        end
    else
        fprintf('%s [\bNothing Saved!]\b\n',indent)
    end
    
    %%
end %file loop
fclose all; %just in case