% convert_OpenBCI
%-------------------------------------------------------------------------
% Converts data recorded by (txt-file)
%
% Exports data to a bin-file (detrended and casted) and also saves a
% binInfo-file with informations of how to read the bin file (load it with
% load(*.binInfo,'-mat'). Channel data will also be saved as (mat-file),
% re-sampled to 1000 Hz ().
%
% Processing of multiple files is possible. It's also possible to get a
% filelist. Use functions export_fileList and import_fileList for that.
%
%
% Thomas Rusterholz, June 2018
%-------------------------------------------------------------------------

clc; clear; close all; fclose all;
% funPaths = genpath(pwd);
% addpath(funPaths)
addpath(genpath('X:\OptoLab_v4.1\function'))
%remove octave function paths
%tmp = regexp(path,';','split');
%ind = contains(tmp,'octavefunc');
%rmpath(strjoin(tmp(ind),';'));


%PARAMETERS
%----------
%FILES; filename, path (opening file browser) or import function
%       examples:
%        - files='OpenBCI-RAW-2018-06-01_15-04-47.txt';
%        - files='OpenBCI-RAW-*.txt';
%        - files={'OpenBCI-RAW....txt','OpenBCI-RAW....txt'};
%        - files=@() import_fileList(C:\filelists\filelist.txt'));
%        - files=@() import_fileList(C:\filelists\*'));
files='X:\1 TIDIS Lab\Tiago\Ephys_PROCESSED\In_vivo\Awake_Sleep\PFC_Processed\20211015_TC_P037_Day2\Evoked_Spikes'; 
%FORCE EXPORTING (if false, omits files that are already exported)
omitExisting=false; %true or false, true will not re-export

%READ INFO REC-FILE (equal for all selected files!)
%format for textscan; read only numeric channels, e.g no timestamp
readFormat='%*f %f %f %f %f %*f %*f %*f %*s';
%channel labels (for all read data by readFormat, same order!)
channels={'EEG1','EEG2','EMG','TTL'};

%OPTION
%bin-file options
opt.precision='int16'; %save precision ('int16' or 'uint16')
opt.fac=1/1000; %factor for saving data in mV (OpenBCI probably in microvolt)
opt.aquisitionRangeMax=10; %real range (from min to max) in mV
opt.chanSTIM={}; %stimulus channels, if any from channels (NOT WORKIN YET)
opt.dateFormat='yyyy-mm-dd_HH-MM-SS'; %date format within filename (empty if non)
%to export only info file without reading/export data
onlyINFO=true; %true of false

%MAIN PROGRAM
%-------------
%bin-file info, initialization
noCHA=numel(channels);
binFile.basedOn='';
binFile.headerLines=[];
binFile.binFile='';
binFile.recStart=[];
binFile.precision=opt.precision;
binFile.SampRate=[];
binFile.channels=channels;
binFile.chanSTIM=opt.chanSTIM;
binFile.acqRange=repmat(opt.aquisitionRangeMax,1,noCHA);
switch lower(binFile.precision)
    case 'int16'
        binFile.offset=zeros(1,noCHA);
    case 'uint16'
        binFile.offset=binFile.acqRange/2;
    otherwise
        error('uups')
end
binFile.bin2data=@(x,acqRange,offset) bsxfun(@plus,...
    bsxfun(@times,x,acqRange(:)/2^16),-offset(:));
binFile.data2bin=@(x,acqRange,offset) bsxfun(@times,...
    bsxfun(@plus,x,offset(:)), 2^16./acqRange(:));

%SELECT FILES
fnames=get_files(files);
if isempty(fnames)
    fprintf('No File Selected!\nProgram Canceled!\n')
    return
end
noFIL=numel(fnames); %files to export

%print out
fprintf('%s\n%s\n',mfilename,repmat('-',size(mfilename)))
%FILE LOOP
NO=numel(num2str(noFIL));
for fil=1:noFIL
    fname=fnames{fil};
    [rPath,id,rExt]=fileparts(fname);
    fprintf('%*i/%i: %s\n',NO,fil,noFIL,[id,rExt])
    if ~strcmpi(rExt,'.txt')
        fprintf('\b[\b - no txt-file]\b\n')
        continue
    end
    %save files
    sFiles=cellfun(@(x) fullfile(rPath,[x,'.mat']),channels,...
        'UniformOutput',false);
    ind=cellfun(@(x) exist(x,'file')==2,sFiles);
    if omitExisting && all(ind)
        fprintf('\b - omitted, already converted\n')
        continue
    end
    %bin-file info
    clear info
    binFile.basedOn=[id,'.txt'];
    binFile.binFile=[id,'.bin'];
    if ~isempty(opt.dateFormat)
        tmp=regexp(id,'\d*','match');
        ind=find(cellfun(@numel,tmp)==4,1,'last');
        ind=regexp(id,tmp{ind},'start');
        binFile.recStart=datevec(id(ind:end),opt.dateFormat);
    end
    
    %READ HEADER/SAMPLING RATE
    %initalization
    fs=[];
    noHLN=0; %number of header lines
    %open file
    fid=fopen(fname,'r');
    str=strtrim(fgetl(fid));
    while str(1)=='%'
        noHLN=noHLN+1;
        %get sampling rate from file
        tmp=regexpi(str,'sample rate');
        if ~isempty(tmp)
            ind=regexp(str,'[0-9.]');
            fs=str2double(str(ind));
        end
        %read next line
        str=strtrim(fgetl(fid));
    end
    if isempty(fs)
        error('uups')
    end
    binFile.SampRate=fs;
    binFile.headerLines=noHLN;
    
    %DATA CHECK
    %number of data channels
    noREC=numel(regexp(str,',','split')); %in rec-file
    noEXP=sum(readFormat=='%');           %expected
    noREA=noEXP-sum(readFormat=='*');     %read
    %one channel label per read channel
    if noCHA~=noREA        
        fprintf(['[\b   OMITTED: channel labels and read channels ',...
            'must be equal (%i vs %i)]\b\n'],noCHA,noREA)
        continue
    end
    %check read format
    if noEXP~=noREC
         fprintf(['[\b   OMITTED: readFormat expects %i channels, ',...
            'but file has %i channels]\b\n'],noEXP,noREC)
        continue
    end
    
    %SAVE BIN_FILE INFO (here for omit saving data)
    %info-file (as mat-file)
    fprintf('   saving: binInfo-file\n');
    sname=fullfile(rPath,[id,'.binINFO']);
    save(sname,'-struct','binFile','-mat')
   
    %READ DATA
    fseek(fid,0,'bof');
    Data=textscan(fid,readFormat,'delimiter',',','HeaderLines',noHLN);
    fclose(fid);
    minSAM = min(cellfun(@numel,Data));
    Data = cellfun(@(x)x(1:minSAM),Data,'uniformoutput',false);
    Data=detrend([Data{:}]*opt.fac)'; %data as matrix in mV
    
    %PREPARE DATA FOR BIN-FILE
    %%% test
    Data(1,1) =  binFile.acqRange(1)/2;
    Data(2,1) = -binFile.acqRange(1)/2;
    DataBIN = binFile.data2bin(Data,binFile.acqRange,binFile.offset);
    Data2   = binFile.bin2data(DataBIN,binFile.acqRange,binFile.offset);
    dif=Data-Data2;
    if ~all(dif==0)
        fprintf('Different, minmax: [%g %g]\n',min(dif(:)),max(dif(:)))
    else
        fprintf('Transformation ok\n')
    end
    if strcmpi(binFile.precision,'int16')
        disp([min(DataBIN(:))+2^15,max(DataBIN(:))-2^15])
    else
        disp([min(DataBIN(:)),max(DataBIN(:))-2^16])
    end
    
    ind=abs(Data)>binFile.acqRange(1)/2;
    if any(ind)
        str=sprintf(['file: %s.txt\n\n%i out of %i data points are ',...
            'outside AcuisitionRangeMax\nThey will be cropped for the ',...
            'bin-file (cast to %s),\nnot for the mat-files'],...
            id,sum(ind(:)),...
            numel(Data),binFile.precision);
        button=questdlg(str,'Warnings','continue','cancel this',...
            'cancel all','continue');
        switch lower(button)
            case 'continue'
            case 'cancel this'
                fprintf('\b[\b - canceled by user]\b\n')
                return
            case 'cancel all'
                fprintf('[\bCanceled ALL by user]\b\n')
                return
            otherwise
                error('uups')
        end
    end
    DataBIN=binFile.data2bin(Data,binFile.acqRange,binFile.offset);
    DataBIN=cast(DataBIN,binFile.precision);
    
    %SAVE DATA FILES
    %bin-file
    fprintf('\b, bin-file\n');
    sname=fullfile(rPath,[id,'.bin']);
    if exist(sname,'file')==2
        delete(sname)
    end
    fid=fopen(sname,'w');
    fwrite(fid,DataBIN,binFile.precision);
    fclose(fid); % close file
    %mat-file (one file per channel)
    fprintf('\b, mat-files\n');
    for cha=1:noCHA
        channel=channels{cha};
        clear tmp
        tmp.SampRate=fs;%fs;
        tmp.channel=channel;
        tmp.id=id;
        % tmp.resampled_data_mV=resample(Data(cha,:),...
        %    tmp.SampRate,binFile.SampRate);
        tmp.resampled_data_mV=Data(cha,:);
        sname=fullfile(rPath,sprintf('%s.mat',channel));
        save(sname,'-struct','tmp')
    end
    
end %file loop