%VM/EEG downsampling for plot
A = DATA.x23705005.Wake.EEG_bout.b1;
n = numel(A);
VM_noSAM = floor(n/1000*200);
A_data = interp1((1:n)/1000,A,(1:VM_noSAM)/200);
close all;
plot (A_data);
%% Extract Basic analysis
zondition = 'cumsumY';
stages = {'Quiet_Wake','NREM'};
rows = 30;
Copy = NaN (rows,4);
Copy_normalized = NaN (rows,4);

for sta = 1: numel (stages)
    stage = stages{sta};
    if sta == 1
        stagez = 'Quiet_Wake';
        data_norm = DZP_STAT.(stagez).VM.(zondition);
    else
        stagez = stage;
        data_norm = DZP_STAT.(stagez).VM.(zondition);
    end
    
    data = ATI_STAT.(stage).VM.(zondition);
    data_norm = nanmean (data_norm);
    data_normalized = (data./data_norm)*100 - 100;
    A = numel (data);
    B = rows - A;
    C = NaN (B, 1);
    data = [data(:); C];
    data_normalized = [data_normalized(:); C];
    Copy(:,sta) = data;
    Copy_normalized (:,sta) = data_normalized;
end

%% Extract UP/DOWN duration
A = RSP_STAT.Quiet_Wake.Up.Peaks.Cell_averages.duration;
B = RSP_STAT.Motion_Wake.Up.Peaks.Cell_averages.duration;
C = RSP_STAT.NREM.Up.Peaks.Cell_averages.duration;
D = RSP_STAT.REM.Up.Peaks.Cell_averages.duration;

%% Extract Spont. SPikes
status1 = 'rate';
status2 = 'isi';
status3 = 'noSPK';
stage = 'Wake';

Avg_Spikes = DZP_STAT.(stage).Up.Spikes.Cell_averages_OnlySpikes.(status3);
Spikes = DZP_STAT.(stage).Up.Spikes.All_Bouts.(status3);
Rate = DZP_STAT.(stage).Up.Spikes.Cell_averages_OnlySpikes.(status1);
ISI = DZP_STAT.(stage).Up.Spikes.Cell_averages_OnlySpikes.(status2);
%% Extract VM/EEG PSDs
rows = 24;
Copy = NaN (rows,4);
stages = {'Quiet_Wake','Motion_Wake','NREM','REM'};
band = 'isi';
zonditions2 = 'Cell_averages_OnlySpikes';

for sta = 1: numel (stages)
    stage = stages{sta};
    data = RSP_STAT.(stage).Up.Spikes.(zonditions2).(band);
    A = numel (data);
    B = rows - A;
    C = NaN (B, 1);
    data = [data(:); C];
    Copy(:,sta) = data;
end

%% Fig2 - Top plots
stage = 'NREM';
cell = 'x21n05005';
bout = 'b36';

downsampling = 100;

A = RSP_DATA.(cell).(stage).EEG_bout.(bout);
B = RSP_DATA.(cell).(stage).Vm_bout.(bout);

n = numel(A);
EEG_noSAM = floor(n/200*downsampling);
A_data = interp1((1:n)/200,A,(1:EEG_noSAM)/downsampling);
n = numel(B);
Vm_noSAM = floor(n/5000*downsampling);
B_data = interp1((1:n)/5000,B,(1:Vm_noSAM)/downsampling);

close all;
subplot (2,1,1)
plot (A_data);
subplot (2,1,2)
plot (B_data);

%% Fig3s - PFC/RSP PSD comparison
% stages = {'Quiet_Wake','Motion_Wake','NREM','REM'};
stages = {'Wake','NREM','REM'};
band = 'spindle';
zonditions2 = 'Bpower_EEG';
rows = 15;
Copy = NaN (rows,4);

for sta = 1: numel (stages)
    stage = stages{sta};
    if sta == 1
        stagez = 'Motion_Wake';
        data_norm = PFC_STAT.(stagez).(zonditions2).(band);
    else
        stagez = stage;
        data_norm = PFC_STAT.(stage).(zonditions2).(band);
    end   
    data_norm = nanmean (data_norm);
    data = DZP_STAT.(stage).(zonditions2).(band);
    data2 = data./data_norm;
    data2 = (data2*100)-100;
    A = numel (data);
    B = rows - A;
    C = NaN (B,1);
    data = [data(:); C];
    Copy(:,sta) = data;
    Copy_norm(:,sta) = data2;
end


%% Extract Pupil/Whisker Motion
zondition1 = 'pupil';
zondition2 = 'mean';
stages = {'Quiet_Wake','Motion_Wake','NREM','REM'};
rows = 30;
Copy = NaN (rows,4);
Copy_normalized = NaN (rows,4);

for sta = 1: numel (stages)
    stage = stages{sta};
    
    data = PFC_STAT.(stage).(zondition1).(zondition2);
    A = numel (data);
    B = rows - A;
    C = NaN (B, 1);
    data = [data(:); C];
    Copy(:,sta) = data;
end


%% Extract UP_DOWN states, spike timing & Spikes
clear Copy
clear clusters

stage = 'NREM';
evoked = 'Evoked3';

rows = numel (PFC_STAT.(evoked).(stage).spike);

Copy = NaN (rows,1);

if strcmp (evoked, 'Evoked1')
    Copy (:,1) = ones (rows,1)*200;
elseif strcmp (evoked, 'Evoked2')
    Copy (:,1) = ones (rows,1)*400;
elseif strcmp (evoked, 'Evoked3')
    Copy (:,1) = ones (rows,1)*600;
end

data_clusters = PFC_STAT.(evoked).(stage).cluster_type;
data_clusters = data_clusters(:);
Copy (:,3) = PFC_STAT.(evoked).(stage).cluster_ts;
Copy (:,4) = PFC_STAT.(evoked).(stage).spike;

%% Evoked spikes quantification
stage = 'NREM';
evoked = 'Evoked3';

clear Copy

data = DZP_STAT.(evoked).(stage).Average_nrSpikes;
Copy (1,1) = nanmean(data);
n_nonNaN = numel(find (~isnan(data)));
Copy(1,3) = n_nonNaN;
std = nanstd(data,[],2);
Copy(1,2) = std/sqrt(n_nonNaN);


%% Extract & Normalize EEG/VM PSDs
rows = 24;
Copy = NaN (rows,4);
stages = {'Quiet_Wake','Motion_Wake','NREM','REM'};
band = 'slowWave';
zonditions2 = 'Bpower_Vm';

for sta = 1: numel (stages)
    stage = stages{sta};
    data = RSP_STAT.(stage).(zonditions2).(band);
    data_norm = PFC_STAT.(stage).(zonditions2).(band);
    data_norm = nanmean (data_norm);
    
    data2 = (data./data_norm);
    data2 = (data2*100) - 100;
    
    A = numel (data2);
    B = rows - A;
    C = NaN (B, 1);
    data_plot = [data2(:); C];
    Copy(:,sta) = data_plot;
end