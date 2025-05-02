clc; clear;
%% Plot Vm during REMs --> Copy/Past from existing traces
REM_traces = [];

%% Save traces for plooting
savePath = 'X:\1 TIDIS Lab\Tiago\Campelo_Manuscript\Fig7s_AllREMtraces\Raw_Data';
saveFile = "dAP5_AllREMtraces.mat";
fullSavePath = fullfile(savePath, saveFile);
save (fullSavePath,'REM_traces');

%% Downsample traces for plotting
REM_traces_interpol =  [];
for a = 1:size (REM_traces,2)
    A = REM_traces(:,a);
    n = numel(A);
    VM_noSAM = floor(n/1000*200);
    A_data = interp1((1:n)/1000,A,(1:VM_noSAM)/200);
    A_data = A_data(:);
    REM_traces_interpol =  [REM_traces_interpol, A_data];
end
%% Plot the downsampled traces
Plotting_nr = size (REM_traces_interpol,2);
figure ()
for a = 1:Plotting_nr
    subplot (6,2,a)
    data2plot = REM_traces_interpol(:,a);
    plot (data2plot)
end