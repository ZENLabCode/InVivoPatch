close all, clc, clear

%LOAD MAT FILES with BASIC ANALYSIS
Path_A = "X:\1 TIDIS Lab\Tiago\Ephys_PROCESSED\In_vivo\Awake_Sleep\PFC_Processed\Pharmacology\NMDAR_Blockers\PFC_NMDAR_SpontaneousSpikes.mat";
ConditionA = load (Path_A); %PFC
Path_B = "X:\1 TIDIS Lab\Tiago\Ephys_PROCESSED\In_vivo\Awake_Sleep\RC_Processed\Pharmacology\NMDAR_Blockers\RC_NMDAR_SpontaneousSpikes.mat";
ConditionB = load (Path_B); %PFC_Benzo
clear Path_**

%% FIGURE I - sleep states comparison for Condition (Basic Analysis)

Data.Condition = 'RC'; %PFC // PFC_Benzo

stages = {'Wake','NREM','REM'};
colorsSTA = [...
    56  77 115;...
    234 174  51;...
    218  93  41;...
    ]/255;

%append corresponding data to structure
for k = 1:numel(stages)
    Data.(stages{k}) = eval(sprintf('%s_%s',Data.Condition,stages{k}));
end

%data type (field) to plot (one firgure per type)
DataTypes = {... {variable, title, ylabel}
    'rate'   , 'Average spontaneous firing rates', 'Spikes (Hz)';...
    'threshold'  , 'Threshold to spike'                   , 'Vm (mV)';...
    'isi', 'Interspike interval'                         , 'secs.';...
    'rate_threshold', 'Spike threshold to firing rates', 'mV'};

figure(3);

DATA_Rate = [];
DATA_threshold = [];

group_Rate = [];
group_threshold = [];

for typ = 1:3
    [datTyp, tit, ylab] = DataTypes{typ,:};
    DATA  = [];
    group = {};
    
    for k = 1:numel(stages)
        %VECTORS DATA, GROUP
        dat = Data.(stages{k}).(datTyp)(:);
        DATA  = [DATA ;dat];
        group = [group;repmat(stages(k),size(dat))];       
    end %loop stages
    
    %Remove NaNs for plotting and stats
    nanIndices = cellfun(@(DATA) isnumeric(DATA) && isnan(DATA), DATA);
    
    DATA = cell2mat (DATA(~nanIndices));
    
    group = group(~nanIndices);
    group = categorical(group);
    
    if typ == 1
        DATA_rate = DATA;
        group_rate = group;
    end
    
    if typ == 2
        DATA_threshold = DATA;
        group_threshold = group;
    end
    
    %STATISTIC
    [p, tbl, stats] = anova1(DATA, group, 'off');
    [c,~,~,gnames]  = multcompare(stats,'display','off');
    for i = 1:size(c, 1)
        fprintf('%s comparison between %s and %s: p = %.4f\n',datTyp,...
            gnames{c(i, 1)}, gnames{c(i, 2)}, c(i, 6));
    end
    
    %PLOT
    subplot(2, 2, typ);
    boxplot(DATA, group,"Colors", colorsSTA, "GroupOrder",stages); hold on
    for i = 1:numel(stages)
        y = DATA(group==stages(i));
        x = i*ones(size(y));
        scatter(x,y, 'Marker', 'o', 'MarkerFaceColor', colorsSTA(i,:), ...
            'MarkerEdgeColor', 'k');
    end
    %text
    title({strrep(Data.Condition,'_','\_'),tit});
    ylabel(ylab);
    
    % Add p-values as text annotations with corresponding lines
    yPos = max(DATA)+0.01*diff(ylim);
    for i = 1:size(c, 1)
        pValue = c(i, 6);
        if pValue < 0.05
            names = gnames(c(i,1:2));
            x = [find(strcmpi(stages,names{1})), find(strcmpi(stages,names{2}))];
            names = cellfun(@(x)strrep(x,'_','\_'),names,'uniformoutput',false);
            line(x, [yPos, yPos], 'Color', 'k', 'LineWidth', 1.5);
            text(mean(x), yPos - 0.08 * (max(ylim) - min(ylim)), sprintf('p = %.4f\n%s - %s', pValue, names{:}), 'Color', 'k', 'HorizontalAlignment', 'center');
        end
    end
end %loop data type


subplot(2, 2, 4);
for sta = 1:length(stages)
    
    stage = stages(sta);
    indices_rate = (group_rate == stage{1,1});
    indices_threshold = (group_threshold == stage{1,1});
    
    DATA_rate_stage = DATA_rate(indices_rate);
    DATA_threshold_stage = DATA_threshold(indices_threshold);
    
    scatter (DATA_threshold_stage,DATA_rate_stage, 'Marker', 'o',...
        'MarkerFaceColor', colorsSTA(sta,:), 'MarkerEdgeColor', 'k');
    hold on
    
    % Calculate the linear fit coefficients
    fit_coefficients = polyfit(DATA_threshold_stage, DATA_rate_stage, 1);
    % Create the linear fit line based on the coefficients
    x_fit = linspace(min(DATA_threshold_stage), max(DATA_threshold_stage), 100);
    y_fit = polyval(fit_coefficients, x_fit);
    
    plot(x_fit, y_fit, 'Color', colorsSTA(sta,:), 'LineWidth', 2);
    
    % Add labels and legend
    ylabel('Spiking rate (Hz)');
    xlabel('Spike threshold (mV)');
%     legend(stages, 'Location', 'Best');
    grid on;
end

%% FIGURE II - effect of Pharmacological intervention (basic analysis) - Wake or NREM

%Conditions to be analyzed
Condition = {'RC', 'RC_NMDAR'};
Data = struct ();
Lin_Fit = struct ();

%SELECT/PREPARE DATA
%Condition to be ploted in Fig. 1
stages = {'Wake','NREM','REM'};

%append corresponding data to structure
for cond = 1:numel(Condition)
    for k = 1:numel(stages)
        Data.(Condition{cond}).(stages{k}) = eval(sprintf('%s_%s',Condition{cond},stages{k}));
    end %stages loop
end %condition loop

colorsSTA = {... {variable, RGB combination}
    'colorSTA1', [56  77 115]/255; ...
    'colorSTA2', [234 174  51]/255;...
    'colorSTA3', [218  93  41]/255; ...
    'colorSTA4', [45  45 45]/255; ...
    };

%data type (field) to plot (one firgure per type)
DataTypes = {... {variable, title, ylabel}
    'rate'   , 'Average spontaneous firing rates', 'Spikes (Hz)';...
    'threshold'  , 'Threshold to spike'                   , 'Vm (mV)';...
    'isi', 'Interspike interval'                         , 'ISI (secs.)';...
    'rate_threshold', 'Spike threshold to firing rates', 'mV'};


for typ = 1:3
    [datTyp, tit, ylab] = DataTypes{typ,:};
    
    for k = 1:numel(stages)
        DATA1  = [];
        group1 = {};
        
        DATA2  = [];
        group2 = {};
        
        figure (k);
        
        %VECTORS DATA, GROUP
        %Condition1
        dat1 = Data.(Condition {1}).(stages{k}).(datTyp)(:);
        DATA1  = [DATA1 ;dat1];
        group1 = [group1; repmat(Condition(1), size(dat1))];
        nanIndices1 = cellfun(@(DATA1) isnumeric(DATA1) && isnan(DATA1), DATA1);
        DATA1 = cell2mat (DATA1(~nanIndices1));
        group1 = group1(~nanIndices1);
        group1 = categorical(group1);
        
        %Condition2
        dat2 = Data.(Condition {2}).(stages{k}).(datTyp)(:);
        DATA2  = [DATA2 ;dat2];
        group2 = [group2; repmat(Condition(2), size(dat2))];
        nanIndices2 = cellfun(@(DATA2) isnumeric(DATA2) && isnan(DATA2), DATA2);
        DATA2 = cell2mat (DATA2(~nanIndices2));
        group2 = group2(~nanIndices2);
        group2 = categorical(group2);
        
        if typ == 1
            Lin_Fit.(Condition {1}).(stages{k}).rate = DATA1;
            Lin_Fit.(Condition {1}).(stages{k}).rate_group = group1;
            Lin_Fit.(Condition {2}).(stages{k}).rate = DATA2;
            Lin_Fit.(Condition {2}).(stages{k}).rate_group = group2;
        end %lin fit loop for figure 4
        
        if typ == 2
            Lin_Fit.(Condition {1}).(stages{k}).threshold = DATA1;
            Lin_Fit.(Condition {1}).(stages{k}).threshold_group = group1;
            Lin_Fit.(Condition {2}).(stages{k}).threshold = DATA2;
            Lin_Fit.(Condition {2}).(stages{k}).threshold_group = group2;
        end %lin fit loop for figure 4
        
        %STATISTIC
        [h, p, ci] = ttest2(DATA1, DATA2, 'Vartype', 'unequal');
        fprintf('%s comparison between %s and %s: p = %.3f\n', datTyp,...
            [Condition{1,1}, ' ', stages{k}], [Condition{1,2}, ' ', stages{k}], p);
        
        %PLOT
        subplot(2, 2, typ);
        DATA = cat (1,DATA1 , DATA2);
        group = cat (1, group1 , group2);
        boxplot (DATA, group, 'Colors', [colorsSTA{k,2}; colorsSTA{4,2}]); hold on;
        
        for a = 1:numel(Condition(1,:))
            
            y = DATA(group ==([Condition{1,a}]));
            x = (ones(size(y)))*a;
            if a == 1
                scatter(x,y, 'Marker', 'o', 'MarkerFaceColor', [colorsSTA{k,2}], 'MarkerEdgeColor', 'k');
            end
            if a ~= 1
                scatter(x,y, 'Marker', 'o', 'MarkerFaceColor', [colorsSTA{4,2}], 'MarkerEdgeColor', 'k');
            end
        end
        
        %Print p-values as text annotations with corresponding lines
        if p<0.05
            %p-value line
            line([1, 2], [max(DATA) + range(DATA) ...
                * 0.05, max(DATA) + range(DATA) * 0.05], 'Color', 'r');
            % y-coordinates for p-values
            yText = max(DATA) + range(DATA) * -0.02;
            % Add p-value text with adjusted y-coordinate
            text(1.5, yText, ['p = ', num2str(p)], 'HorizontalAlignment', 'center');
        end
        
        %text
        title(tit);
        ylabel(ylab);
        
    end %loop sleep stages
    
end %loop data typ

for k = 1: length(stages)
    figure (k);
    typ=4;
    [datTyp, tit, ylab] = DataTypes{typ,:};
    subplot (2,2,4);
    Y = Lin_Fit.(Condition {1}).(stages{k}).rate;
    X = Lin_Fit.(Condition {1}).(stages{k}).threshold;
    scatter (X, Y, ...
        'Marker', 'o','MarkerFaceColor', colorsSTA{k,2}, 'MarkerEdgeColor', 'k'); hold on
    
    % Calculate the linear fit coefficients
    fit_coefficients = polyfit(X, Y, 1);
    % Create the linear fit line based on the coefficients
    x_fit = linspace(min(X), max(X), 100);
    y_fit = polyval(fit_coefficients, x_fit);
    
    plot(x_fit, y_fit, 'Color', colorsSTA{k,2}, 'LineWidth', 2);
    
    Y2 = Lin_Fit.(Condition {2}).(stages{k}).rate;
    X2 = Lin_Fit.(Condition {2}).(stages{k}).threshold;
    scatter (X2, Y2, ...
        'Marker', 'o','MarkerFaceColor', colorsSTA{4,2}, 'MarkerEdgeColor', 'k');
    
    % Calculate the linear fit coefficients
    fit_coefficients = polyfit(X2, Y2, 1);
    % Create the linear fit line based on the coefficients
    x_fit = linspace(min(X2), max(X2), 100);
    y_fit = polyval(fit_coefficients, x_fit);
    
    plot(x_fit, y_fit, 'Color', colorsSTA{4,2}, 'LineWidth', 2);
    
    grid on
    title (tit);
    xlabel ('Spike Threshold');
    ylabel ('Spike rate');
    
end

% subplot(2, 2, 4);
% for sta = 1:length(stages)
%     
%     stage = stages(sta);
%     indices_rate = (group_rate == stage{1,1});
%     indices_threshold = (group_threshold == stage{1,1});
%     
%     DATA_rate_stage = DATA_rate(indices_rate);
%     DATA_threshold_stage = DATA_threshold(indices_threshold);
%     
%     scatter (DATA_threshold_stage,DATA_rate_stage, 'Marker', 'o',...
%         'MarkerFaceColor', colorsSTA(sta,:), 'MarkerEdgeColor', 'k');
%     hold on
%     
%     % Calculate the linear fit coefficients
%     fit_coefficients = polyfit(DATA_threshold_stage, DATA_rate_stage, 1);
%     % Create the linear fit line based on the coefficients
%     x_fit = linspace(min(DATA_threshold_stage), max(DATA_threshold_stage), 100);
%     y_fit = polyval(fit_coefficients, x_fit);
%     
%     plot(x_fit, y_fit, 'Color', colorsSTA(sta,:), 'LineWidth', 2);
%     
%     % Add labels and legend
%     ylabel('Spiking rate (Hz)');
%     xlabel('Spike threshold (mV)');
% %     legend(stages, 'Location', 'Best');
%     grid on;
% end
