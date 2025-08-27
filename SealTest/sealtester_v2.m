% Optimized Seal Test Analysis Script (MATLAB)
% -------------------------------------------------------------
% * Loads a .mat file containing `resampled_data_mV`
% * Down‑samples from 100 kHz to 5 kHz (anti‑alias filtering)
% * Extracts the seal‑test segment (23.4–43 ms)
% * Fits a single‑exponential decay
% * Calculates membrane resistance for V‑steps of 30 mV and 10 mV
% -------------------------------------------------------------

close all; clear; clc;

%% USER‑CONFIGURABLE PARAMETERS -----------------------------------------
filePath        = 'C:\Users\U1306265\Desktop\NewFigures\SealTest\PFC_Baseline/21o01000_IN0.mat';
startTime_s     = 0.0234;   % seal‑test window start (s)
endTime_s       = 0.0430;   % seal‑test window end   (s)
Vsteps_mV       = [30, 10]; % voltage steps to evaluate (mV)

origFS_Hz       = 100000;   % original sampling rate (Hz)
targetFS_Hz     = 5000;     % desired sampling rate  (Hz)

%% ----------------------------------------------------------------------
% Derived parameters
DSfactor = origFS_Hz / targetFS_Hz;
assert(mod(DSfactor,1)==0, 'Down‑sample factor must be integer.');

%% Load data -------------------------------------------------------------
S = load(filePath);     % assumes variable `resampled_data_mV`
if ~isfield(S,'resampled_data_mV')
    error('Variable "resampled_data_mV" not found in %s',filePath);
end
signal_mV = S.resampled_data_mV(:);           % force column vector
numSamples = numel(signal_mV);

%% Build original time vector -------------------------------------------
t_original = (0:numSamples-1).' / origFS_Hz;   % column vector (s)

%% Down‑sample with anti‑alias filter -----------------------------------
signal_ds = decimate(signal_mV, DSfactor);     % antialiasing FIR inside decimate
t_ds      = t_original(1:DSfactor:end);        % keep synchronous time base

%% Extract seal‑test segment -------------------------------------------
idxWin = t_ds >= startTime_s & t_ds <= endTime_s;
if nnz(idxWin) < 10
    error('Less than 10 samples in seal‑test window – adjust parameters.');
end
t_seal = t_ds(idxWin);
y_seal = signal_ds(idxWin);

% Shift time so window starts at t=0 for numerical stability

t0 = t_seal - t_seal(1);

%% Exponential fit -------------------------------------------------------
expFun = @(p,x) p(1).*exp(p(2).*x) + p(3);    % p = [a b c]

% Initial guess heuristics
p0 = [max(y_seal)-min(y_seal), -100, min(y_seal)];
opts = optimoptions('lsqcurvefit', 'Display','off');
[pFit, resnorm] = lsqcurvefit(expFun, p0, t0, y_seal, [], [], opts);

%% Plot data & fit -------------------------------------------------------
figure('Name','Seal‑test Fit','Color','w');
plot(t0*1e3, y_seal, 'b', 'DisplayName','Data'); hold on;
plot(t0*1e3, expFun(pFit, t0), 'r--', 'LineWidth',1.2, 'DisplayName','Exp fit');
xlabel('Time (ms)'); ylabel('Signal (mV)');
legend('Location','best'); grid on;
title('Seal‑test exponential fit');

%% Membrane resistance calculation --------------------------------------
I_steady_pA = pFit(3);               % fitted steady‑state current (units assumed pA)

% Convert V (mV -> V) and I (pA -> A) then compute R (Ohm)
R_MOhm = (Vsteps_mV(:) * 1e-3) ./ (I_steady_pA * 1e-12) / 1e6; % (vector, MΩ)

%% Display results -------------------------------------------------------
fprintf('\nSteady‑state current  (pA): %.3f\n', I_steady_pA);
for k = 1:numel(Vsteps_mV)
    fprintf('V_step = %2d mV  -->  R_m = %.2f MΩ\n', Vsteps_mV(k), R_MOhm(k));
end

%% (Optional) Return results to workspace
results.R_MOhm       = R_MOhm; %#ok<STRNU>
results.Vsteps_mV    = Vsteps_mV;
results.I_steady_pA  = I_steady_pA;
results.pFit         = pFit;
results.resnorm      = resnorm;
assignin('base','sealtest_results',results);

% End of script
