%% Manage interpreters
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');
%% Signals to construct
% Case 1 - Closely spaced frequencies.
% NOTE: Cascade through increasing mode orders to observe how amplitude,
% decay factors, energies in the signal converge with appropriate model
% order selection.
clear;
clc;

fs = 10000;
dt = 1/fs;
nSamples = 1000;

f_true = [300, 300, 310, 500, 600, 1000, 2390];

amp_true = [1.5, 0.3, 0.75, 1, 1, 1, 1];

alpha_true = 50 : 10 : 110;
zeta_true = alpha_true ./ (f_true * 2 * pi);

phase_true = deg2rad(0 : 20 : 120);

nSignals = numel(f_true);

[allSignals, outputSignal, timeVector, debug] = ...
    createImpulseResponse(nSignals, nSamples, fs, SNR=0, Frequency=f_true, ...
    Alpha=alpha_true, Amplitude=amp_true, Phase=phase_true);

opts = eraConfig();
opts.returnDebug = true;
opts.poleScaling = false;
svdTolerance = 7;

% Try changing the number of modes to see how the system reconstructs a
% subspace of the original signal.
[eraSignal, fEra, dfEra, ampEra, phaseEra, modesEra, debug] = ...
    eigensystemRealisation(outputSignal, svdTolerance, fs, opts);

poles = debug.Modal.poles;
residues = debug.Modal.residues;

inputImpulse = zeros(1, nSamples);
inputImpulse(1) = 1;

[B, A] = buildFilterCoefficients(poles, residues, outputSignal(1));
parallelFilterSignal = parallelFilter(B, A, inputImpulse);

f1 = figure;
tl = tiledlayout(1, 2);
xlabel(tl, 'Time (s)', 'Interpreter', 'Latex')
ylabel(tl, 'Amplitude', 'Interpreter', 'Latex')

t1 = nexttile;
hold(t1, "on")

t2 = nexttile;
hold(t2, "on")

parallelFilterNMSE = sum((outputSignal(:) - ...
    parallelFilterSignal(:)) .^ 2) / sum(outputSignal(:) .^2);

% t1 — True vs ERA
plot(t1, timeVector, eraSignal,    '-',  'Color', [0.84 0.37 0.00 0.55], ...
     'LineWidth', 3, 'DisplayName', 'ERA signal')
plot(t1, timeVector, outputSignal, '-',  'Color', [0.00 0.45 0.70], ...
     'LineWidth', 1.2, 'DisplayName', 'True signal')
legend(t1, 'Location', 'southeast', 'Box', 'off')
text(0.8, 0.9, "Reconstruction error: " + string(debug.Reconstruction.impulseNMSE))


% t2 — ERA vs Parallel Filter
plot(t2, timeVector, parallelFilterSignal, '-',  'Color', [0.00 0.62 0.45 0.55], ...
     'LineWidth', 3, 'DisplayName', 'Parallel Filter signal')
plot(t2, timeVector, eraSignal,            '-',  'Color', [0.84 0.37 0.00], ...
     'LineWidth', 1.2, 'DisplayName', 'ERA signal')
legend(t2, 'Location', 'southeast', 'Box', 'off')

modeIdx = 1:7;

phaseEra = phaseEra+pi/2;

phase_true = round(phase_true, 2, "significant");
alpha_true = round(alpha_true, 2, "significant");
amp_true = round(amp_true, 2, "significant");
f_true = round(f_true, 2, "significant");

phaseEra = round(phaseEra, 2, "significant");
alpha_true = round(alpha_true, 2, "significant");
amp_true = round(amp_true, 2, "significant");
f_true = round(f_true, 2, "significant");

phaseEra([2, 1]) = phaseEra([1, 2]);
fEra([2, 1]) = fEra([1, 2]);
ampEra([2, 1]) = ampEra([1, 2]);
dfEra([2, 1]) = dfEra([1, 2]);

T = table(modeIdx', ...
    f_true(:),     fEra(:), ...
    zeta_true(:),  dfEra(:), ...
    amp_true(:),   ampEra(:), ...
    phase_true(:), phaseEra(:), ...
    'VariableNames', {'Mode', ...
    'True Frequency','ERA Frequency', ...
    'True Alpha','ERA Alpha', ...
    'True Amp','Era Amp', ...
    'True Phase','ERA Phase'});

format short g
newline
disp(T)


