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

Fs = 10000;
dt = 1/Fs;
Ns = 1000;
t = (0:(Ns-1))*dt;
frequencyAxis = (0 : (Ns / 2)) * (Fs / Ns);

f_true = [300, 300, 310, 500, 600, 1000, 2390];

amp_true = [1.5, 0.3, 0.75, 1, 1, 1, 1];

alpha_true = 50 : 10 : 110;
zeta_true = alpha_true ./ (f_true * 2 * pi);

phase_true = deg2rad(0 : 10 : 60);

M = length(f_true);
sigs = zeros(M,Ns);

for m = 1 : M
    sigs(m, :) = amp_true(m) * exp(-alpha_true(m) * t) .* ...
        sin(2 * pi * f_true(m) * t + phase_true(m));
    modeEnergyTrue(m) = sum(sigs(m, :) .^ 2);
end
modeIndex = sort(modeEnergyTrue, 'descend');
sig = sum(sigs);

opts = eraConfig();
opts.returnDebug = true;
opts.poleScaling = true;
svdTolerance = 40;

% Try changing the number of modes to see how the system reconstructs a
% subspace of the original signal.
[outputSignal, f_era, zeta_era, amp_era, phase_era, modes_era, debug] = ...
    eigensystemRealisation(sig, svdTolerance, Fs, opts);

poles = diag(debug.Modal.At);
residues = debug.Modal.Ct' .* debug.Modal.Bt;

inputImpulse = zeros(1, Ns);
inputImpulse(1) = 1;
[B, A] = buildFilterCoefficients(poles, residues, sig(1));
parallelFilterSignal = parallelFilter(B, A, inputImpulse);

% --- Match ERA modes to true modes by nearest frequency ---
n = numel(f_true);
idx = zeros(n,1);
for k = 1:n
    [~, idx(k)] = min(abs(debug.Reconstruction.modeEnergies - modeEnergyTrue(k)));
end

f_e     = f_era(idx);
zeta_e  = zeta_era(idx);
amp_e   = amp_era(idx);
phase_e = wrapToPi(phase_era(idx) + pi/2);

debug.Reconstruction.modalImpulses = debug.Reconstruction.modalImpulses(idx, :);
debug.Modal.decayRates = debug.Modal.decayRates(idx);

f = figure;
til = tiledlayout(3, 3);

tf = nexttile(1);
title('7 Individual impulse responses', 'Interpreter', 'latex')
hold(tf, "on")
xlabel('Time (ms)', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
ylim([-3 3])

tf2 = nexttile(2);
title('Frequency response of individual impulse responses', 'Interpreter', 'latex')
xlabel('Frequency (kHz)', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
hold(tf2, "on")
ylim([0 0.4])

tf3 = nexttile(4);
title("Reconstructed impulse responses with a singular value tolerance " + ...
    "of "  +  string(svdTolerance), 'Interpreter', 'latex')
hold(tf3, "on")
xlabel('Time (ms)', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
ylim([-3 3])

tf4 = nexttile(5);
title('Frequency response of reconstructed individual impulse responses', 'Interpreter', 'latex')
xlabel('Frequency (kHz)', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
hold(tf4, "on")
ylim([0 0.4])

for m = 1 : M
    plot(tf, t/1000, sigs(m, :), 'LineWidth', 1.1, 'DisplayName', ...
        "Frequency: " + string(f_true(m)))


    frequencyResponse = fft(sigs(m, :));
    P2 = abs(frequencyResponse / Ns);
    P1 = P2(1:Ns/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    plot(tf2, frequencyAxis / 1000, P1, 'LineWidth', 1.1, 'DisplayName', ...
        "Frequency: " + string(f_true(m)))
end

for m = 1 : M
    plot(tf3, t/1000, debug.Reconstruction.modalImpulses(m, :), 'LineWidth', 1.1, 'DisplayName', ...
        "Frequency: " + string(f_true(m)))

    frequencyResponse = fft(debug.Reconstruction.modalImpulses(m, :));
    P2 = abs(frequencyResponse / Ns);
    P1 = P2(1:Ns/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    plot(tf4, frequencyAxis / 1000, P1, 'LineWidth', 1.1, 'DisplayName', ...
        "Frequency: " + string(f_true(m)))

end

legend(tf)
legend(tf2)
legend(tf3)
legend(tf4)

nexttile(7)
title(['Composite impulse response and reconstructed impulse response' ...
    ' using ERA'], 'Interpreter', 'latex')
hold on
plot(t / 1000, sig, 'k', 'LineWidth', 1.5, 'DisplayName', 'Original signal')
plot(t / 1000, outputSignal, '--r', 'LineWidth', 1.5, 'DisplayName', ...
    'Reconstructed signal')
plot(t / 1000, parallelFilterSignal, '-.b', 'LineWidth', 1.4, ...
    'DisplayName', 'Parallel Filter signal')
legend
xlabel('Time (ms)', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
reconError = sprintf("NMSE: %.4d", debug.Reconstruction.impulseNMSE);
text(0.7, 0.3, reconError, 'Units', 'normalized')

nexttile(8)
hold on
title(['Associated frequency response and reconstructed frequency ' ...
    'response'], 'Interpreter', 'latex')

frequencyResponse = fft(sig);
P2 = abs(frequencyResponse / Ns);
P1 = P2(1:Ns/2+1);
P1(2:end-1) = 2*P1(2:end-1);

plot(frequencyAxis / 1000, P1, 'k', 'LineWidth', 1.3, 'DisplayName', ...
    'Composite Frequency Response')

frequencyResponse = fft(outputSignal);
P2 = abs(frequencyResponse / Ns);
P1 = P2(1:Ns/2+1);
P1(2:end-1) = 2*P1(2:end-1);

plot(frequencyAxis / 1000, P1, '--r', 'LineWidth', 1.3, 'DisplayName', ...
    'Reconstructed Frequency Response')
xlabel('Frequency (kHz)', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
title('One sided frequency response of impulse function', 'Interpreter', 'latex')
grid;

nexttile([3, 1])
semilogy(f_true, alpha_true, 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, ...
    'DisplayName', ...
    'Target decay factors per frequency')
hold on
semilogy(f_true, debug.Modal.decayRates, 'rx', 'MarkerSize', 10, ...
    'LineWidth', 1.5, 'DisplayName', 'Recalculated decay factors')
xlabel('Frequency (Hz)', 'Interpreter', 'latex')
ylabel('Decay rate - $\alpha: (s^{-1})$', 'Interpreter', 'latex')
xlim([0 2.5E3])
ylim([ 0 150])
legend
title(['Comparison of decay rate for simulated signals ' ...
    'vs. ERA reconstruction'], 'Interpreter', 'latex')


modeIdx = 1 : 7;
figure
tiledlayout(2, 2)
nexttile
semilogy(modeIdx, f_true, 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, ...
    'DisplayName', ...
    'Target frequency')
hold on
semilogy(modeIdx, f_e, 'rx', 'MarkerSize', 10, ...
    'LineWidth', 1.5, 'DisplayName', 'Extracted Frequency')
xlabel('Mode Index')
ylabel('Frequency (Hz)')
legend
title(['Comparison of Frequency for Simulated and ERA extracted values'])

nexttile
semilogy(modeIdx, zeta_true, 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, ...
    'DisplayName', ...
    'Target decay factor')
hold on
semilogy(modeIdx, zeta_e, 'rx', 'MarkerSize', 10, ...
    'LineWidth', 1.5, 'DisplayName', 'Extracted decay factor')
xlabel('Mode Index')
ylabel('Nepers')
legend
title(['Comparison of decay factor for Simulated and ERA extracted values'])

nexttile
plot(modeIdx, amp_true, 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, ...
    'DisplayName', ...
    'Target amplitude')
hold on
plot(modeIdx, amp_e, 'rx', 'MarkerSize', 10, ...
    'LineWidth', 1.5, 'DisplayName', 'Extracted amplitude')
xlabel('Mode Index')
ylabel('Amplitude')
legend
title(['Comparison of amplitude for Simulated and ERA extracted values'])

nexttile
polarplot(phase_true, 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, ...
    'DisplayName', ...
    'Target phase')
hold on
polarplot(phase_e, 'rx', 'MarkerSize', 10, ...
    'LineWidth', 1.5, 'DisplayName', 'Extracted phase')
legend
title(['Comparison of phase for Simulated and ERA extracted values'])
 
% Errors
f_err     = round(100 * (f_e(:) - f_true(:)) ./ f_true(:), 2);
zeta_err  = round(100 * (zeta_e(:) - zeta_true(:)) ./ zeta_true(:), 2);
amp_err   = round(100 * (amp_e(:)  - amp_true(:))  ./ amp_true(:), 2);
phase_err = round(wrapToPi(phase_e(:) - phase_true(:)), 2);   % rad
 
% Table

T = table(modeIdx', ...
    f_true(:),     f_e(:),     f_err, ...
    zeta_true(:),  zeta_e(:),  zeta_err, ...
    amp_true(:),   amp_e(:),   amp_err, ...
    phase_true(:), phase_e(:), phase_err, ...
    'VariableNames', {'Mode', ...
    'Freq_true_Hz','Freq_ERA_Hz','Freq_err_pct', ...
    'Alpha_true','Alpha_ERA','Alpha_err_pct', ...
    'Amp_true','Amp_ERA','Amp_err_pct', ...
    'Phase_true_rad','Phase_ERA_rad','Phase_err_rad'});
 
format short g
disp(T)

H1 = debug.Hankel.shiftedHankelMatrix;
H0 = debug.Hankel.hankelMatrix;
nColumns = debug.Hankel.nColumns;
nRows = debug.Hankel.nRows;
m = svdTolerance;

[U1, S1, V1] = svd(H1);

Um1 = U1(:, 1 : 2 * m);
Sm1 = S1(1 : 2 * m, 1 : 2 * m);
Vm1 = V1(:, 1 : 2 * m);

% Half power and inverse matrices. Transform to avoid computing over the
% entire S matrix. Only saves time for large S (assumed).
singularValues  = diag(Sm1);
sHalf = diag(sqrt(singularValues));
sHalfInv = diag(1 ./ sqrt(singularValues)); 

% Eq. (22-24).
A1 = sHalfInv * Vm1' * H0 * Um1 * sHalfInv;

[psi, lambda] = eig(A1);

% Discrete to continuous time conversion
continuousEigenValues = log(lambda) ./ dt;
continuousEigenValues = diag(continuousEigenValues);

naturalFrequencies = abs(continuousEigenValues);

% Extract frequencies from the eigenvalues.
shiftedFrequencies = naturalFrequencies ./ (2*pi);

% Extract the decay rate.
shiftedDecayRates = -real(continuousEigenValues);

% Extract damping factor from the eigenvalues.
shiftedDampingFactors = -shiftedDecayRates ./ naturalFrequencies;

figure
plot(f_true, zeta_true, 'rsq', 'LineWidth', 1.4, 'DisplayName', ...
    'True eigenvalues', 'MarkerSize', 10)
hold on
plot(f_era, zeta_era, 'ko', 'LineWidth', 1.4, ...
    'DisplayName', 'Standard frequencies')
plot(shiftedFrequencies, shiftedDampingFactors, 'gx', 'LineWidth', 1.4, ...
    'DisplayName', 'Shifted Frequencies')
ylabel('Damping factor')
xlabel('Frequency')
legend

msi = modeSimilarityIndex(debug.Modal.allFrequencies, shiftedFrequencies, ...
    debug.Modal.allDampingFactors, shiftedDampingFactors);

figure
bar(1 : length(msi), msi)
xlabel('mode index')
ylabel('msi')

%% Case 2 - Signals with noise.
clear;
clc;

Fs = 10000;
dt = 1/Fs;
Ns = 1000;
t = (0:(Ns-1))*dt;
frequencyAxis = (0 : (Ns / 2)) * (Fs / Ns);

SNR = 20;

f_true = [500, 1000, 1500];

amp_true = [1.5, 0.3, 0.75];

alpha_true = 50 : 30 : 110;
zeta_true = alpha_true ./ (f_true * 2 * pi);

phase_true = wrapToPi(deg2rad(0 : 30 : 60));

M = length(f_true);
sigs = zeros(M,Ns);

for m = 1 : M
    sigs(m, :) = amp_true(m) * exp(-alpha_true(m) * t) .* ...
        sin(2 * pi * f_true(m) * t + phase_true(m));
        modeEnergyTrue(m) = sum(sigs(m, :) .^ 2);
end

sig = sum(sigs);

%%%%%%%%%%%%%%%% ADD NOISE %%%%%%%%%%%%%%%%%%%%%%
Es = sum(sig.^2);
En = Es/(10^(SNR/10));
sign0 = 1 - 2*rand(1,Ns);
En0 = sum(sign0.^2);
an = sqrt(En/En0);
sign = an*sign0;
sig = sig + sign;

opts = eraConfig();
opts.returnDebug = true;
opts.poleScaling = true;
svdTolerance = 10;

% Try changing the number of modes to see how the system reconstructs a
% subspace of the original signal.
[outputSignal, f_era, zeta_era, amp_era, phase_era, modes_era, debug] = ...
    eigensystemRealisation(sig, svdTolerance, Fs, opts);

% --- Match ERA modes to true modes by nearest frequency ---
n = numel(f_true);
idx = zeros(n,1);
for k = 1:n
    [~, idx(k)] = min(abs(debug.Reconstruction.modeEnergies - modeEnergyTrue(k)));
end

f_e     = f_era(idx);
zeta_e  = zeta_era(idx);
amp_e   = amp_era(idx);
phase_e = wrapToPi(phase_era(idx) + pi/2);

debug.Reconstruction.modalImpulses = debug.Reconstruction.modalImpulses(idx, :);
debug.Modal.decayRates = debug.Modal.decayRates(idx);

f = figure;
til = tiledlayout(3, 3);

tf = nexttile(1);
title([num2str(M) ' Individual impulse responses'], 'Interpreter', 'latex')
hold(tf, "on")
xlabel('Time (ms)', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
ylim([-3 3])

tf2 = nexttile(2);
title('Frequency response of individual impulse responses', 'Interpreter', 'latex')
xlabel('Frequency (kHz)', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
hold(tf2, "on")
ylim([0 0.4])

tf3 = nexttile(4);
title("Reconstructed impulse responses with a singular value tolerance " + ...
    "of "  +  string(svdTolerance), 'Interpreter', 'latex')
hold(tf3, "on")
xlabel('Time (ms)', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
ylim([-3 3])

tf4 = nexttile(5);
title('Frequency response of reconstructed individual impulse responses', 'Interpreter', 'latex')
xlabel('Frequency (kHz)', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
hold(tf4, "on")
ylim([0 0.4])

for m = 1 : M
    plot(tf, t/1000, sigs(m, :), 'LineWidth', 1.1, 'DisplayName', ...
        "Frequency: " + string(f_true(m)))


    frequencyResponse = fft(sigs(m, :));
    P2 = abs(frequencyResponse / Ns);
    P1 = P2(1:Ns/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    plot(tf2, frequencyAxis / 1000, P1, 'LineWidth', 1.1, 'DisplayName', ...
        "Frequency: " + string(f_true(m)))
end

for m = 1 : M
    plot(tf3, t/1000, debug.Reconstruction.modalImpulses(m, :), 'LineWidth', 1.1, 'DisplayName', ...
        "Frequency: " + string(f_true(m)))

    frequencyResponse = fft(debug.Reconstruction.modalImpulses(m, :));
    P2 = abs(frequencyResponse / Ns);
    P1 = P2(1:Ns/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    plot(tf4, frequencyAxis / 1000, P1, 'LineWidth', 1.1, 'DisplayName', ...
        "Frequency: " + string(f_true(m)))

end

legend(tf)
legend(tf2)
legend(tf3)
legend(tf4)

nexttile(7)
title(['Composite impulse response and reconstructed impulse response' ...
    ' using ERA'], 'Interpreter', 'latex')
hold on
plot(t / 1000, sig, 'k', 'LineWidth', 1.5, 'DisplayName', 'Original signal')
plot(t / 1000, outputSignal, '--r', 'LineWidth', 1.5, 'DisplayName', ...
    'Reconstructed signal')
legend
xlabel('Time (ms)', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
reconError = sprintf("NMSE: %.4d", debug.Reconstruction.impulseNMSE);
text(0.7, 0.3, reconError, 'Units', 'normalized')

nexttile(8)
hold on
title(['Associated frequency response and reconstructed frequency ' ...
    'response'], 'Interpreter', 'latex')

frequencyResponse = fft(sig);
P2 = abs(frequencyResponse / Ns);
P1 = P2(1:Ns/2+1);
P1(2:end-1) = 2*P1(2:end-1);

plot(frequencyAxis / 1000, P1, 'k', 'LineWidth', 1.3, 'DisplayName', ...
    'Composite Frequency Response')

frequencyResponse = fft(outputSignal);
P2 = abs(frequencyResponse / Ns);
P1 = P2(1:Ns/2+1);
P1(2:end-1) = 2*P1(2:end-1);

plot(frequencyAxis / 1000, P1, '--r', 'LineWidth', 1.3, 'DisplayName', ...
    'Reconstructed Frequency Response')
xlabel('Frequency (kHz)', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
title('One sided frequency response of impulse function', 'Interpreter', 'latex')
grid;

nexttile([3, 1])
semilogy(f_true, alpha_true, 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, ...
    'DisplayName', ...
    'Target decay factors per frequency')
hold on
semilogy(f_true, debug.Modal.decayRates, 'rx', 'MarkerSize', 10, ...
    'LineWidth', 1.5, 'DisplayName', 'Recalculated decay factors')
xlabel('Frequency (Hz)', 'Interpreter', 'latex')
ylabel('Decay rate - $\alpha: (s^{-1})$', 'Interpreter', 'latex')
xlim([0 2.5E3])
ylim([ 0 150])
legend
title(['Comparison of decay rate for simulated signals ' ...
    'vs. ERA reconstruction'], 'Interpreter', 'latex')


modeIdx = 1 : m;
figure
tiledlayout(2, 2)
nexttile
semilogy(modeIdx, f_true, 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, ...
    'DisplayName', ...
    'Target frequency')
hold on
semilogy(modeIdx, f_e, 'rx', 'MarkerSize', 10, ...
    'LineWidth', 1.5, 'DisplayName', 'Extracted Frequency')
xlabel('Mode Index')
ylabel('Frequency (Hz)')
legend
title(['Comparison of Frequency for Simulated and ERA extracted values'])

nexttile
semilogy(modeIdx, zeta_true, 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, ...
    'DisplayName', ...
    'Target decay factor')
hold on
semilogy(modeIdx, zeta_e, 'rx', 'MarkerSize', 10, ...
    'LineWidth', 1.5, 'DisplayName', 'Extracted decay factor')
xlabel('Mode Index')
ylabel('Nepers')
legend
title(['Comparison of decay factor for Simulated and ERA extracted values'])

nexttile
plot(modeIdx, amp_true, 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, ...
    'DisplayName', ...
    'Target amplitude')
hold on
plot(modeIdx, amp_e, 'rx', 'MarkerSize', 10, ...
    'LineWidth', 1.5, 'DisplayName', 'Extracted amplitude')
xlabel('Mode Index')
ylabel('Amplitude')
legend
title(['Comparison of amplitude for Simulated and ERA extracted values'])

nexttile
polarplot(phase_true, 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, ...
    'DisplayName', ...
    'Target phase')
hold on
polarplot(phase_e, 'rx', 'MarkerSize', 10, ...
    'LineWidth', 1.5, 'DisplayName', 'Extracted phase')
legend
title(['Comparison of phase for Simulated and ERA extracted values'])
 
% Errors
f_err     = round(100 * (f_e(:) - f_true(:)) ./ f_true(:), 2);
zeta_err  = round(100 * (zeta_e(:) - zeta_true(:)) ./ zeta_true(:), 2);
amp_err   = round(100 * (amp_e(:)  - amp_true(:))  ./ amp_true(:), 2);
phase_err = round(wrapToPi(phase_e(:) - phase_true(:)), 2);   % rad
 
% Table

T = table(modeIdx', ...
    f_true(:),     f_e(:),     f_err, ...
    zeta_true(:),  zeta_e(:),  zeta_err, ...
    amp_true(:),   amp_e(:),   amp_err, ...
    phase_true(:), phase_e(:), phase_err, ...
    'VariableNames', {'Mode', ...
    'Freq_true_Hz','Freq_ERA_Hz','Freq_err_pct', ...
    'Alpha_true','Alpha_ERA','Alpha_err_pct', ...
    'Amp_true','Amp_ERA','Amp_err_pct', ...
    'Phase_true_rad','Phase_ERA_rad','Phase_err_rad'});
 
format short g
disp(T)

H1 = debug.Hankel.shiftedHankelMatrix;
H0 = debug.Hankel.hankelMatrix;
nColumns = debug.Hankel.nColumns;
nRows = debug.Hankel.nRows;
m = svdTolerance;

[U1, S1, V1] = svd(H1);

Um1 = U1(:, 1 : 2 * m);
Sm1 = S1(1 : 2 * m, 1 : 2 * m);
Vm1 = V1(:, 1 : 2 * m);

% Half power and inverse matrices. Transform to avoid computing over the
% entire S matrix. Only saves time for large S (assumed).
singularValues  = diag(Sm1);
sHalf = diag(sqrt(singularValues));
sHalfInv = diag(1 ./ sqrt(singularValues)); 

% Eq. (22-24).
A1 = sHalfInv * Um1' * H0 * Vm1 * sHalfInv;

[psi, lambda] = eig(A1);

% Discrete to continuous time conversion
continuousEigenValues = log(lambda) ./ dt;

naturalFrequencies = abs(continuousEigenValues);

% Extract frequencies from the eigenvalues.
shiftedFrequencies = naturalFrequencies ./ (2*pi);

% Extract the decay rate.
shiftedDecayRates = -real(continuousEigenValues);

% Extract damping factor from the eigenvalues.
shiftedDampingFactors = -real(continuousEigenValues) ./ ...
    naturalFrequencies;

figure
plot(f_era, zeta_era, 'o', 'LineWidth', 1.4, ...
    'DisplayName', 'Standard frequencies')
hold on
plot(shiftedFrequencies, -shiftedDampingFactors, 'gx', 'LineWidth', 1.4, ...
    'DisplayName', 'Shifted Frequencies')
ylabel('Damping factor')
xlabel('Frequency')
legend

%% Noisy signals, more modes
clear;
clc;

SNR = 20;

Fs = 10000;
dt = 1/Fs;
Ns = 1000;
t = (0:(Ns-1))*dt;
frequencyAxis = (0 : (Ns / 2)) * (Fs / Ns);

fmin = 20;
fmax = 500;
M = 15;
f_true = fmin + (fmax - fmin)*rand(1,M);
alpha_true = 10 + 10*rand(1,M);
zeta_true = alpha_true ./ (f_true * 2 * pi);
amp_true = rand(1,M);
phase_true = wrapToPi(deg2rad(-180 + 360*rand(1,M)));

%%%%%%%% BUILD TARGET SIGNAL %%%%
M = length(f_true);
sigs = zeros(M,Ns);
for m=1:M
    sigs(m, :) = amp_true(m)*exp(-alpha_true(m)*t).*sin(2*pi*f_true(m)*t + phase_true(m));
    modeEnergyTrue(m) = sum(sigs(m, :) .^ 2);
end

sig = sum(sigs);

%%%%%%%% ADD NOISE %%%%%%%%%%%%%%%%%%%%%%
Es = sum(sig.^2);
En = Es/(10^(SNR/10));
sign0 = 1 - 2*rand(1,Ns);
En0 = sum(sign0.^2);
an = sqrt(En/En0);
sign = an*sign0;
sig = sig + sign;

opts = eraConfig();
opts.returnDebug = true;
opts.poleScaling = true;
opts.nColumns = [];
opts.nRows = [];
svdTolerance = 30;

% Try changing the number of modes to see how the system reconstructs a
% subspace of the original signal.
[outputSignal, f_era, zeta_era, amp_era, phase_era, modes_era, debug] = ...
    eigensystemRealisation(sig, svdTolerance, Fs, opts);

% --- Match ERA modes to true modes by nearest frequency ---
n = numel(f_true);
idx = zeros(n,1);
for k = 1:n
    [~, idx(k)] = min(abs(debug.Reconstruction.modeEnergies - modeEnergyTrue(k)));
end

f_e     = f_era(idx);
zeta_e  = zeta_era(idx);
amp_e   = amp_era(idx);
phase_e = wrapToPi(phase_era(idx) + pi/2);

debug.Reconstruction.modalImpulses = debug.Reconstruction.modalImpulses(idx, :);
debug.Modal.decayRates = debug.Modal.decayRates(idx);

f = figure;
til = tiledlayout(3, 3);

tf = nexttile(1);
title([num2str(M) ' Individual impulse responses'], 'Interpreter', 'latex')
hold(tf, "on")
xlabel('Time (ms)', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
ylim([-3 3])

tf2 = nexttile(2);
title('Frequency response of individual impulse responses', 'Interpreter', 'latex')
xlabel('Frequency (kHz)', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
hold(tf2, "on")
ylim([0 0.5])
xlim([0 0.6])

tf3 = nexttile(4);
title("Reconstructed impulse responses with a singular value tolerance " + ...
    "of "  +  string(svdTolerance), 'Interpreter', 'latex')
hold(tf3, "on")
xlabel('Time (ms)', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
ylim([-3 3])

tf4 = nexttile(5);
title('Frequency response of reconstructed individual impulse responses', 'Interpreter', 'latex')
xlabel('Frequency (kHz)', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
hold(tf4, "on")
ylim([0 0.5])
xlim([0 0.6])

for m = 1 : M
    plot(tf, t/1000, sigs(m, :), 'LineWidth', 1.1, 'DisplayName', ...
        "Frequency: " + string(f_true(m)))


    frequencyResponse = fft(sigs(m, :));
    P2 = abs(frequencyResponse / Ns);
    P1 = P2(1:Ns/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    plot(tf2, frequencyAxis / 1000, P1, 'LineWidth', 1.1, 'DisplayName', ...
        "Frequency: " + string(f_true(m)))
end

for m = 1 : M
    plot(tf3, t/1000, debug.Reconstruction.modalImpulses(m, :), 'LineWidth', 1.1, 'DisplayName', ...
        "Frequency: " + string(f_true(m)))

    frequencyResponse = fft(debug.Reconstruction.modalImpulses(m, :));
    P2 = abs(frequencyResponse / Ns);
    P1 = P2(1:Ns/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    plot(tf4, frequencyAxis / 1000, P1, 'LineWidth', 1.1, 'DisplayName', ...
        "Frequency: " + string(f_true(m)))

end

nexttile(7)
title(['Composite impulse response and reconstructed impulse response' ...
    ' using ERA'], 'Interpreter', 'latex')
hold on
plot(t / 1000, sig, 'k', 'LineWidth', 1.5, 'DisplayName', 'Original signal')
plot(t / 1000, outputSignal, '--r', 'LineWidth', 1.5, 'DisplayName', ...
    'Reconstructed signal')
legend
xlabel('Time (ms)', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
reconError = sprintf("NMSE: %.4d", debug.Reconstruction.impulseNMSE);
text(0.7, 0.3, reconError, 'Units', 'normalized')

nexttile(8)
hold on
title(['Associated frequency response and reconstructed frequency ' ...
    'response'], 'Interpreter', 'latex')

frequencyResponse = fft(sig);
P2 = abs(frequencyResponse / Ns);
P1 = P2(1:Ns/2+1);
P1(2:end-1) = 2*P1(2:end-1);

plot(frequencyAxis / 1000, P1, 'k', 'LineWidth', 1.3, 'DisplayName', ...
    'Composite Frequency Response')
xlim([0 0.6])

frequencyResponse = fft(outputSignal);
P2 = abs(frequencyResponse / Ns);
P1 = P2(1:Ns/2+1);
P1(2:end-1) = 2*P1(2:end-1);

plot(frequencyAxis / 1000, P1, '--r', 'LineWidth', 1.3, 'DisplayName', ...
    'Reconstructed Frequency Response')
xlabel('Frequency (kHz)', 'Interpreter', 'latex')
ylabel('Amplitude', 'Interpreter', 'latex')
title('One sided frequency response of impulse function', 'Interpreter', 'latex')
grid;
xlim([0 0.6])

nexttile([3, 1])
semilogy(f_true, alpha_true, 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, ...
    'DisplayName', ...
    'Target decay factors per frequency')
hold on
semilogy(f_true, debug.Modal.decayRates, 'rx', 'MarkerSize', 10, ...
    'LineWidth', 1.5, 'DisplayName', 'Recalculated decay factors')
xlabel('Frequency (Hz)', 'Interpreter', 'latex')
ylabel('Decay rate - $\alpha: (s^{-1})$', 'Interpreter', 'latex')
xlim([0 600])
ylim([0 150])
legend
title(['Comparison of decay rate for simulated signals ' ...
    'vs. ERA reconstruction'], 'Interpreter', 'latex')


modeIdx = 1 : m;
figure
tiledlayout(2, 2)
nexttile
semilogy(modeIdx, f_true, 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, ...
    'DisplayName', ...
    'Target frequency')
hold on
semilogy(modeIdx, f_e, 'rx', 'MarkerSize', 10, ...
    'LineWidth', 1.5, 'DisplayName', 'Extracted Frequency')
xlabel('Mode Index')
ylabel('Frequency (Hz)')
legend
title(['Comparison of Frequency for Simulated and ERA extracted values'])

nexttile
semilogy(modeIdx, zeta_true, 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, ...
    'DisplayName', ...
    'Target decay factor')
hold on
semilogy(modeIdx, zeta_e, 'rx', 'MarkerSize', 10, ...
    'LineWidth', 1.5, 'DisplayName', 'Extracted decay factor')
xlabel('Mode Index')
ylabel('Nepers')
legend
title(['Comparison of decay factor for Simulated and ERA extracted values'])

nexttile
plot(modeIdx, amp_true, 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, ...
    'DisplayName', ...
    'Target amplitude')
hold on
plot(modeIdx, amp_e, 'rx', 'MarkerSize', 10, ...
    'LineWidth', 1.5, 'DisplayName', 'Extracted amplitude')
xlabel('Mode Index')
ylabel('Amplitude')
legend
title(['Comparison of amplitude for Simulated and ERA extracted values'])

nexttile
polarplot(phase_true, 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, ...
    'DisplayName', ...
    'Target phase')
hold on
polarplot(phase_e, 'rx', 'MarkerSize', 10, ...
    'LineWidth', 1.5, 'DisplayName', 'Extracted phase')
legend
title(['Comparison of phase for Simulated and ERA extracted values'])
 
% Errors
f_err     = round(100 * (f_e(:) - f_true(:)) ./ f_true(:), 2);
zeta_err  = round(100 * (zeta_e(:) - zeta_true(:)) ./ zeta_true(:), 2);
amp_err   = round(100 * (amp_e(:)  - amp_true(:))  ./ amp_true(:), 2);
phase_err = round(wrapToPi(phase_e(:) - phase_true(:)), 2);   % rad
 
% Table

T = table(modeIdx', ...
    f_true(:),     f_e(:),     f_err, ...
    zeta_true(:),  zeta_e(:),  zeta_err, ...
    amp_true(:),   amp_e(:),   amp_err, ...
    phase_true(:), phase_e(:), phase_err, ...
    'VariableNames', {'Mode', ...
    'Freq_true_Hz','Freq_ERA_Hz','Freq_err_pct', ...
    'Alpha_true','Alpha_ERA','Alpha_err_pct', ...
    'Amp_true','Amp_ERA','Amp_err_pct', ...
    'Phase_true_rad','Phase_ERA_rad','Phase_err_rad'});

meanErrors = grpstats(T, [], "mean", "DataVars", ["Freq_err_pct", "Alpha_err_pct", ...
    "Amp_err_pct", "Phase_err_rad"])
 
format('default')
disp(T)

H1 = debug.Hankel.shiftedHankelMatrix;
H0 = debug.Hankel.hankelMatrix;
nColumns = debug.Hankel.nColumns;
nRows = debug.Hankel.nRows;
m = svdTolerance;

[U1, S1, V1] = svd(H1);

Um1 = U1(:, 1 : 2 * m);
Sm1 = S1(1 : 2 * m, 1 : 2 * m);
Vm1 = V1(:, 1 : 2 * m);

% Half power and inverse matrices. Transform to avoid computing over the
% entire S matrix. Only saves time for large S (assumed).
singularValues  = diag(Sm1);
sHalf = diag(sqrt(singularValues));
sHalfInv = diag(1 ./ sqrt(singularValues)); 

% Eq. (22-24).
A1 = sHalfInv * Vm1' * H0 * Um1 * sHalfInv;

[psi, lambda] = eig(A1);

% Discrete to continuous time conversion
continuousEigenValues = log(lambda) ./ dt;

naturalFrequencies = abs(continuousEigenValues);

% Extract frequencies from the eigenvalues.
shiftedFrequencies = diag(naturalFrequencies ./ (2*pi));

% Extract the decay rate.
shiftedDecayRates = -real(continuousEigenValues);

% Extract damping factor from the eigenvalues.
shiftedDampingFactors = -diag(real(continuousEigenValues) ./ ...
    naturalFrequencies);

figure
plot(f_true, zeta_true, 'rsq', 'LineWidth', 1.4, 'DisplayName', ...
    'True eigenvalues', 'MarkerSize', 10)
hold on
plot(f_era, zeta_era, 'ko', 'LineWidth', 1.4, ...
    'DisplayName', 'Standard frequencies')
plot(shiftedFrequencies, -shiftedDampingFactors, 'gx', 'LineWidth', 1.4, ...
    'DisplayName', 'Shifted Frequencies')
ylabel('Damping factor')
xlabel('Frequency')
legend


%% Subfunctions

function msi = modeSimilarityIndex(freq1, freq2, df1, df2)

eigenfrequencyWeight = 0.5;

dampingFactorWeight = 0.5;

freqRelError = 0.01;
dampingRelError = 0.05;

for k = 1 : numel(freq1)

    P1 = (eigenfrequencyWeight / freqRelError) * ...
        (abs(freq1(k) - freq2(k)) / max(freq1(k), freq2(k)));

    P2 = (dampingFactorWeight / dampingRelError) * ...
        (abs(df1(k) - df2(k)) / max(df1(k), df2(k)));

    msi(k) = P1 + P2;

end

end