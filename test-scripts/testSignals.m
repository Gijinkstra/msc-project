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

freqs = [300, 300, 310, 500, 600, 1000, 2390];

amplitudes = [1.5, 0.3, 0.75, 1, 1, 1, 1];

alpha = 50 : 10 : 110;

phases = pi * (0.1 : 0.1 : 0.7);

M = length(freqs);
sigs = zeros(M,Ns);

for m = 1 : M
    sigs(m, :) = amplitudes(m) * exp(-alpha(m) * t) .* ...
        sin(2 * pi * freqs(m) * t + (pi / 180) * phases(m));
end

sig = sum(sigs);

opts = eraConfig();
opts.returnDebug = true;
svdTolerance = 0.2;

% Try changing the number of modes to see how the system reconstructs a
% subspace of the original signal.
[outputSignal, frequencies, decayFactors, modes, debug] = ...
    eigensystemRealisation(sig, svdTolerance, Fs, opts);

f = figure;
til = tiledlayout(3, 3);

tf = nexttile(1);
title('7 Individual impulse responses')
hold(tf, "on")
xlabel('Time (ms)')
ylabel('Amplitude')

tf2 = nexttile(2);
title('Frequency response of individual impulse responses')
xlabel('Frequency (kHz)')
ylabel('Amplitude')
hold(tf2, "on")

tf3 = nexttile(4);
title("Reconstructed impulse responses with a singular value tolerance " + ...
    "of "  +  string(svdTolerance))
hold(tf3, "on")
xlabel('Time (ms)')
ylabel('Amplitude')

tf4 = nexttile(5);
title('Frequency response of reconstructed individual impulse responses')
xlabel('Frequency (kHz)')
ylabel('Amplitude')
hold(tf4, "on")

for m = 1 : M
    plot(tf, t/1000, sigs(m, :), 'LineWidth', 1.1, 'DisplayName', ...
        "Frequency: " + string(freqs(m)))


    frequencyResponse = fft(sigs(m, :));
    P2 = abs(frequencyResponse / Ns);
    P1 = P2(1:Ns/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    plot(tf2, frequencyAxis / 1000, P1, 'LineWidth', 1.1, 'DisplayName', ...
        "Frequency: " + string(freqs(m)))
end

for m = 1 : numel(frequencies)
    plot(tf3, t/1000, debug.Reconstruction.modalImpulses(m, :), 'LineWidth', 1.1, 'DisplayName', ...
        "Frequency: " + string(frequencies(m)))

    frequencyResponse = fft(debug.Reconstruction.modalImpulses(m, :));
    P2 = abs(frequencyResponse / Ns);
    P1 = P2(1:Ns/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    plot(tf4, frequencyAxis / 1000, P1, 'LineWidth', 1.1, 'DisplayName', ...
        "Frequency: " + string(frequencies(m)))

end

legend(tf)
legend(tf2)
legend(tf3)
legend(tf4)

nexttile(7)
title(['Composite impulse response and reconstructed impulse response' ...
    ' using ERA'])
hold on
plot(t / 1000, sig, 'k', 'LineWidth', 1.5, 'DisplayName', 'Original signal')
plot(t / 1000, outputSignal, '--r', 'LineWidth', 1.5, 'DisplayName', ...
    'Reconstructed signal')
legend
xlabel('Time (ms)')
ylabel('Amplitude')
reconError = sprintf("NMSE: %.4d", debug.Reconstruction.impulseNMSE);
text(0.7, 0.3, reconError, 'Units', 'normalized')

nexttile(8)
hold on
title('Associated frequency response and reconstructed frequency response')

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
xlabel('Frequency (kHz)')
ylabel('Amplitude')
title('One sided frequency response of impulse function')
grid;

nexttile([3, 1])
semilogy(freqs, alpha, 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, ...
    'DisplayName', ...
    'Target decay factors per frequency')
hold on
semilogy(frequencies, debug.Modal.decayRates, 'rx', 'MarkerSize', 10, ...
    'LineWidth', 1.5, 'DisplayName', 'Recalculated decay factors')
xlabel('Frequency (Hz)')
ylabel('Decay rate: \alpha: [s^{-1}]')
xlim([0 2.5E3])
ylim([ 0 150])
legend
title(['Comparison of decay rate for simulated signals ' ...
    'vs. ERA reconstruction'])

%% Case 2 - Signals with noise.
clear;
clc;

Fs = 10000;
dt = 1/Fs;
Ns = 1000;
t = (0:(Ns-1))*dt;
frequencyAxis = (0 : (Ns / 2)) * (Fs / Ns);

SNR = 10;

freqs = [500, 1000, 1500];

amplitudes = [1.5, 0.3, 0.75];

alpha = 50 : 30 : 110;

phases = pi * (0.1 : 0.5 : 1.1);

M = length(freqs);
sigs = zeros(M,Ns);

for m = 1 : M
    sigs(m, :) = amplitudes(m) * exp(-alpha(m) * t) .* ...
        sin(2 * pi * freqs(m) * t + (pi / 180) * phases(m));
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

% Try changing the number of modes to see how the system reconstructs a
% subspace of the original signal.
[outputSignal, frequencies, dampingFactors, modes, debug] = ...
    eigensystemRealisation(sig, 3, Fs, opts);

debug.impulseNMSE = sum((sig - outputSignal).^2) ./ sum(sig.^2);

f = figure;
til = tiledlayout(3, 3);

tf = nexttile(1);
title(string(M) +  "Individual impulse responses")
hold(tf, "on")
xlabel('Time (ms)')
ylabel('Amplitude')

tf2 = nexttile(2);
title('Frequency response of individual impulse responses')
xlabel('Frequency (kHz)')
ylabel('Amplitude')
hold(tf2, "on")

tf3 = nexttile(4);
title("Reconstructed impulse responses with a singular value tolerance " + ...
    "of "  +  string(opts.svdTolerance))
hold(tf3, "on")
xlabel('Time (ms)')
ylabel('Amplitude')

tf4 = nexttile(5);
title('Frequency response of reconstructed individual impulse responses')
xlabel('Frequency (kHz)')
ylabel('Amplitude')
hold(tf4, "on")

for m = 1 : M
    plot(tf, t/1000, sigs(m, :), 'LineWidth', 1.1, 'DisplayName', ...
        "Frequency: " + string(freqs(m)))


    frequencyResponse = fft(sigs(m, :));
    P2 = abs(frequencyResponse / Ns);
    P1 = P2(1:Ns/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    plot(tf2, frequencyAxis / 1000, P1, 'LineWidth', 1.1, 'DisplayName', ...
        "Frequency: " + string(freqs(m)))
end

for m = 1 : numel(frequencies)
    plot(tf3, t/1000, debug.Reconstruction.modalImpulses(m, :), 'LineWidth', 1.1, 'DisplayName', ...
        "Frequency: " + string(frequencies(m)))

    frequencyResponse = fft(debug.Reconstruction.modalImpulses(m, :));
    P2 = abs(frequencyResponse / Ns);
    P1 = P2(1:Ns/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    plot(tf4, frequencyAxis / 1000, P1, 'LineWidth', 1.1, 'DisplayName', ...
        "Frequency: " + string(frequencies(m)))

end

legend(tf)
legend(tf2)
legend(tf3)
legend(tf4)

nexttile(7)
title(['Composite impulse response and reconstructed impulse response' ...
    ' using ERA'])
hold on
plot(t / 1000, sig, 'k', 'LineWidth', 1.5, 'DisplayName', 'Original signal')
plot(t / 1000, outputSignal, '--r', 'LineWidth', 1.5, 'DisplayName', ...
    'Reconstructed signal')
legend
xlabel('Time (ms)')
ylabel('Amplitude')
reconError = sprintf("NMSE: %.4d", debug.impulseNMSE);
text(0.7, 0.3, reconError, 'Units', 'normalized')

nexttile(8)
hold on
title('Associated frequency response and reconstructed frequency response')

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
xlabel('Frequency (kHz)')
ylabel('Amplitude')
title('One sided frequency response of impulse function')
grid;

nexttile([3, 1])
semilogy(freqs, alpha, 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, ...
    'DisplayName', ...
    'Target decay factors per frequency')
hold on
semilogy(frequencies, debug.Modal.decayRates, 'rx', 'MarkerSize', 10, ...
    'LineWidth', 1.5, 'DisplayName', 'Recalculated decay factors')
xlabel('Frequency (Hz)')
ylabel('Decay rate: \alpha: [s^{-1}]')
xlim([0 2.5E3])
ylim([ 0 150])
legend
title(['Comparison of decay rate for simulated signals ' ...
    'vs. ERA reconstruction'])
%% Noisy signals, more modes
clear;
clc;

SNR = 30;

Fs = 10000;
dt = 1/Fs;
Ns = 1000;
t = (0:(Ns-1))*dt;
frequencyAxis = (0 : (Ns / 2)) * (Fs / Ns);

fmin = 20;
fmax = 500;
M = 20;    
freqs = fmin + (fmax - fmin)*rand(1,M);
alpha = 10 + 10*rand(1,M);
amplitudes = rand(1,M);
phases = -180 + 360*rand(1,M);

%%%%%%%% BUILD TARGET SIGNAL %%%%
M = length(freqs);
sigs = zeros(M,Ns);
for m=1:M
    sigs(m, :) = amplitudes(m)*exp(-alpha(m)*t).*sin(2*pi*freqs(m)*t + (pi/180)*phases(m));
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

% Try changing the number of modes to see how the system reconstructs a
% subspace of the original signal.
[outputSignal, frequencies, dampingFactors, modes, debug] = ...
    eigensystemRealisation(sig, 0.1, Fs, opts);

debug.impulseNMSE =  sum(((sig - outputSignal).^2)) ./ ...
    sum(sig.^2);

f = figure;
til = tiledlayout(3, 3);

tf = nexttile(1);
title(string(M) +  "Individual impulse responses")
hold(tf, "on")
xlabel('Time (ms)')
ylabel('Amplitude')

tf2 = nexttile(2);
title('Frequency response of individual impulse responses')
xlabel('Frequency (kHz)')
ylabel('Amplitude')
hold(tf2, "on")

tf3 = nexttile(4);
title("Reconstructed impulse responses with a singular value tolerance " + ...
    "of "  +  string(opts.svdTolerance))
hold(tf3, "on")
xlabel('Time (ms)')
ylabel('Amplitude')

tf4 = nexttile(5);
title('Frequency response of reconstructed individual impulse responses')
xlabel('Frequency (kHz)')
ylabel('Amplitude')
hold(tf4, "on")

for m = 1 : M
    plot(tf, t/1000, sigs(m, :), 'LineWidth', 1.1, 'DisplayName', ...
        "Frequency: " + string(freqs(m)))


    frequencyResponse = fft(sigs(m, :));
    P2 = abs(frequencyResponse / Ns);
    P1 = P2(1:Ns/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    plot(tf2, frequencyAxis / 1000, P1, 'LineWidth', 1.1, 'DisplayName', ...
        "Frequency: " + string(freqs(m)))
end

for m = 1 : numel(frequencies)
    plot(tf3, t/1000, debug.Reconstruction.modalImpulses(m, :), 'LineWidth', 1.1, 'DisplayName', ...
        "Frequency: " + string(frequencies(m)))

    frequencyResponse = fft(debug.Reconstruction.modalImpulses(m, :));
    P2 = abs(frequencyResponse / Ns);
    P1 = P2(1:Ns/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    plot(tf4, frequencyAxis / 1000, P1, 'LineWidth', 1.1, 'DisplayName', ...
        "Frequency: " + string(frequencies(m)))

end

legend(tf)
legend(tf2)
legend(tf3)
legend(tf4)

nexttile(7)
title(['Composite impulse response and reconstructed impulse response' ...
    ' using ERA'])
hold on
plot(t / 1000, sig, 'k', 'LineWidth', 1.5, 'DisplayName', 'Original signal')
plot(t / 1000, outputSignal, '--r', 'LineWidth', 1.5, 'DisplayName', ...
    'Reconstructed signal')
legend
xlabel('Time (ms)')
ylabel('Amplitude')
reconError = sprintf("NMSE: %.4d", debug.impulseNMSE);
text(0.7, 0.3, reconError, 'Units', 'normalized')

nexttile(8)
hold on
title('Associated frequency response and reconstructed frequency response')

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
xlabel('Frequency (kHz)')
ylabel('Amplitude')
title('One sided frequency response of impulse function')
grid;

nexttile([3, 1])
semilogy(freqs, alpha, 'ko', 'MarkerSize', 10, 'LineWidth', 1.5, ...
    'DisplayName', ...
    'Target decay factors per frequency')
hold on
semilogy(frequencies, debug.Modal.decayRates, 'rx', 'MarkerSize', 10, ...
    'LineWidth', 1.5, 'DisplayName', 'Recalculated decay factors')
xlabel('Frequency (Hz)')
ylabel('Decay rate: \alpha: [s^{-1}]')
xlim([0 2.5E3])
ylim([ 0 150])
legend
title(['Comparison of decay rate for simulated signals ' ...
    'vs. ERA reconstruction'])