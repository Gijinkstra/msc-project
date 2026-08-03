%% Manage interpreters
set(groot, 'defaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultLegendInterpreter', 'latex');

%% Cello dominant modes
close all
clear
clc

% These values have been selected from Thomas D. Rossings "The Science of
% String Instruments".
frequencies = [57, 102, 144, 170, 195, 203, 219, 277, 302];
nSignals = length(frequencies);

fs = 44.1E3;
totalTime = 0.5;
timeStep = 1 / fs;
nSamples = totalTime / timeStep - 1;
freqStep = fs / nSamples;
frequencyAxis = (0 : (nSamples / 2)) * freqStep;

SNR = 5; % [dB]

[modes, signal, timeVector, debugImpulse] = ...
    createImpulseResponse(nSignals, nSamples, fs, ...
    Frequency=frequencies, SNR=SNR, Debug=true);

frequencyResponse = fft(signal);
P2 = abs(frequencyResponse / nSamples);
P1 = P2(1 : nSamples / 2+1);
P1(2 : end - 1) = 2*P1(2 : end - 1);

figure
tiledlayout(3, 1)

nexttile
plot(timeVector, modes, 'LineWidth', 1.1)
xlabel('Time (s)')
ylabel('Amplitude')
title('Modes')
xlim([0 0.1])

nexttile
plot(timeVector, signal, 'LineWidth', 1.5)
xlabel('Time (s)')
ylabel('Amplitude')
title('Simulated cello impulse response')
xlim([0 0.1])

nexttile
plot(frequencyAxis, mag2db(P1), 'LineWidth', 1.5);
xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Frequency Response Function')
xlim([0 500])

opts = eraConfig();
opts.returnDebug = true;
opts.poleScaling = true;
svdTolerance = 12;

endIndex = 6000;

[eraImpulse, f_era, zeta_era, amp_era, phase_era, modes_era, debug] = ...
    eigensystemRealisation(signal(1 : endIndex), svdTolerance, fs, opts);

%%

sv = debug.SVD.singularValues;
normSv = sv ./ sv(1);
diffSv = diff(normSv);
numel(diffSv(diffSv < -0.01))

medianSv = median(normSv);
cutoff = 2.858 * medianSv;
retainedSv = normSv(normSv > cutoff);

svKurt = (normSv(1 : end - 1) - normSv(2 : end)) ./ normSv(2 : end);



% Convert to sinusoidal phase.
phase_era = phase_era + pi/2;

[eraModes, eraReconstructedImpulse, timeVector] = ...
    createImpulseResponse(length(f_era), nSamples, fs, Frequency=f_era, ...
    Amplitude=amp_era, Phase=phase_era, Alpha=debug.Modal.decayRates);

figure
tiledlayout(2, 1)
nexttile
plot(timeVector(1 : endIndex), signal(1 : endIndex), 'LineWidth', 1.5, 'DisplayName', 'True signal')
hold on
plot(timeVector, eraReconstructedImpulse, 'LineWidth', 1.5, 'DisplayName', 'ERA reconstructed signal')
xlabel('Time (s)')
ylabel('Amplitude')
title('Simulated cello impulse response')
xlim([0 endIndex * timeStep])
legend

nexttile
hold on
plot(normSv, 'k.', 'DisplayName', 'Singular values')
plot(retainedSv, 'or', 'DisplayName', 'Singular values > 3 * Median Absolute Deviation')
ylabel('Normalised singular values')
xlabel('Singular value index')
title("Normalised singular values of cello impulse response - SNR: " + string(SNR))
xlim([0 100])
legend
text(0.2, 0.9, "Number of retained MAD modes: " + string(ceil(length(retainedSv) / 2)), 'Units', 'normalized')

figure
plot(diff(normSv), 'ko', 'LineWidth', 1.3)

f1 = figure;
t1 = tiledlayout(f1);
f2 = figure;
t2 = tiledlayout(f2);
f3 = figure;
t3 = tiledlayout(f3);

U = debug.SVD.Um;
V = debug.SVD.Vm;

[r, c] = size(U);

for i = 1 : c

    tt = nexttile(t1);
    thisSig = U(:, i);
    plot(tt, thisSig)

    smoothness(i) = sum(abs(diff(thisSig))) / ...
        (r - 1) * (max(thisSig) - min(thisSig));
    title(tt, "Smoothness: " + string(smoothness(i)));

    ttt = nexttile(t2);
    thisSigV = V(:, i);
    plot(ttt, thisSigV);

    tttt = nexttile(t3);
    nsamp = r;
    freqAx = (0 : nsamp / 2) * freqStep;
    frequencyResponse = fft(thisSig);
    P2 = abs(frequencyResponse / nsamp);
    P1 = P2(1 : nsamp / 2+1);
    P1(2 : end - 1) = 2*P1(2 : end - 1);
    plot(tttt, freqAx, P1);
    xlim([0 500])

    if i > 20
        break
    end

end

xlabel(t2, 'Frequency (Hz)')
ylabel(t2, 'Magniutde')

figure
bar(smoothness)
