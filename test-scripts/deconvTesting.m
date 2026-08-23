clear
close all
clc
%% Load files.
filepath = "C:\Users\ademp\OneDrive - Queen's University Belfast\" + ...
    "Maarten Van Walstijn's files - Alastair\AD FILES\Measurements\" + ...
    "hand-measurements_05-08";

dataInfo = dir(filepath);
filesToInspect = {dataInfo(3 : end).name}';
filePaths = fullfile(filepath, filesToInspect);
nFiles = length(filePaths);

% LS convolution loading factor.
lambda = 1e-4;

fig = figure;
tld = tiledlayout(fig, 2, 3);
title(tld, "LS Loading factor \lambda: " + string(lambda))

fig2 = figure;
tld2 = tiledlayout(fig2, 2, 3);
title(tld2, "LS Loading factor \lambda: " + string(lambda))

fig3 = figure;
tld3 = tiledlayout(fig3, 2, 3);

impulseStartIndex = [220254;
                     172992;
                     256364;
                     1;
                     187798];

impulseEndIndex = impulseStartIndex + 44100;

opts = eraConfig();
opts.returnDebug = true;
opts.zeroFrequencyFilter = true;
svdTolerance = 50;

for iFile = 1 : nFiles
    
    thisFilepath = filePaths(iFile);
    [audioData, fs] = audioread(thisFilepath);
    thisImpulseStartIndex = impulseStartIndex(iFile);
    thisImpulseEndIndex = impulseEndIndex(iFile);

    celloOutput = audioData(thisImpulseStartIndex : thisImpulseEndIndex, ...
        1);
    hammerInput = audioData(thisImpulseStartIndex : thisImpulseEndIndex, ...
        2);

    celloOutput = celloOutput - mean(celloOutput);
    hammerInput = hammerInput - mean(hammerInput);

    celloOutput = resample(celloOutput, 1, 4);
    hammerInput = resample(hammerInput, 1, 4);
    fs = 44100 / 4;
    
    [celloSpectrum, freqAxis] = singleSidedFft(celloOutput, fs);
    [hammerSpectrum, ~] = singleSidedFft(hammerInput, fs);

    celloFft = fft(celloOutput);
    hammerFft = fft(hammerInput);

    N        = numel(hammerFft);
    binLo    = round(1000 * N / fs);
    binHi    = round(5000 * N / fs);
    noiseFloor = mean(abs(hammerFft(binLo:binHi)).^2);
    % lambda   = 0.001 * noiseFloor;

    % impulseResponse = real(ifft(frf));

    % Least squares
    nTaps = 1000;
    X = convmtx(hammerInput(:), nTaps);   % convolution matrix
    X = X(1 : length(celloOutput), :);    % trim to output length
    hLS = (X' * X + lambda * eye(nTaps)) \ (X' * celloOutput(:));

    
    [frfAxis, frf] = singleSidedFft(hLS, fs);

    [eraSignal, frequencies, amps, phi, decayFactors, modes, debugEra] = ...
        eigensystemRealisation(hLS, svdTolerance, fs, opts);

    [eraAxis, eraFRF] = singleSidedFft(eraSignal, fs);

    residues = debugEra.Modal.residues;
    poles = debugEra.Modal.poles;
    inputSignal = zeros(1, length(celloOutput));
    inputSignal(1) = 1;

    [B, A] = buildFilterCoefficients(poles, residues, eraSignal(1));
    parallelFilterSignal = parallelFilter(B, A, inputSignal);

    til = nexttile(tld);
    hold(til, "on")
    plot(til, celloOutput, 'LineWidth', 1.5, 'DisplayName', 'Cello admittance')
    plot(til, hammerInput, 'LineWidth', 1.5, 'DisplayName', 'Hammer Impulse')
    plot(til, hLS, '--', 'LineWidth', 1.4, 'DisplayName', 'Impulse response')
    % plot(til, eraSignal, 'LineWidth', 1.5, 'DisplayName', "ERA signal")
    plot(til, parallelFilterSignal, '--', 'LineWidth', 1.5, DisplayName='Parallel Filter signal')
    legend

    til2 = nexttile(tld2);
    hold(til2, "on")
    plot(til2, freqAxis, 20*log10(celloSpectrum), 'LineWidth', 1.5, 'DisplayName', 'Cello admittance')
    plot(til2, freqAxis, 20*log10(hammerSpectrum), 'LineWidth', 1.5, 'DisplayName', 'Hammer Impulse')
    % yline(noiseFloor, 'r--', 'LineWidth', 1.5)
    yline(20*log10(lambda), 'b--', 'LineWidth', 1.5)
    plot(til2, freqAxis, 20*log10(abs(frf(1 : 5514)) ./ length(freqAxis)), '--', 'LineWidth', 1.5, 'DisplayName', 'Hammer Impulse')
    xlim([0 1000])
    legend

    til3 = nexttile(tld3);


end