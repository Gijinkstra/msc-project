close all
clear 
clc

%% Setup
dataPath = "Project\data\bridge_LDV";
dataInfo = dir(dataPath);
filesToInspect = {dataInfo(3 : end).name}';
filePaths = fullfile(dataPath, filesToInspect);
nFiles = length(filePaths);
allEraFreqs = cell(1, nFiles);
allEraPhase = cell(1, nFiles);
allEraAlpha = cell(1, nFiles);
allEraAmplitude = cell(1, nFiles);

% Define signal impulse start samples and noise periods for wiener filter.
impulseStartSample = [2481
                      2446
                      2476
                      2460
                      2518
                      2538
                      2503
                      2454
                      2501
                      2489];

noiseEndIndex = 2000;

originalFs = 100000;
resampleFactor = 10;
newFs = originalFs / resampleFactor;
nTaps = 500;

% ERA parameters;
opts = eraConfig();
opts.returnDebug = true;
% opts.residueScaling = true;
svdTolerance = 50;

% Modal resynthesis parameters
SNR = 0;

for iFile = 1 : nFiles

    f1 = figure;
    tl1 = tiledlayout(f1, 3, 2);

    % Plot 1 - Original impulse and modal resynthesis.
    a1 = nexttile(tl1, [1, 2]);
    hold(a1, "on");
    xlabel(a1, 'Time (ms)')
    ylabel(a1, 'Amplitude')
    title(a1, "Comparison of original impulse response and modal " ...
         + "resynthesis via ERA, svdTolerance: " + string(svdTolerance))
    legend(a1, "AutoUpdate", "on")
    xlim([0 0.2])

    % Plot 2 - Normalised singular values with MAD outliers.
    a2 = nexttile;
    hold(a2, "on");
    ylabel(a2, 'Normalised singular values')
    xlabel(a2, 'Singular value index')
    title(a2, "Normalised singular values of cello impulse response - SNR: " ...
        + string(SNR))
    xlim(a2, [0 100])
    legend(a2, "AutoUpdate", "on")

    % Plot 3 - Derivative of singular values.
    a3 = nexttile;
    hold(a3, "on");
    ylabel(a3, '-\Delta\sigma')
    xlabel(a3, 'Singular Value index')
    title(a3, 'Difference between contiguous singular values')
    xlim(a3, [0 100])

    % Plot 4 - Frequency Spectrum of original signal and reconstruction.
    a4 = nexttile;
    hold(a4, "on");
    ylabel('Magnitude')
    xlabel('Frequency (Hz)')
    title(['Comparison of frequency response of original impulse and ' ...
        'modal resynthesis signal'])
    grid(a4, "on")
    legend(a4, "AutoUpdate", 'on')
    xlim([0 1000])

    % Plot 5 - Left singular value smoothness.
    a5 = nexttile;
    hold(a5, "on");
    xlabel(a5, 'Mode Index')
    ylabel(a5, 'Left singular vector smoothness value')
    title('Smoothness factor of left singular vectors')
    
    % Define file parameters.
    thisFilePath = filePaths(iFile);
    thisFileData = readmatrix(thisFilePath);
    nSignals = size(thisFileData, 2);
    
    % The signal parameters for this file (Signal 5 only is defined as it
    % is assumed to be the most similar to the setup that will be worked
    % on).
    thisStartSample = impulseStartSample(iFile);
    thisSignal = thisFileData(:, 5);
    thisSignal = thisSignal - mean(thisSignal);
    thisSignal = resample(thisSignal, 1, resampleFactor);

    % Split up noise and filter.
    thisNoiseSignal = thisSignal(1 : noiseEndIndex);
    % Truncate the signal to only the impulse response data.
    thisSignal = thisSignal(thisStartSample : end);
    nSamples = length(thisSignal);
    thisFilteredSignal = wienerFilter(thisSignal, thisNoiseSignal, nTaps);

    % Perform ERA.
    [eraImpulse, freqEra, zetaEra, ampEra, phaseEra, modesEra, debug] = ...
        eigensystemRealisation(thisFilteredSignal, svdTolerance, newFs, opts);

    % Convert to sinusoidal phase.
    phaseEra = phaseEra + pi/2;

    % Modal reconstruction with ERA defined frequencies.
    [eraModes, eraReconstructedImpulse, timeVector] = ...
        createImpulseResponse(length(freqEra), nSamples, newFs, ...
        Frequency=freqEra, Amplitude=ampEra, Phase=phaseEra, ...
        Alpha=debug.Modal.decayRates, SNR=0);

    % Frequency spectrum.
    [impulseSpectrum, frequencyAxis] = ...
        singleSidedFft(thisFilteredSignal, newFs);
    [reconstructedSpectrum, ~] = ...
        singleSidedFft(eraReconstructedImpulse, newFs);
    [eraSignalSpectrum, ~] = ...
        singleSidedFft(eraImpulse, newFs);
    
    % Singular value processing.
    singularValues = debug.SVD.singularValues;
    normSingularValues = singularValues ./ singularValues(1);
    diffSingularValues = diff(normSingularValues(1 : 2 : end));
    numel(diffSingularValues(diffSingularValues < -0.01))
    medianSv = median(normSingularValues);
    cutoff = 2.858 * medianSv;
    retainedSv = normSingularValues(normSingularValues > cutoff);
    svKurt = (normSingularValues(1 : end - 1) - ...
        normSingularValues(2 : end)) ./ normSingularValues(2 : end);

    % Plot 1 - Original impulse and modal resynthesis.
    plot(a1, timeVector, thisFilteredSignal, 'k', 'LineWidth', 1.5, ...
        'DisplayName', 'True signal')
    plot(a1, timeVector, eraReconstructedImpulse, 'r', 'LineWidth', 1.3, ...
        'DisplayName', 'Resynthesized signal')
    plot(a1, timeVector, eraImpulse, 'b', 'LineWidth', 2, ...
        'DisplayName', 'ERA reconstructed impulse')

    % Plot 2 - Normalised singular values with MAD outliers.
    plot(a2, normSingularValues, 'k.', 'DisplayName', 'Singular values', ...
        'MarkerSize', 15)
    plot(a2, normSingularValues(1 : 2*svdTolerance), 'or', 'DisplayName', ...
        'Retained singular values via manual truncation', 'LineWidth', ...
        1.5, 'MarkerSize', 6)
    plot(a2, retainedSv, 'y.', 'DisplayName', ...
        'Singular values > 3 * singular values MAD', 'LineWidth', 1.5, ...
        'MarkerSize', 8)
    text(a2, 0.2, 0.9, "Number of retained MAD modes: " + ...
        string(ceil(length(retainedSv) / 2)), 'Units', 'normalized')
    
    % Plot 3 - Derivative of singular values.
    plot(a3, -diff(normSingularValues), 'ko', 'LineWidth', 1.3)

    % Plot 4 - Frequency Spectrum.
    plot(a4, frequencyAxis, 20*log10(impulseSpectrum), 'k',  'LineWidth', 2, ...
    'DisplayName', 'True signal')
    plot(a4, frequencyAxis, 20*log10(reconstructedSpectrum), 'r', 'LineWidth', 2, ...
        'DisplayName', 'Resynthesized signal')
    plot(a4, frequencyAxis, 20*log10(eraSignalSpectrum), 'b', 'LineWidth', 2, ...
        'DisplayName', 'ERA constructed signal')

    allEraFreqs{iFile} = freqEra;
    allEraPhase{iFile} = phaseEra;
    allEraAlpha{iFile} = debug.Modal.decayRates;
    allEraAmplitude{iFile} = ampEra;

end

%% Subfunctions
function filteredSignal = wienerFilter(signal, noiseSignal, nTaps)

ryy = xcorr(signal, nTaps - 1);
rnn = xcorr(noiseSignal, nTaps - 1);
ryy = ryy(nTaps : end);
rnn = rnn(nTaps : end);

Ryy = toeplitz(ryy);
rxy = ryy - rnn;
rxy = rxy(:);
B = Ryy \ rxy;

B = B(:);

filteredSignal = fftfilt(B, signal);

end