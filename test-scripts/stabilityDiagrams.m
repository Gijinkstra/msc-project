close all
clear 
clc

%% Setup
dataPath = "Project\data\bridge_LDV";
dataInfo = dir(dataPath);
filesToInspect = {dataInfo(3 : end).name}';
filePaths = fullfile(dataPath, filesToInspect);
nFiles = length(filePaths);
% nFiles = 1;
allEraFreqs = cell(1, nFiles);
allEraPhase = cell(1, nFiles);
allEraAlpha = cell(1, nFiles);
allEraAmplitude = cell(1, nFiles);

% Define signal impulse start samples and noise periods for wiener filter.
impulseStartSample = [2482
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
svdTolerance = 2 : 2 : 60;
impulseErrors = 1 : length(svdTolerance);
frfError = 1 : length(svdTolerance);

fig = figure;
til = tiledlayout(ceil(nFiles / 4), ceil(nFiles / 4));
ylabel(til, 'Mode order')
xlabel(til, 'Frequency (Hz)')

fig2 = figure;
til2 = tiledlayout(ceil(nFiles / 4), ceil(nFiles / 4));
xlabel(til2, 'Mode order')
ylabel(til2, 'NMSE')

for iFile = 1 : nFiles

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

    t1 = nexttile(til);
    hold(t1, "on")

    for k = 1 : length(svdTolerance)

        [eraImpulse, freqEra, zetaEra, ampEra, phaseEra, modesEra, debug] = ...
            eigensystemRealisation(thisFilteredSignal, svdTolerance(k), newFs, opts);

        scatter(t1, freqEra, svdTolerance(k), '.r')

        impulseErrors(k) = debug.Reconstruction.impulseNMSE;
        frfError(k) = nmseFreq(thisFilteredSignal, eraImpulse);

    end

    t2 = nexttile(til2);
    plot(t2, 1 : 30, impulseErrors, 'ko')

end

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

function freqError = nmseFreq(originalSignal, eraSignal)

fftSig = fft(originalSignal);
fftEra = fft(eraSignal);
fftSig = fftSig(:);
fftEra = fftEra(:);

freqError = sum((abs(fftSig) - abs(fftEra)).^2) / sum(abs(fftSig).^2);

end