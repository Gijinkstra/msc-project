close all
clear
clc
%%
fs = 1000;
time = 2;
nSamples = (time * fs);
freqs = [20 5 2 60 80];
nSignals = length(freqs);
alpha = 5 * ones(1, nSignals);
inputSignal = zeros(1, nSamples);
inputSignal(1) = 1;

opts = eraConfig();
opts.returnDebug = true;
% opts.residueScaling = true;
SNR = 10;

[allSignals, outputSignal, timeVector, debug] = ...
    createImpulseResponse(nSignals, nSamples, fs, Alpha=alpha, SNR=SNR);

svdTolerance = 5;

[eraSignal, frequencies, amps, phi, decayFactors, modes, debugEra] = ...
    eigensystemRealisation(outputSignal, svdTolerance, fs, opts);

figure
tiledlayout
nexttile
plot(timeVector, outputSignal, 'LineWidth', 1.5, 'DisplayName', 'Original Signal')
hold on
plot(timeVector, eraSignal, '--k', 'LineWidth', 1.5, 'DisplayName', 'ERA Signal')
legend("AutoUpdate", "on")

residues = debugEra.Modal.residues;
poles = debugEra.Modal.poles;

[B, A] = buildFilterCoefficients(poles, residues, outputSignal(1));
parallelFilterSignal = parallelFilter(B, A, inputSignal);

nexttile
plot(timeVector, eraSignal, '--k', 'LineWidth', 2, 'DisplayName', 'ERA signal')
hold on
plot(timeVector, parallelFilterSignal, '.-r', 'LineWidth', 1.3, 'DisplayName', 'Parallel filter function signal')
legend


% Filter the poles and residues
filteredPoles = diag(debugEra.Modal.At);
filteredResidues = debugEra.Modal.Ct' .* debugEra.Modal.Bt;

[BT, AT] = buildFilterCoefficients(filteredPoles, filteredResidues, outputSignal(1));
parallelFilterSignalFiltered = parallelFilter(BT, AT, inputSignal);

nexttile
plot(timeVector, eraSignal, '--k', 'LineWidth', 1.5, 'DisplayName', 'ERA signal')
hold on
plot(timeVector, parallelFilterSignalFiltered, '.-r', 'LineWidth', 1.5, 'DisplayName', 'Parallel filter function signal (filtered poles)')
legend


%% Setup
% dataPath = '/Users/ali/Dev/msc-project/data/bridge_LDV'; % MAC filepath.
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
                      2447
                      2477
                      2460
                      2518
                      2538
                      2503
                      2454
                      2501
                      2489];
impulseEndIndex = impulseStartSample + 3000;

noiseEndIndex = 2000;

originalFs = 100000;
resampleFactor = 10;
newFs = originalFs / resampleFactor;
nTaps = 500;

% ERA parameters;
opts = eraConfig();
opts.returnDebug = true;
opts.zeroFrequencyFilter = false;
svdTolerance = 55;
impulseErrors = 1 : length(svdTolerance);
frfError = 1 : length(svdTolerance);

fig = figure;
til = tiledlayout(3, 3);
ylabel(til, 'Amplitude')
xlabel(til, 'Time (s)')

fig2 = figure;
til2 = tiledlayout(3, 3);
xlabel(til2, 'Frequency (Hz)')
ylabel(til2, 'Magnitude (dB)')

for iFile = 1 : nFiles
 
    % Define file parameters.
    thisFilePath = filePaths(iFile);
    thisFileData = readmatrix(thisFilePath{1});
    nSignals = size(thisFileData, 2);

    % The signal parameters for this file (Signal 5 only is defined as it
    % is assumed to be the most similar to the setup that will be worked
    % on).
    thisStartSample = impulseStartSample(iFile);
    thisEndSample = impulseEndIndex(iFile);
    thisSignal = thisFileData(:, 5);
    thisSignal = thisSignal - mean(thisSignal);
    thisSignal = resample(thisSignal, 1, resampleFactor);

    % Split up noise and filter.
    thisNoiseSignal = thisSignal(1 : noiseEndIndex);
    % Truncate the signal to only the impulse response data.
    thisSignal = thisSignal(thisStartSample : thisEndSample);
    nSamples = length(thisSignal);
    thisFilteredSignal = wienerFilter(thisSignal, thisNoiseSignal, nTaps);
    timeVector = (0 : nSamples - 1) / newFs;

    % MPR
    [~, thisFilteredSignal] = rceps(thisFilteredSignal);

    % Build the impulse signal.
    inputSignal = zeros(1, nSamples);
    inputSignal(1) = 1;

    % ERA.
    t1 = tic;
    [eraSignal, frequencies, dampingFactors, amplitudes, ...
        phases, modes, debug] = eigensystemRealisation(thisFilteredSignal, ...
        svdTolerance, newFs, opts);
    t2 = toc;
    eraDur = t2 - t1;
    sprintf("ERA time duration: " + string(eraDur))

    poles = debug.Modal.poles;
    residues = debug.Modal.residues;

    % Build the filter coefficients from ERA.
    [B, A] = buildFilterCoefficients(poles, ...
        residues, thisFilteredSignal(1));

    % Parallel filter the impulse signal.
    thisParallelFilteredSignal = parallelFilter(B, A, inputSignal);

    [filteredSignalFft, freqAxis] = singleSidedFft(thisFilteredSignal, newFs);
    [eraSignalFft, ~] = singleSidedFft(eraSignal, newFs);
    [parallelSignalFft, ~] = singleSidedFft(thisParallelFilteredSignal, newFs);

    t1 = nexttile(til);
    hold(t1, "on")

    plot(t1, timeVector, thisFilteredSignal, 'k', 'LineWidth', 1.5, 'DisplayName', 'Original impulse')
    plot(t1, timeVector, eraSignal, '--r', 'LineWidth', 1.5, 'DisplayName', 'ERA reconstructed impulse')
    plot(t1, timeVector, thisParallelFilteredSignal, '-.g', 'LineWidth', 1.5, 'DisplayName', 'Modal form signal')
    legend

    t2 = nexttile(til2);
    hold(t2, "on")

    plot(t2, freqAxis, 20*log10(filteredSignalFft), 'k', 'LineWidth', 1.5, 'DisplayName', 'Original impulse')
    plot(t2, freqAxis, 20*log10(eraSignalFft), '--r', 'LineWidth', 1.5, 'DisplayName', 'ERA reconstructed impulse')
    plot(t2, freqAxis, 20*log10(parallelSignalFft), '-.g', 'LineWidth', 1.5, 'DisplayName', 'Modal form signal')
    legend

    fmin = 0;
    fmax = 3000;
    
    t3 = tic;
    [freq, growth, amp, phase]=fdm_FAST(thisFilteredSignal', newFs, 0, 3000, 5);

    ampE = 2*amp;      % correct for 2*cos(X) = exp(+j*X) + exp(-j*X)
    alpE = -growth;       % flip growth to obtain attenuation
    phaE = phase + 90;   % correct for sin rather than cos
    ii = find(phaE > 180);  % ensure phase lying in [-180,180]
    phaE(ii) = phaE(ii) - 360;

    % REMOVING SPURIOUS COMPONENTS %%%%%
    freES = freq; ampES = ampE; alpES = alpE; phaES = phaE; 
    ii = find(freES > fmin & freES < fmax); % within analysis band
    freES = freES(ii); ampES = ampES(ii); alpES = alpES(ii); phaES = phaES(ii); 
    % ii = find(ampES > 1E-6);                % significant amplitude
    % freES = freES(ii); ampES = ampES(ii); alpES = alpES(ii); phaES = phaES(ii);
    nModesRetained = length(amplitudes);
    ii = find(alpES > 0);                   % not growing
    freES = freES(ii); ampES = ampES(ii); alpES = alpES(ii); phaES = phaES(ii);

    % To draw comparison between ERA and FDM, first perform ERA
    % analysis to determine the number of frequencies to inspect. Once this
    % has been completed, retain that number of frequencies from the FDM
    % method.
    try
        freES = freES(1 : nModesRetained);
        ampES = ampES(1 : nModesRetained);
        alpES = alpES(1 : nModesRetained);
        phaES = phaES(1 : nModesRetained);
    catch
        print('unlucky')
    end

    nSigs = length(freES);
    phaES = deg2rad(phaES);

    [~, fdmSignal, ~, ~] = createImpulseResponse(nSigs, nSamples, newFs, Frequency=freES, ...
        alpha=alpES, amplitude=ampES, phase=phaES, SNR=0);
    t4 = toc;

    fdmDur = t4-t3;
    sprintf("FDM time duration: " + string(fdmDur))

    [fdmFft, fdmFreqAxis] = singleSidedFft(fdmSignal, newFs);

    plot(t1, timeVector, fdmSignal, '--m', 'LineWidth', 1.3, 'DisplayName', 'FDM Modal Resynthesis')
    plot(t2, fdmFreqAxis, 20*log10(fdmFft), '--m', 'LineWidth', 1.3, 'DisplayName', 'FDM Modal Resynthesis')

    allTimes(:, iFile) = [eraDur, fdmDur];
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