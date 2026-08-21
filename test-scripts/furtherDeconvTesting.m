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
impulseStartSample = [24820
                      24470
                      24770
                      24600
                      25180
                      25380
                      25030
                      24540
                      25010
                      24890];

impulseEndIndex = impulseStartSample + 30000;

% Define impulse hammer start sampled and noise periods for wiener filter.
hammerStartSample = [24677
                     24321
                     24622
                     24470
                     25042
                     25241
                     24889
                     24408
                     24871
                     24750];

hammerEndSample = hammerStartSample + 40;

hammerStartSample = hammerStartSample - 20;

noiseEndIndex = 20000;

% Resampling, wiener filter and deconvolution parameters.
originalFs = 100000;
resampleFactor = 10;
newFs = originalFs / resampleFactor;
nTaps = 500;
lambda = 1e-2;

% ERA parameters;
opts = eraConfig();
opts.returnDebug = true;
opts.zeroFrequencyFilter = false;
svdTolerance = 80;
impulseErrors = 1 : length(svdTolerance);
frfError = 1 : length(svdTolerance);

% Figures.
fig = figure;
til = tiledlayout(5, 2);
ylabel(til, 'Amplitude')
xlabel(til, 'Time (s)')

fig2 = figure;
til2 = tiledlayout(5, 2);
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
    thisHammerStartSample = hammerStartSample(iFile);
    thisHammerEndSample = hammerEndSample(iFile);

    thisSignal = thisFileData(:, 5);
    thisSignal = thisSignal - mean(thisSignal);
    % Split up noise and filter.
    thisNoiseSignal = thisSignal(1 : noiseEndIndex);
    % Truncate the signal to only the impulse response data.
    thisSignal = thisSignal(thisStartSample : thisEndSample);
    thisFilteredSignal = wienerFilter(thisSignal, thisNoiseSignal, nTaps);

    % Apply the same process to the hammer impact signal.
    thisImpactSignal = thisFileData(:, 1);
    thisImpactSignal = thisImpactSignal - mean(thisImpactSignal);
    figure
    plot(thisImpactSignal)
    thisImpactNoiseSignal = thisImpactSignal(1 : noiseEndIndex);
    thisImpactSignal = thisImpactSignal(thisHammerStartSample : thisHammerEndSample);
    thisFilteredImpactSignal = wienerFilter(thisImpactSignal, ...
        thisImpactNoiseSignal, nTaps);
    
    % Frequency domain deconvolution.
    thisFilteredSignal = deconvolveFreqDomain(thisFilteredImpactSignal, ...
        thisFilteredSignal, lambda);
    [thisFRF, freqAx] = singleSidedFft(thisFilteredSignal, 100000);

    % Resample the filtered signal for processing.
    thisFilteredSignal = resample(thisFilteredSignal, 1, resampleFactor);
    nSamples = length(thisFilteredSignal);
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
    ii = find(ampES > 1E-6);                % significant amplitude
    freES = freES(ii); ampES = ampES(ii); alpES = alpES(ii); phaES = phaES(ii);
    nModesRetained = length(amplitudes);
    ii = find(alpES > 0);                   % not growing
    freES = freES(ii); ampES = ampES(ii); alpES = alpES(ii); phaES = phaES(ii);
    phaES = deg2rad(phaES);
    % To draw comparison between ERA and FDM, first perform ERA
    % analysis to determine the number of frequencies to inspect. Once this
    % has been completed, retain that number of frequencies from the FDM
    % method.
    % try
    %     freES = freES(1 : nModesRetained);
    %     ampES = ampES(1 : nModesRetained);
    %     alpES = alpES(1 : nModesRetained);
    %     phaES = phaES(1 : nModesRetained);
    % catch
    %     print('unlucky')
    % end

    [b, a] = generateIIMcoefficients(freES, alpES, ampES, phaES, 1/newFs);

    parallelFilterFDM = parallelFilter(b, a, inputSignal);

    t4 = toc;

    fdmDur = t4-t3;
    sprintf("FDM time duration: " + string(fdmDur))

    [fdmFft, fdmFreqAxis] = singleSidedFft(parallelFilterFDM, newFs);

    plot(t1, timeVector, parallelFilterFDM, '--m', 'LineWidth', 1.3, 'DisplayName', 'FDM parallel filter')
    plot(t2, fdmFreqAxis, 20*log10(fdmFft), '--m', 'LineWidth', 1.3, 'DisplayName', 'FDM parallel filter')

    allTimes(:, iFile) = [eraDur, fdmDur];
    clear t1 t2 t3 t4;

    eraNMSEError(iFile) = sum((thisFilteredSignal(:) - ...
        thisParallelFilteredSignal(:)) .^ 2) / sum(thisFilteredSignal(:) .^ 2);
    fdmNMSEError(iFile) = sum((thisFilteredSignal(:) - ...
        parallelFilterFDM(:)) .^ 2) / sum(thisFilteredSignal(:) .^ 2);
    sprintf("ERA Parallel Filter NMSE Error: " + string(round(eraNMSEError(iFile), 4)) ...
        + "\nFDM Parallel Filter NMSE Error: " + string(round(fdmNMSEError(iFile), 4)))

end

t = table(eraNMSEError', fdmNMSEError', 'VariableNames', {'ERA NMSE Error', 'FDM NMSE Error'})


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

function freqError = nmseFreq(originalSignal, eraSignal)

fftSig = fft(originalSignal);
fftEra = fft(eraSignal);
fftSig = fftSig(:);
fftEra = fftEra(:);

freqError = sum((abs(fftSig) - abs(fftEra)).^2) / sum(abs(fftSig).^2);

end

function [b, a] = generateIIMcoefficients(frequency, alpha, amplitude, ...
    phase, dt)

a0 = ones(length(frequency), 1);

R = exp(-alpha .* dt);
S = sin(2*pi*frequency.*dt);
C = cos(2*pi*frequency.*dt);

b0 = amplitude .* sin(phase);
b1 = amplitude .* R .* (cos(phase) .* S - sin(phase) .* C);
a1 = -2 .* R .* C;
a2 = R .^2;

b = [b0 b1];
a = [a0 a1 a2];

end

function admittanceIR = deconvolveFreqDomain(hammerInput, bridgeVelocity, lambda)
    hammerInput    = hammerInput(:)    - mean(hammerInput);
    bridgeVelocity = bridgeVelocity(:) - mean(bridgeVelocity);

    % Zero-pad hammer to match bridge input.
    N = numel(bridgeVelocity);
    nfft = 2 ^ nextpow2(2 * N);

    H = fft(hammerInput, nfft);
    Y = fft(bridgeVelocity, nfft);

    % Tikhonov deconvolution.
    frf = (Y .* conj(H)) ./ (abs(H) .^ 2 + lambda);
    admittanceIR = real(ifft(frf));
end

function [admittanceIR, convergence, yHat] = lmsDeconvolution( ...
    hammerInput, bridgeVelocity, mu, nTaps, nEpochs)

    nSamples = numel(bridgeVelocity);

    H           = zeros(nTaps, 1);   % filter weights = impulse response
    convergence = zeros(nEpochs, 1);
    yHat        = zeros(nSamples, 1);

    for epoch = 1 : nEpochs

        xBuffer  = zeros(nTaps, 1);   
        sqError  = 0;

        for n = 1 : nSamples

            % Shift new input sample into the delay line
            xBuffer = [hammerInput(n); xBuffer(1:end-1)];

            % Model output: current weights applied to input history
            yHat(n) = H' * xBuffer;

            % Error between measured velocity and model output
            e = bridgeVelocity(n) - yHat(n);

            % LMS weight update
            H = H + 2 * mu * e * xBuffer;

            sqError = sqError + e^2;

        end

        convergence(epoch) = sqError / nSamples;

    end

    % The converged weights ARE the admittance impulse response
    admittanceIR = H;

end