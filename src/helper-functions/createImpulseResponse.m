function [allSignals, outputSignal, debug] = ...
    createImpulseResponse(nSignals, frequencies, amplitudes, phases, ...
    alpha, timeVector, SNR, debugFlag)
% Vectorised method for creating impulse responses. Will generate random
% amplitudes, phases, and decay rates if input arguments are empty.

% Frequency [Hz]
% Amplitudes [1]
% Phases [Rad]
% Alpha [Nepers] {s^-1}
% Time vector [s]

if isempty(nSignals) && isempty(frequencies)

    error('Input either a number of signals, or specify frequencies')

end

% When amplitude is empty, create a random amplitude between 0.5 and 1.
if isempty(amplitudes)

    amplitudes = 0.5 + (1 - 0.5) .* rand(nSignals, 1);

end

% When phase is empty, create a random array of either 0 or pi/2 phase.
if isempty(phases)

    phases = rand(nSignals, 1);
    phases(phases >= 0.5) = pi / 2;
    phases(phases < 0.5) = 0;
   
end

% When empty create a decay value between 20 and 100.
if isempty(alpha)

    alpha = 20 + (100 - 20) .* rand(nSignals, 1);

end

if isempty(debugFlag)

    debugFlag = false;

end

frequencies = frequencies(:);
alpha = alpha(:);
amplitudes = amplitudes(:);
phases = phases(:);

nSamples = numel(timeVector);
allPhases = repmat(phases, 1, nSamples);

allSignals = amplitudes .* exp(-alpha .* timeVector) .* ...
    sin((2 * pi .* frequencies .* timeVector) + allPhases);

outputSignal = sum(allSignals, 1);

if ~isempty(SNR)

    Es = sum(outputSignal .^ 2);
    En = Es / (10 ^ (SNR / 10));
    sign0 = 1 - 2 * rand(1, nSamples);
    En0 = sum(sign0 .^ 2);
    an = sqrt(En / En0);
    sign = an * sign0;
    outputSignal = outputSignal + sign;

end

%% Debug info.
if debugFlag

    debug = struct();
    
    modeIndex = 1 : nSignals;

    debug.table = table(modeIndex', frequencies, alpha, ...
        amplitudes, phases, 'VariableNames', {'Mode index', ...
        'Frequency (Hz)', 'Decay rate (s^-1)', 'Amplitude (1)', ...
        'Phase (rads)'});

end