function [allSignals, outputSignal, timeVector, debug] = ...
    createImpulseResponse(nSignals, nSamples, fs, options)
% createImpulseResponse Vectorised method for creating impulse responses. 
% Will generate uniformly distributed random frequencies, amplitudes, and 
% decay rates if input arguments are empty. Phase will be set at either 0
% or pi / 2.

% Frequency [Hz]
% Amplitudes [1]
% Phases [Rad]
% Alpha [Nepers] {s^-1}
% Time vector [s]

%% Input checking.
arguments
    nSignals   (1, 1) double {mustBeInteger, mustBePositive, mustBeFinite}
    nSamples   (1, 1) double {mustBeInteger, mustBePositive, mustBeFinite}
    fs         (1, 1) double {mustBePositive, mustBeFinite}
    options.Frequency  (:, 1) double {mustBeReal, mustBeFinite, ...
                             mustBeNonnegative, ...
                             mustBeNumelOrEmpty(options.Frequency, ...
                             nSignals)} = []
    options.Amplitude  (:, 1) double {mustBeReal, mustBeFinite, ...
                             mustBeNumelOrEmpty(options.Amplitude, ...
                             nSignals)} = []
    options.Phase      (:, 1) double {mustBeReal, mustBeFinite, ...
                             mustBeNumelOrEmpty(options.Phase, ...
                             nSignals)} = []
    options.Alpha      (:, 1) double {mustBeReal, mustBeFinite, ...
                       mustBeNonnegative, ...
                       mustBeNumelOrEmpty(options.Alpha, nSignals)} = []
    options.SNR        (1 ,1) double {mustBeReal, mustBeFinite} = 0
    options.Debug      (1, 1) logical = false
end

%% Defaults.

% Note: come back and make these array input options. This is fine for now.
MIN_FREQUENCY = 100;
MAX_FREQUENCY = 2500;

MIN_AMPLITUDE = 0.5;
MAX_AMPLITUDE = 1;

MAX_PHASE = pi / 2;
PHASE_LIMIT = 0.5;

MIN_ALPHA = 20;
MAX_ALPHA = 100;

DB_SCALING = 10;
POWER_FACTOR = 10;

%% Random signal parameter generation.

% Reassign for neatness.
frequency = options.Frequency;
amplitude = options.Amplitude;
phase = options.Phase;
alpha = options.Alpha;
SNR = options.SNR;
debugFlag = options.Debug;

% Randomly generate parameters when not provided.
if isempty(frequency)
    frequency = scaledRand(MIN_FREQUENCY, MAX_FREQUENCY, nSignals);
end

if isempty(amplitude)
    amplitude = scaledRand(MIN_AMPLITUDE, MAX_AMPLITUDE, nSignals);
end

if isempty(phase)
    phase = (rand(nSignals, 1) >= PHASE_LIMIT) * MAX_PHASE;
end

if isempty(alpha)
    alpha = scaledRand(MIN_ALPHA, MAX_ALPHA, nSignals);
end

%% Build signals.

% Time vector must be a row to align with column vector parameters above.
timeVector = (0 : nSamples - 1) / fs;

allSignals = amplitude .* exp(-alpha .* timeVector) .* ...
    sin((2 * pi .* frequency .* timeVector) + phase);

outputSignal = sum(allSignals, 1);

% Add gaussian noise.
if SNR ~= 0
    signalEnergy = sum(outputSignal .^ 2);
    targetEnergy = signalEnergy / (DB_SCALING ^ (SNR / POWER_FACTOR));
    
    noiseSignal = randn(1, nSamples);
    noiseEnergy = sum(noiseSignal .^ 2);
    noiseAmplitudeScaling = sqrt(targetEnergy / noiseEnergy);
    finalNoiseSignal = noiseAmplitudeScaling * noiseSignal;
    outputSignal = outputSignal + finalNoiseSignal;
end

%% Debug info.
debug = struct();
debug.frequency     = [];
debug.amplitude     = [];
debug.phase         = [];
debug.alpha         = [];
debug.signalEnergy  = [];
debug.targetEnergy  = [];
debug.noiseEnergy   = [];
debug.noiseSignal   = [];

if debugFlag
    % In the event any of the parameters are randomly created, assign to
    % the debug struct to inspect.
    debug.frequency     = frequency;
    debug.amplitude     = amplitude;
    debug.phase         = phase;
    debug.alpha         = alpha;
    debug.signalEnergy  = signalEnergy;
    debug.targetEnergy  = targetEnergy;
    debug.noiseEnergy   = noiseEnergy;
    debug.noiseSignal   = finalNoiseSignal;

end

end

%% Subfunctions.
function output = scaledRand(lowerLimit, upperLimit, nElements)

output = lowerLimit + (upperLimit - lowerLimit) .* rand(nElements, 1);

end

function mustBeNumelOrEmpty(value, n)

if ~isempty(value) && numel(value) ~= n
    error('createImpulseResponse:wrongNumel', ...
        'Must be empty or have %d elements.', n);
end

end