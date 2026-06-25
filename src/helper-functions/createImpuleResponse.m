function outputSignal = createImpuleResponse(frequencies, amplitudes, ...
    phases, alpha, timeVector)
% Vectorised method for creating impulse responses. Will generate random
% amplitudes, phases, and decay rates if input arguments are empty.
nFrequencies = numel(frequencies);
nSamples = numel(timeVector);
allPhases = repmat(phases', 1, nSamples);

allSignals = amplitudes' .* exp(-alpha' * timeVector) .* ...
    sin((2 * pi .* frequencies' * t) + ((pi / 180) .* allPhases));

outputSignal = sum(allSignals, 1);

end