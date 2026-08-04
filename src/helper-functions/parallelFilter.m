function outputSignal = parallelFilter(B, A, signal)

[nRows, ~] = size(B);
nSamples = numel(signal);

outputSignal = zeros(1, nSamples);

for iRow = 1 : nRows

    outputSignal = outputSignal + filter(B(iRow, :), ...
        A(iRow, :), signal);

end

end