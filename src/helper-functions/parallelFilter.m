function outputSignal = parallelFilter(B, A, signal)

[nRows, ~] = size(B);
nSamples = numel(signal);

outputSignal = zeros(nRows, nSamples);

parfor iRow = 1 : nRows

    outputSignal(iRow, :) = filter(B(iRow, :), A(iRow, :), signal)

end

outputSignal = sum(outputSignal, 1);

end