function [outputSignal, allSignals] = parallelFilter(B, A, signal)

[nRows, ~] = size(B);
nSamples = numel(signal);

outputSignal = zeros(1, nSamples);
allSignals = zeros(nRows, nSamples);

for iRow = 1 : nRows
    
    thisSignal = filter(B(iRow, :), A(iRow, :), signal);
    outputSignal = outputSignal + thisSignal;
    allSignals(iRow, :) = thisSignal;

end

end