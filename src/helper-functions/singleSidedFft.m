function [frequencySpectrum, frequencyAxis] = singleSidedFft(signal, fs)


nSamples = numel(signal);
fftSignal = fft(signal);
P2 = abs(fftSignal / nSamples);
P1 = P2(1 : (nSamples / 2 + 1));
P1(2 : end - 1) = 2 * P1(2 : end - 1);
frequencySpectrum = P1;
frequencyAxis = (0 : nSamples / 2) * (fs / nSamples);

end