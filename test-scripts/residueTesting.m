clear;
clc;

Fs = 10000;
dt = 1/Fs;
Ns = 1000;
t = (0:(Ns-1))*dt;
frequencyAxis = (0 : (Ns / 2)) * (Fs / Ns);

freqs = 100;

amplitudes = 1;

alpha = 80;

phases = 60;

sig = amplitudes * exp(-alpha * t) .* ...
    sin(2 * pi * freqs * t + deg2rad(phases));

opts = eraConfig();
opts.returnDebug = true;
svdTolerance = 1;

[outputSignal, frequencies, amps, phi, decayFactors, modes, debug] = ...
    eigensystemRealisation(sig, svdTolerance, Fs, opts);

figure
tiledlayout(3, 1)
nexttile
hold on
plot(sig, 'LineWidth', 1.5)
plot(outputSignal, '--', 'LineWidth', 1.5)
title('Simulated single mode impulse and ERA reconstructed single mode')
%% True impulse reconstruction from simulated parameters

Ri = exp(-alpha * dt);
Si = sin(freqs * 2 * pi * dt);
Ci = cos(freqs * 2 * pi * dt);
phi = deg2rad(phases);
Ai = amplitudes;

b0 = Ai * sin(phi);
b1 = Ai * Ri * (cos(phi) * Si - sin(phi)*Ci);

a1 = 2 * Ri * Ci;
a2 = Ri ^ 2;

den = [b0 b1];
num = [1 -a1 a2];

inputSignal = [1 zeros(1, Ns - 1)];
test = filter(den, num, inputSignal);

% Plot every 5th sample so it looks better
nexttile
stem(test(1 : 5 : end), '.-', 'LineWidth', 1.5, 'DisplayName', 'True impulse response');
hold on

%% Use phase and amplitude directly from residues.

% Take the phase and angle directly from the residue terms and use the ERA
% parameters.
ang = angle(debug.Modal.Ct(1) * debug.Modal.Bt(1));
amp = 2 * abs(debug.Modal.Ct(1) * debug.Modal.Bt(1));

phaseCos = ang;
phaseSin = phaseCos + pi/2;
phaseDeg = rad2deg(phaseSin);

Ri = exp(-debug.Modal.allDecayRates(1) * dt);
Si = sin(frequencies * 2 * pi * dt);
Ci = cos(frequencies * 2 * pi * dt);

b0 = amp * sin(phaseSin);
b1 = amp * Ri * (cos(phaseSin) * Si - sin(phaseSin)*Ci);

a1 = 2 * Ri * Ci;
a2 = Ri ^ 2;

b = [b0 b1];
a = [1 -a1 a2];

test = filter(b, a, inputSignal);

stem(test(1 : 5 : end), 'LineWidth', 1.5, 'DisplayName', 'Residue impulse response')


%% Derivation of phase and ampitude from impulse invariant method.

Ri = exp(-debug.Modal.allDecayRates(1) * dt);
Si = sin(frequencies * 2 * pi * dt);
Ci = cos(frequencies * 2 * pi * dt);

b0i = 2 * real(debug.Modal.Ct(1) * debug.Modal.Bt(1));
b1i = -2 * real(debug.Modal.Ct(1) * debug.Modal.Bt(1) * debug.Modal.At(2,2));

phi = acot((b1 / (b0 * Ri * Si)) + Ci / Si);
newAmp = b0 / sin(phi);

a1i = 2 * Ri * Ci;
a2i = Ri ^ 2;

bi = [b0i b1i];
ai = [1 -a1i a2i];

test2 = filter(bi, ai, inputSignal);
stem(test2(1 : 5 : end), '--', 'LineWidth', 1.5, 'DisplayName', 'IIM impulse response')
title('Comparison of true impulse response, residue + pole impulse response and IIM derived phase and amplitude impulse response')
legend

%%
nexttile
stem(test(1 : 100))
hold on
stem(test2(1 : 100))
title('Comparison of ERA parameter impulse response and Residue + Pole method')

%% Divide the residue by a factor of one pole.

resi = debug.Modal.Ct(1) * debug.Modal.Bt(1);

newResi = resi / debug.Modal.At(1, 1);

testPhase = rad2deg(angle(newResi) + pi / 2);

testAmp = 2 * abs(newResi);

%% Test with factor of d built-in

d = sig(1);

b0new = d;

b1new = (2 * real(debug.Modal.Ct(1) * debug.Modal.Bt(1))) - 2*d*real(debug.Modal.At(1, 1));

b2new = d * abs(debug.Modal.At(1, 1))^2 - 2*real(resi * debug.Modal.At(2, 2));

a0new = 1;

a1new = -2*real(debug.Modal.At(1, 1));

a2new = abs(debug.Modal.At(1, 1))^2;

bnew = [b0new b1new b2new];
anew = [a0new a1new a2new];

test3 = filter(bnew, anew, inputSignal);

figure
plot(sig, 'LineWidth', 2)
hold on
plot(test3, 'r--', 'LineWidth', 1.5)

%% Test with higher order of modes.

fs = 1000;
time = 2;
nSamples = (time * fs);
freqs = [2 5 10];
nSignals = length(freqs);
alpha = [5 5 5];

opts = eraConfig();
opts.returnDebug = true;
% opts.residueScaling = true;

[allSignals, outputSignal, timeVector, debug] = ...
    createImpulseResponse(nSignals, nSamples, fs, Frequency=freqs, Alpha=alpha, SNR=0);

svdTolerance = 3;

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

nFilters = length(poles) / 2;
filterIndex = 1 : 0.5 : nFilters;
d = outputSignal(1);

b = zeros(nFilters, 3);
a = zeros(nFilters, 3);
modes = zeros(nFilters, nSamples);
inputSignal = [1 zeros(1, nSamples - 1)];

% Try two alternative methods - 1. factor d into the first modes'
% difference equation.
for i = 1 : 2 : length(poles)

    thisFilter = filterIndex(i);

    if i == 1
    
        b0i = d;
        b1i = (2 * real(residues(i))) - 2*d*real(poles(i));
        b2i = d * abs(poles(i))^2 - 2*real(residues(i) * poles(i + 1));

    else

        b0i = 0;
        b1i = 2 * real(residues(i));
        b2i = -2 * real(residues(i) * poles(i + 1));

    end

    a0new = 1;
    a1new = -2 * real(poles(i));    
    a2new = abs(poles(i)) ^ 2;

    thisB = [b0i, b1i, b2i];
    thisA = [a0new, a1new, a2new];

    b(thisFilter, :) = thisB;
    a(thisFilter, :) = thisA;
    modes(thisFilter, :) = filter(thisB, thisA, inputSignal);

end

% sosFilters = tf2sos(b, a, nFilters);
filterSignalOut = modes(1, :) + modes(2, :) + modes(3, :);

nexttile
plot(timeVector, eraSignal, 'k', 'LineWidth', 1.5, 'DisplayName', 'ERA signal')
hold on
plot(timeVector, filterSignalOut, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Filter sections signal')
legend

[zb, za] = residue(residues, poles, d);

zrs = roots(zb);
poles_check = roots(za);

nexttile
zplane(zeros, poles_check)

[B, A] = buildFilterCoefficients(poles, residues, outputSignal(1));
parallelFilterSignal = parallelFilter(B, A, inputSignal);

nexttile
plot(timeVector, filterSignalOut, 'r', 'LineWidth', 1.5, 'DisplayName', 'Filter sections signal')
hold on
plot(timeVector, parallelFilterSignal, '.-b', 'LineWidth', 1.5, 'DisplayName', 'Parallel filter function signal')


%% Equal factor of d for all filters.

bc = zeros(nFilters, 3);
ac = zeros(nFilters, 3);
modes = zeros(nFilters, nSamples);
inputSignal = [1 zeros(1, nSamples - 1)];

% d = d / nFilters;

% Try two alternative methods - 1. factor d into all modes.
for i = 1 : 2 : length(poles)

    thisFilter = filterIndex(i);
    
    b0i = d;
    b1i = (2 * real(residues(i))) - 2*d*real(poles(i));
    b2i = d * abs(poles(i))^2 - 2*real(residues(i) * poles(i + 1));

    a0new = 1;
    a1new = -2 * real(poles(i));    
    a2new = abs(poles(i)) ^ 2;

    thisB = [b0i, b1i, b2i];
    thisA = [a0new, a1new, a2new];

    bc(thisFilter, :) = thisB;
    ac(thisFilter, :) = thisA;
    modes(thisFilter, :) = filter(thisB, thisA, inputSignal);

end

% sosFilters = tf2sos(b, a, nFilters);
filterSignalOut = modes(1, :) + modes(2, :) + modes(3, :);

nexttile
plot(timeVector, eraSignal, 'k', 'LineWidth', 1.5, 'DisplayName', 'ERA signal')
hold on
plot(timeVector, filterSignalOut, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Filter sections signal')
legend

[b, a] = residue(residues, poles, d);

zeros = roots(bc);
poles_check = roots(ac);

nexttile
zplane(zeros, poles_check)

%% 

[btest, atest] = buildFilterCoefficients(poles, residues, d);

outputSignal = parallelFilter(btest, atest, inputSignal);

nexttile
plot(timeVector, eraSignal, 'k', 'LineWidth', 1.5, 'DisplayName', 'ERA signal')
hold on
plot(timeVector, outputSignal, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Filter sections signal')
legend