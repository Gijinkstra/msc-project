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

%% Mode Similarity Index

H1 = debug.Hankel.shiftedHankelMatrix;
nColumns = debug.Hankel.nColumns;
nRows = debug.Hankel.nRows;
m = svdTolerance;

[U1, S1, V1] = svd(H1);

Um1 = U1(:, 1 : m);
Sm1 = S1(1 : m, 1 : m);
Vm1 = V1(:, 1 : m);

columnBasisVector = zeros(nColumns, 1);
columnBasisVector(1) = 1;

rowBasisVector = zeros(nRows, 1);
rowBasisVector(1) = 1;

% Half power and inverse matrices. Transform to avoid computing over the
% entire S matrix. Only saves time for large S (assumed).
singularValues  = diag(Sm1);
sHalf = diag(sqrt(singularValues));       
sHalfInv = diag(1 ./ sqrt(singularValues)); 

% Eq. (22-24).
A1 = sHalfInv * Um1' * H1 * Vm1 * sHalfInv;

[psi, lambda] = eig(A1);

% Discrete to continuous time conversion
continuousEigenValues = log(lambda) ./ dt;

naturalFrequencies = abs(continuousEigenValues);

% Extract frequencies from the eigenvalues.
shiftedFrequencies = naturalFrequencies ./ (2*pi);

% Extract the decay rate.
shiftedDecayRates = -real(continuousEigenValues);

% Extract damping factor from the eigenvalues.
shiftedDampingFactors = -real(continuousEigenValues) ./ ...
    naturalFrequencies;

figure
plot(frequencies, dampingFactors, 'o', 'LineWidth', 1.4, ...
    'DisplayName', 'Standard frequencies')
hold on
plot(shiftedFrequencies, shiftedDampingFactors, 'gx', 'LineWidth', 1.4, ...
    'DisplayName', 'Shifted Frequencies')
ylabel('Damping factor')
xlabel('Frequency')
legend