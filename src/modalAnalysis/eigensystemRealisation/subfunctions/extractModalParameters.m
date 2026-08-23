function [filteredFrequencies, filteredDampingFactors, ...
    filteredAmplitudes, filteredPhase, filteredModes, ...
    filteredDecayRates, filteredAt, filteredBt, filteredCt, varargout] = ...
    extractModalParameters(A, B, C, timeStep, opts)
% EXTRACTMODALPARAMETERS  Extracts modal filtered parameters and  filtered
%                         modal matrices from ERA state space matrices.
%
% DESCRIPTION:
%   Recovers natural frequencies, decay rates, damping ratios, and mode
%   shapes from the discrete time state transition matrix A, input matrix
%   B, and output matrix C identified by ERA.
%
% INPUTS:
%   A        - [double, n x n]  Discrete time state transition matrix.
%   B        - [double, n x 1]  Input vector.
%   C        - [double, 1 x n]  Output vector.
%   timeStep - [double, scalar] Sampling interval in seconds (1/fs).
%   opts     - [struct]         Configuration struct from eraConfig().
%
% OPTIONS (opts fields consumed by this function):
%   filterUnstable     - [logical] Remove poles with |lambda_d| >= 1. [true]
%   filterNegativeDamp - [logical] Remove modes with damping < 0. [true]
%                        Negative damping corresponds to growing
%   filterNegativeFreq - [logical] Remove modes with imag(lambda) < 0. [true]
%                        Retains only positive imaginary half of each
%                        conjugate pair to avoid duplicate physical modes.
%                        Note that this is not applied to the state space
%                        matrices At, Bt, Ct.
%                        oscillations and is non-physical.
%   dampingThreshold   - [double]  Upper damping ratio bound. [0.99]
%                        Modes with damping > threshold are considered
%                        overdamped and non-resonant. Set to 1.0 to
%                        retain critically and overdamped modes.
%   returnDiagnostics  - [logical] Return diagnostics struct as 8th
%                        output. [false] Controlled by caller for debug
%                        collection.
%
% OUTPUTS:
%   filteredFrequencies    - [double, k x 1]  Natural frequencies in Hz.
%                            k <= n/2 after filtering conjugate pairs and
%                            computational modes.
%   filteredDecayRates     - [double, k x 1]  Decay rates in Nepers/s.
%   filteredDampingFactors - [double, k x 1]  Damping ratios (non-dimensional).
%                            zeta = decayRate / naturalFrequency.
%   filteredModes          - [double, k x 1]  Observed mode shapes.
%   filteredAt             - [double, 2k x 2k]  Diagonal eigenvalue matrix 
%                            (modal A).
%   filteredBt             - [double, 2k x 1]  Modal input matrix.
%   filteredCt             - [double, 1 x 2k]  Modal output matrix.
%   varargout{1}           - [struct]          Diagnostics struct (if 
%                            requested via opts.returnDiagnostics and
%                            nargout > 7). See NOTES.
%
% NOTES:
%   At, Bt, Ct are returned filtered, but with complex conjugate pairs
%   intact. This is necessary to calculate the impulse response later, but
%   will result in double the appropriate dimension compared to extracted
%   physical information.
%
%   Diagnostics struct fields (when returned):
%     .At                     - Unfiltered diagonalised A matrix 
%     .Bt                     - Unfiltered transformed B matrix (psi^{-1}B)
%     .Ct                     - Unfiltered transformed C matrix (C * psi)
%     .totalEigenvalues       - Total number of eigenvalues identified (n)
%     .retainedCount          - Number of modes after filtering
%     .filteredCount          - Number removed by filtering
%     .filteredOutEigenvalues - Discrete eigenvalues removed
%     .allFrequencies         - All frequencies before filtering [n x 1]
%     .allDecayRates          - All decay rates before filtering [n x 1]
%     .allDampingFactors      - All damping ratios before filtering [n x 1]
%     .filterContributions    - Struct with per-filter removal counts:
%                               .unstable, .negativeFreq, .negativeDamp,
%                               .overdamped, .zeroFreq
%
% REFERENCES:
%   [1] Juang, J.-N. & Pappa, R.S. (1985). An eigensystem realization
%       algorithm for modal parameter identification and model reduction.
%       Journal of Guidance, Control, and Dynamics, 8(5), 620-627.
%   [2] Juang, J.-N. (1994). Applied System Identification.
%       Prentice Hall, Englewood Cliffs, NJ. Chapters 5-7.
%
% AUTHOR:   -
% DATE:     June 2026
% VERSION:  1.2
%
% See also: eigensystemRealisation, stateSpaceMatrixConstruction,
%           buildImpulseResponse, buildModalImpulseResponse
%% Input validation
n = size(A, 1);
assert(size(A, 2) == n, ...
    'extractModalParameters:nonSquareA', 'A must be square.');
assert(size(B, 1) == n, ...
    'extractModalParameters:BDimensionMismatch', ...
    'B rows must match A dimension.');
assert(size(C, 2) == n, ...
    'extractModalParameters:CDimensionMismatch', ...
    'C columns must match A dimension.');

%% Main code
% Eigenvalues and eigenvectors of A.
[psi, At] = eig(A);

% Co-ordinate transformations for modal forms of B and C given A is now
% diagonalised.
Bt = psi \ B;
Ct = C * psi;

% Extract frequencies from the eigenvalues.
[sortedAt, sortedBt, sortedCt, frequencies, dampingFactors, amplitudes, ...
    phases, decayRates, modes, poles, residues] = ...
    extractModeComponents(At, Bt, Ct, timeStep, opts);

%% Stage 1 — Non-conjugate filters (applied to both matrices and parameters)
stage1Mask = true(size(poles));

% Track individual filter contributions for diagnostics
unstableMask   = false(size(poles));
negDampMask    = false(size(poles));
overdampedMask = false(size(poles));

if opts.filterUnstable
    unstableMask = (abs(poles) >= 1);
    stage1Mask   = stage1Mask & ~unstableMask;
end

if opts.filterNegativeDamp
    negDampMask = (dampingFactors < 0);
    stage1Mask  = stage1Mask & ~negDampMask;
end

if opts.filterOverdamped
    overdampedMask = (dampingFactors > opts.dampingThreshold);
    stage1Mask     = stage1Mask & ~overdampedMask;
end

% Remove the 0 frequency components.

if opts.zeroFrequencyFilter
    zeroFreqMask = (imag(poles) == 0);
    stage1Mask = stage1Mask & ~zeroFreqMask;
else
    zeroFreqMask = true(size(poles));
end

% Reshape the mask into a 2*n column vector to check both conjugate pairs
% are removed, then rebuild the original column vector. Note that this
% assumes the eigenvalue decomposition returns the complex conjugate pairs
% in adjacent elements. This should generally be the case for real
% matrices.
conjugateMask = reshape(stage1Mask, 2, []);
conjugateRetained = all(conjugateMask, 1);
stage1Mask = reshape(repmat(conjugateRetained, 2, 1), [], 1);

%% Apply Stage 1 filter to matrices — preserve conjugate pairs
filteredAt = sortedAt(stage1Mask, stage1Mask);
filteredBt = sortedBt(stage1Mask);
filteredCt = sortedCt(stage1Mask);

%% Stage 1.5 - Mode energy order - REMOVED
sortedPoles = poles(stage1Mask);
sortedResidues = residues(stage1Mask);

%% Stage 2 — Conjugate pair reduction (applied only to parameter outputs)
stage2Mask = stage1Mask;
negFreqMask = false(size(poles));

if opts.filterNegativeFreq
    negFreqMask = (imag(poles) < 0);
    stage2Mask  = stage2Mask & ~negFreqMask;
end

%% Apply Stage 2 filter to parameter outputs only
filteredFrequencies = frequencies(stage2Mask);
filteredDecayRates = decayRates(stage2Mask);
filteredDampingFactors = dampingFactors(stage2Mask);
filteredModes = modes(stage2Mask);
filteredPhase = phases(stage2Mask);
filteredAmplitudes = amplitudes(stage2Mask);

%% Optional debug info.
if nargout > 7 && isfield(opts, 'returnDebug') && opts.returnDebug

    diagnostics.At                     = At;
    diagnostics.Bt                     = Bt;
    diagnostics.Ct                     = Ct;
    diagnostics.sortedAt               = sortedAt;
    diagnostics.sortedBt               = sortedBt;
    diagnostics.sortedCt               = sortedCt;
    % diagnostics.totalEigenvalues       = numel(eigenValues);
    diagnostics.retainedCount          = sum(stage2Mask);
    diagnostics.filteredCount          = sum(~stage2Mask);
    % diagnostics.filteredOutEigenvalues = eigenValues(~stage2Mask);
    diagnostics.residues               = sortedResidues;
    diagnostics.poles                  = sortedPoles;
    diagnostics.allFrequencies         = frequencies;
    diagnostics.allDecayRates          = decayRates;
    diagnostics.allDampingFactors      = dampingFactors;

    diagnostics.filterContributions.unstable     = sum(unstableMask);
    diagnostics.filterContributions.negativeFreq = sum(negFreqMask);
    diagnostics.filterContributions.negativeDamp = sum(negDampMask);
    diagnostics.filterContributions.overdamped   = sum(overdampedMask);
    diagnostics.filterContributions.zeroFreq     = sum(zeroFreqMask);
    varargout{1} = diagnostics;

end

end

%% Subfunctions
function [sortedAt, sortedBt, sortedCt, frequency, dampingFactor, ...
    amplitude, phase, decayRate, modes, poles, residues] = ...
    extractModeComponents(At, Bt, Ct, timeStep, opts)
% Input check.
if isfield(opts, 'residueScaling')
    residueScalingFlag = opts.residueScaling;
else
    residueScalingFlag = false;
end

% Constants.
TWO_PI = 2 * pi;

% Necessary pole scaling factor to recover the proper phase and amplitudes.
POLE_SCALING_FACTOR = opts.hankelSignalIndex - 1;

% Used to scale amplitudes for complex conjugate pairs.
AMPLTIUDE_SCALING_FACTOR = 2;

%% Main code.
eigenValues = diag(At);

% Sort the eigenvalues based on their imaginary components.
[~, sortIndex] = sort(abs(imag(eigenValues)));
poles = eigenValues(sortIndex);

sortedAt = At(sortIndex, sortIndex);
sortedBt = Bt(sortIndex);
sortedCt = Ct(sortIndex);

% Discrete to continuous time conversion
continuousEigenValues = log(poles) ./ timeStep;

naturalFrequencies = abs(continuousEigenValues);

frequency = naturalFrequencies ./ TWO_PI;

% Extract the decay rate.
decayRate = -real(continuousEigenValues);

% Extract damping factor from the eigenvalues.
dampingFactor = -real(continuousEigenValues) ./ ...
    naturalFrequencies;

% Modes are the Ct matrix. Reassign here just for clarity and transpose for
% consistency of output dimensions.
modes = sortedCt.';

% Determine the amplitude and phase of the mode via the residues of the
% signal.
residues = sortedCt.' .* sortedBt;

if residueScalingFlag
    % Divide through by a multiple of the pole, depending on Hankel matrix
    % index.
    correctedResidues = residues ./ ...
        (poles .^ POLE_SCALING_FACTOR);
else
    correctedResidues = residues;
end

% Extract the phase and the magnitude of the residues.
phase = angle(correctedResidues);

% Amplitude is doubled due to complex conjugate.
amplitude =  AMPLTIUDE_SCALING_FACTOR * abs(correctedResidues);

end