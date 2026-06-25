function [filteredFrequencies, filteredDampingFactors, ...
    filteredModes, filteredDecayRates, filteredAt, filteredBt, ...
    filteredCt, varargout] = extractModalParameters(A, B, C, timeStep, opts)
% EXTRACTMODALPARAMETERS  Extracts modal filtered parameters and  filtered
%                         modal matrices from ERA state space matrices.
%
% DESCRIPTION:
%   Recovers natural frequencies, decay rates, damping ratios, and mode
%   shapes from the discrete time state transition matrix A, input matrix
%   B, and output matrix C identified by ERA.
%
%   Discrete time eigenvalues are mapped to continuous time via:
%       s_r = log(lambda_r) / timeStep
%
%   Modal parameters are extracted as:
%       frequency  = |s_r| / (2*pi)        [Hz]
%       decayRate  = -real(s_r)            [Nepers/s]
%       damping    = -real(s_r) / |s_r|    [non-dimensional]
%       mode shape = C * eigenvectors      [observed at output]
%
%   The coordinate transformation to modal form diagonalises A:
%       At = psi^{-1} * A * psi = Lambda  [diagonal eigenvalue matrix]
%       Bt = psi^{-1} * B                 [modal input matrix]
%       Ct = C * psi                      [modal output matrix / mode shapes]
%
%   Two-stage filtering is applied:
%
%   Stage 1 — Non-physical mode removal (applied to all outputs):
%     - Unstable poles (|lambda_d| >= 1)
%     - Negative damping (growing oscillations)
%     - Overdamped modes (damping > threshold)
%     - Zero-frequency modes
%
%   Stage 2 — Conjugate pair reduction (applied only to parameter outputs):
%     - Negative imaginary part eigenvalues removed
%     - Parameters become 1 per physical mode (not 2)
%
%   The modal matrices At, Bt, Ct retain both halves of each conjugate
%   pair after stage 1 filtering.
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
%   The principal branch of the complex logarithm is used in the
%   discrete-to-continuous conversion. This maps frequencies above the
%   Nyquist boundary into aliased estimates. For ERA on impulse responses
%   sampled above Nyquist this is not an issue, but it is worth noting
%   for non-standard applications.
%
% REFERENCES:
%   [1] Juang, J.-N. & Pappa, R.S. (1985). An eigensystem realization
%       algorithm for modal parameter identification and model reduction.
%       Journal of Guidance, Control, and Dynamics, 8(5), 620-627.
%   [2] Juang, J.-N. (1994). Applied System Identification.
%       Prentice Hall, Englewood Cliffs, NJ. Chapters 5-7.
%   [3] Pappa, R.S., Elliott, K.B. & Schenk, A. (1993). Consistent-mode
%       indicator for the eigensystem realization algorithm. Journal of
%       Guidance, Control, and Dynamics, 16(5), 852-858.
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

%% Constants
TWO_PI = 2 * pi;

%% Main code
% Eigenvalues and eigenvectors of a linear transformation of A are equal to
% the eigenvalues and eigenvectors of A.
[psi, At] = eig(A);
eigenValues = diag(At);

% Co-ordinate transformations for modal forms of B and C given A is now
% diagonalised.
Bt = psi \ B;
Ct = C * psi;

% Discrete to continuous time conversion
% lambda_d = exp(-lambda_c*timeStep);
% lambda_c = log(lambda_d) / timeStep
continuousEigenValues = log(eigenValues) ./ timeStep;

naturalFrequencies = abs(continuousEigenValues);

% Extract frequencies from the eigenvalues.
frequencies = naturalFrequencies ./ TWO_PI;

% Extract the decay rate.
decayRates = -real(continuousEigenValues);

% Extract damping factor from the eigenvalues.
dampingFactors = -real(continuousEigenValues) ./ ...
    naturalFrequencies;

% Modes are the Ct matrix. Reassign here just for clarity and transpose for
% consistency of output dimensions.
modes = Ct.';

%% Stage 1 — Non-conjugate filters (applied to both matrices and parameters)
stage1Mask = true(size(eigenValues));

% Track individual filter contributions for diagnostics
unstableMask   = false(size(eigenValues));
negDampMask    = false(size(eigenValues));
overdampedMask = false(size(eigenValues));
zeroFreqMask   = false(size(eigenValues));

if opts.filterUnstable
    unstableMask = (abs(eigenValues) >= 1);
    stage1Mask   = stage1Mask & ~unstableMask;
end

if opts.filterNegativeDamp
    negDampMask = (dampingFactors < 0);
    stage1Mask  = stage1Mask & ~negDampMask;
end

if opts.dampingThreshold
    overdampedMask = (dampingFactors > opts.dampingThreshold);
    stage1Mask     = stage1Mask & ~overdampedMask;
end

zeroFreqMask = (frequencies <= 0);
stage1Mask   = stage1Mask & ~zeroFreqMask;

%% Apply Stage 1 filter to matrices — preserve conjugate pairs
filteredAt = At(stage1Mask, stage1Mask);
filteredBt = Bt(stage1Mask);
filteredCt = Ct(stage1Mask);

%% Stage 2 — Conjugate pair reduction (applied only to parameter outputs)
stage2Mask = stage1Mask;
negFreqMask = false(size(eigenValues));

if opts.filterNegativeFreq
    negFreqMask = (imag(eigenValues) < 0);
    stage2Mask  = stage2Mask & ~negFreqMask;
end

%% Apply Stage 2 filter to parameter outputs only
filteredFrequencies    = frequencies(stage2Mask);
filteredDecayRates     = decayRates(stage2Mask);
filteredDampingFactors = dampingFactors(stage2Mask);
filteredModes          = modes(stage2Mask);

%% Optional debug info.
if nargout > 7 && isfield(opts, 'returnDebug') && opts.returnDebug

    diagnostics.At                     = At;
    diagnostics.Bt                     = Bt;
    diagnostics.Ct                     = Ct;
    diagnostics.totalEigenvalues       = numel(eigenValues);
    diagnostics.retainedCount          = sum(stage2Mask);
    diagnostics.filteredCount          = sum(~stage2Mask);
    diagnostics.filteredOutEigenvalues = eigenValues(~stage2Mask);
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