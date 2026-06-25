function [Um, Sm, Vm, varargout] = svdTruncation(hankelMatrix, ...
    svdTolerance, opts)
% SVDTRUNCATION  Computes the truncated singular value decomposition
%                of the Hankel matrix for ERA modal subspace identification.
%
% DESCRIPTION:
%   Performs the singular value decomposition of the Hankel matrix and
%   returns the truncated factors corresponding to the signal subspace.
%   The truncation order is determined by the svdTolerance parameter,
%   which supports two complementary modes:
%
%   MODE 1 — Threshold mode (svdTolerance < 1):
%       Retains all normalised singular values above the threshold:
%           sigma_i / sigma_1 > svdTolerance
%       Useful when the number of modes is unknown but a noise floor
%       can be estimated. An even number of modes will be retained in the
%       event the noise floor is off (exclusive).
%
%   MODE 2 — Count mode (svdTolerance >= 1):
%       Retains exactly 2 * svdTolerance singular values, treating
%       svdTolerance as the number of physical modes expected. Each
%       physical mode requires two singular values for its conjugate
%       pair, so the truncation order is doubled.
%
%   In both modes the final truncation order is clamped to the smaller
%   of the requested count and the minimum dimension of S, ensuring
%   the truncation does not exceed available singular values.
%
% INPUTS:
%   hankelMatrix - [double, M x N]   Hankel matrix from
%                  constructHankelMatrices().
%   svdTolerance - [double, scalar]  Truncation control parameter.
%                  Two interpretations depending on magnitude:
%                  - If < 1: normalised threshold (e.g. 0.01 = 1%).
%                  - If >= 1: number of physical modes (doubled
%                    internally for conjugate pairs).
%                  Set by opts.svdTolerance in
%                  eigensystemRealisation().
%   opts         - [struct]          Configuration struct.
%
% OPTIONS (opts fields consumed by this function):
%   validateInputs    - [logical]  Run input validation. [true]
%   returnDiagnostics - [logical] Return diagnostics struct as 4th
%                       output. [false] Controlled by caller for debug
%                       collection.
%
% OUTPUTS:
%   Um           - [double, M x n]  Truncated left singular vectors.
%   Sm           - [double, n x n]  Truncated diagonal singular value
%                  matrix.
%   Vm           - [double, N x n]  Truncated right singular vectors.
%   varargout{1} - [struct]         Diagnostics struct.
%
% NOTES:
%   Diagnostics struct fields (when returned):
%     .singularValues     - Full singular value vector before truncation
%     .truncationOrder    - Final n after doubling and clamping
%     .toleranceMode      - 'threshold' or 'count'
%     .reconstructionErr  - Frobenius norm relative reconstruction error
%     .conditionNumber    - Condition number of Sm (sigma_1 / sigma_n)
%     .gapRatio           - sigma_n / sigma_(n+1) at truncation boundary
%                           (NaN if n equals full rank)
%     .effectiveRank      - Numerical rank from rank(hankelMatrix)
%     .requestedModes     - User-requested mode count before clamping
%     .clampedToBoundary  - True if truncation was clamped to S size
%
% REFERENCES:
%   [1] Juang, J.-N. & Pappa, R.S. (1985). An eigensystem realization
%       algorithm for modal parameter identification and model reduction.
%       Journal of Guidance, Control, and Dynamics, 8(5), 620-627.
%   [2] Juang, J.-N. (1994). Applied System Identification.
%       Prentice Hall, Englewood Cliffs, NJ. Ch. 5.
%   [3] pyyeti ERA implementation (reference):
%       https://github.com/twmacro/pyyeti/blob/master/pyyeti/era.py
%
% AUTHOR:   -
% DATE:     June 2026
% VERSION:  1.1
%
% See also: eigensystemRealisation, constructHankelMatrices,
%           stateSpaceMatrixConstruction
%% Input check.

assert(~isempty(hankelMatrix), ...
    'svdTruncation:emptyHankel', 'hankelMatrix must not be empty.');

assert(svdTolerance > 0, ...
    'svdTruncation:invalidTolerance', 'svdTolerance must be positive.');

%% Fixed parameters
% The cutoff value for svdTolerance to change mode truncation method.
SVD_TOLERANCE_CUTOFF = 1;

% Scaling factor applied to the svdTolerance when >1.
MODE_SCALING_FACTOR = 2;

% This is here because I'm pedantic about magic numbers.
MAX_SV_INDEX = 1;

SVD_TOLERANCE_FLAG = {'threshold', 'count'};

%% Main code.
[U, S, V] = svd(hankelMatrix);

% Minimum dimension for the singular value matrix. In the event S is not
% square.
minSingularValueIndex = min(size(S, 1), size(S, 2));

% Avoid calculations on non-diagonal elements.
sigma = diag(S);

if svdTolerance < SVD_TOLERANCE_CUTOFF

    % Internal method flag.
    toleranceMode = SVD_TOLERANCE_FLAG{1};
    
    % SVD always returns the singular value with the maximum energy in the
    % first element.
    indexToKeep = find(sigma ./ sigma(MAX_SV_INDEX) > svdTolerance);

    % Internal error check.
    assert(~isempty(indexToKeep), ...
        'svdTruncation:noModesAboveThreshold', ...
        ['No singular values exceed the normalised threshold of %.4g. ' ...
        'Reduce svdTolerance or check Hankel matrix conditioning.'], ...
        svdTolerance);

    nModesRequested = indexToKeep(end);

else
    
    toleranceMode = SVD_TOLERANCE_FLAG{2};
    nModesRequested = svdTolerance;

end

if strcmp(toleranceMode, SVD_TOLERANCE_FLAG{2})
    nModesRequested = nModesRequested * MODE_SCALING_FACTOR;
end

% Ensure the number of requested modes has a smaller dimension than the
% smallest dimension of the singular value matrix.
nModes = min(nModesRequested, minSingularValueIndex);
clampedToBoundary = nModes > minSingularValueIndex;

% Truncation.
Um = U(:, 1 : nModes);
Sm = S(1 : nModes, 1 : nModes);
Vm = V(:, 1 : nModes);

if nargout > 3 && isfield(opts, 'returnDebug') && opts.returnDebug

    diagnostics.singularValues    = sigma;
    diagnostics.truncationOrder   = nModes;
    diagnostics.toleranceMode     = toleranceMode;
    diagnostics.requestedModes    = nModesRequested;
    diagnostics.clampedToBoundary = clampedToBoundary;
    diagnostics.reconstructionErr = ...
        norm(hankelMatrix - Um * Sm * Vm', 'fro') / ...
        norm(hankelMatrix, 'fro');
    diagnostics.conditionNumber = sigma(1) / sigma(nModes);

    if nModes < length(sigma)
        diagnostics.gapRatio = sigma(nModes) / sigma(nModes + 1);
    else
        diagnostics.gapRatio = NaN;
    end

    diagnostics.effectiveRank = rank(hankelMatrix);
    varargout{1}              = diagnostics;

end

end