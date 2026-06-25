function [outputSignal, physicalModes, varargout] = ...
    buildModalImpulseResponse(At, Bt, Ct, D, nSamples, opts)
% BUILDMODALIMPULSERESPONSE  Reconstructs impulse response and individual
%                            modal contributions.
%
% DESCRIPTION:
%   Computes the impulse response of an using the diagonalised modal state 
%   space matrices. Uses the diagonal structure of At:
%
%       x_t(k) = At^(k-1) * Bt
%       y(k)   = Ct * x_t(k)
%
%   The direct feedthrough term D is prepended to the output signal.
%   A zero column is similarly prepended to physical modes to align
%   sample indexing.
%
%   IMPORTANT — Caller responsibilities:
%   (1) At, Bt, Ct must be filtered to contain only physical modes,
%       with conjugate pairs preserved. Filtering is performed in
%       extractModalParameters().
%   (2) Conjugate pairs must be stored at adjacent indices. MATLAB's
%       eig() satisfies this for real input matrices by default.
%
% INPUTS:
%   At       - [double, m x m]   Diagonal modal eigenvalue matrix from
%              extractModalParameters(). m must be even.
%   Bt       - [double, m x 1]   Modal input vector.
%   Ct       - [double, 1 x m]   Modal output vector.
%   D        - [double, scalar]  Direct feedthrough term y(0).
%   nSamples - [integer, scalar] Number of output samples.
%   opts     - [struct]          Configuration struct.
%
% OPTIONS (opts fields consumed by this function):
%   returnDiagnostics - [logical] Return diagnostics struct as 3rd
%                       output. [false] Controlled by caller for debug
%                       collection.
%
% OUTPUTS:
%   outputSignal  - [double, 1 x nSamples] Reconstructed impulse response
%                   including D at sample 1.
%   physicalModes - [double, (m/2) x nSamples] Individual physical mode
%                   trajectories.
%   varargout{1}  - [struct] Diagnostics struct (if opts.returnDiagnostics
%                   and nargout > 2).
%
% NOTES:
%   The conjugate pair recombination assumes pairs are stored at adjacent
%   indices.
%
%   Diagnostics struct fields (when returned):
%     .modalCoordinates       - Modal coordinate trajectories before
%                               output projection [m x nSamples]
%     .modeEnergies           - Energy per physical mode [(m/2) x 1]
%     .modePeakAmplitudes     - Peak absolute amplitude per mode
%     .totalEnergy            - Sum of mode energies (parseval check)
%     .reconstructionEnergy   - Energy of full reconstructed signal
%
% REFERENCES:
%   [1] Juang, J.-N. (1994). Applied System Identification.
%       Prentice Hall, Englewood Cliffs, NJ. Ch. 5-7.
%
% AUTHOR:   -
% DATE:     June 2026
% VERSION:  1.1
%
% See also: eigensystemRealisation, extractModalParameters,
%           buildImpulseResponse
%% Lightweight contract assertions
assert(isequal(size(At, 1), size(At, 2)), ...
    'buildModalImpulseResponse:nonSquareAt', 'At must be square.');
assert(mod(size(At, 1), 2) == 0, ...
    'buildModalImpulseResponse:oddDimension', ...
    ['At dimension must be even — conjugate pairs expected. ' ...
     'Got dimension %d.'], size(At, 1));
assert(size(Bt, 1) == size(At, 1), ...
    'buildModalImpulseResponse:BtDimensionMismatch', ...
    'Bt rows must match At dimension.');
assert(size(Ct, 2) == size(At, 1), ...
    'buildModalImpulseResponse:CtDimensionMismatch', ...
    'Ct columns must match At dimension.');

%% Fixed values.

% Used to index every other mode as they are conjugate pairs.
PHYSICAL_MODE_STRIDE = 2;

% exp(wt) = 2*cos(wt).
AMPLTIUDE_SCALING_FACTOR = 2;

% Select the negative frequency components of the signal so the output
% modes are sinusoids instead of cosines.
PHYSICAL_MODE_INDEX = 2;

%% Main code.
% At is the diagonalised form of the state space matrices. Transformation
% into a vector allows is more efficient.
At = diag(At);

% The relation for x_t[k] = At ^ (k - 1) * Bt. A vector of powers for the
% calculation. nSamples must be reduced by 2 as indexing starts from 0, and
% there must be a sample reserved for the feedthrough term D {y(0)}.
nPowers = 0 : nSamples - 2;
eigenValuePowers = At .^ nPowers;

modalCoordinates = eigenValuePowers .* Bt;

outputSignal = real(Ct * modalCoordinates);

physicalModes = real(Ct.' .* modalCoordinates);

% Prepend the feedthrough term to the start of the signal and prepend 0 to
% the physical modes.
outputSignal = [D outputSignal];

physicalModes = [zeros(size(physicalModes, 1), 1) physicalModes];

% Odd indexed modes correspond to the cosine components of the modes, and
% even indexed modes correspond to the sin components. Selecting the
% sinusoids this nicer to plot, though there isn't really any difference.
% The amplitude is then scaled by 2 to compensate for losing the 
% conjugate pair.
% Note: Assumes conjugate pairs are adjacent and the loss of floating point
% precision is minimal from discarding the conjugate pair.
physicalModes = physicalModes(PHYSICAL_MODE_INDEX : ...
    PHYSICAL_MODE_STRIDE : end, :) * AMPLTIUDE_SCALING_FACTOR;

%% Optional diagnostics output
if nargout > 2 && isfield(opts, 'returnDebug') && opts.returnDebug

    diagnostics.modalCoordinates     = modalCoordinates;

    % Energy and amplitude metrics per physical mode
    diagnostics.modeEnergies       = sum(physicalModes .^ 2, 2);
    diagnostics.modePeakAmplitudes = max(abs(physicalModes), [], 2);
    diagnostics.totalEnergy        = sum(diagnostics.modeEnergies);
    diagnostics.reconstructionEnergy = sum(outputSignal .^ 2);

    varargout{1} = diagnostics;

end

end