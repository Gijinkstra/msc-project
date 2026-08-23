function [outputSignal, frequencies, dampingFactors, amplitudes, ...
    phases, modes, varargout] = eigensystemRealisation(impulseResponse, ...
    svdTolerance, fs, opts)
% EIGENSYSTEMREALISATION  Identifies modal parameters via the Eigensystem
%                         Realisation Algorithm.
%
% DESCRIPTION:
%   Constructs Hankel matrices from a measured impulse response, performs
%   truncated singular value decomposition, recovers a discrete state
%   space realisation, and extracts modal frequencies, damping factors,
%   decay rates, and individual modal impulse responses.
%
%   The pipeline:
%     1. Slice impulse response from opts.hankelSignalIndex onwards.
%     2. Construct H(0) and H(1) Hankel matrices.
%     3. SVD truncation per svdTolerance.
%     4. Recover (A, B, C, D) via the ERA realisation formulae.
%     5. Extract modal parameters and filter non-physical modes.
%     6. Reconstruct impulse response via modal superposition.
%
%   All user input validation is performed once at the entry point. All
%   subfunctions trust their inputs and perform only lightweight contract
%   assertions.
%
% INPUTS:
%   impulseResponse - [double, N x 1]    Measured impulse response samples.
%   svdTolerance    - [double, scalar]   SVD truncation control.
%                     If < 1: normalised threshold (e.g. 0.01 = 1% of
%                             largest singular value).
%                     If >= 1: number of physical modes (doubled
%                              internally for conjugate pairs).
%   fs              - [double, scalar]   Sample rate in Hz.
%   opts            - [struct, optional] Configuration struct from
%                     eraConfig(). If omitted, internal defaults are used.
%
% OPTIONS (opts fields consumed across the pipeline):
%   Pipeline control:
%     validateInputs        - Run input validation.           [true]
%     returnDebug           - Return debug struct.            [false]
%
%   Hankel construction:
%     hankelSignalIndex     - Start index for Hankel signal.  [2]
%     nRows                 - Hankel row count. []=auto-size. [[]]
%     nColumns              - Hankel column count. []=auto.   [[]]
%
%   Modal parameter extraction:
%     filterUnstable        - Remove unstable poles.          [true]
%     filterNegativeDamp    - Remove negative damping modes.  [true]
%     filterNegativeFreq    - Remove negative frequency.      [true]
%     filterOverdamped      - Remove overdamped modes.        [true]
%     dampingThreshold      - Upper damping bound.            [0.99]
%
% OUTPUTS:
%   outputSignal    - [double, 1 x N]  Reconstructed impulse response.
%   frequencies     - [double, k x 1]  Natural frequencies in Hz.
%   dampingFactors  - [double, k x 1]  Modal damping ratios.
%   modalImpulses   - [double, k x N]  Per-mode impulse contributions.
%   varargout{1}    - [struct]         Debug struct (if opts.returnDebug
%                                      and nargout > 4).
%
% USAGE:
%   % Minimal — uses internal defaults
%   [sig, f, z, phi] = eigensystemRealisation(ir, 0.01, fs);
%
%   % Mode count mode
%   [sig, f, z, phi] = eigensystemRealisation(ir, 5, fs);
%
%   % With explicit configuration
%   opts = eraConfig();
%   opts.returnDebug = true;
%   [sig, f, z, phi, dbg] = eigensystemRealisation(ir, 0.01, fs, opts);
%
% REFERENCES:
%   [1] Juang, J.-N. & Pappa, R.S. (1985). An eigensystem realization
%       algorithm for modal parameter identification and model reduction.
%       Journal of Guidance, Control, and Dynamics, 8(5), 620-627.
%   [2] Juang, J.-N. (1994). Applied System Identification.
%       Prentice Hall, Englewood Cliffs, NJ.
%
% AUTHOR:   -
% DATE:     June 2026
% VERSION:  2.0
%
% See also: eraConfig, constructHankelMatrices, svdTruncation,
%           stateSpaceMatrixConstruction, extractModalParameters,
%           buildModalImpulseResponse

%% Default configuration
defaults.validateInputs     = true;
defaults.returnDebug        = false;

defaults.hankelSignalIndex  = 2;
defaults.nRows              = [];
defaults.nColumns           = [];

defaults.residueScaling        = false;
defaults.filterUnstable     = true;
defaults.filterNegativeDamp = true;
defaults.filterNegativeFreq = true;
defaults.filterOverdamped   = true;
defaults.dampingThreshold   = 0.99;

%% Input checks and optional arguments.
maxArgs = 4;
minArgs = 3;
narginchk(minArgs, maxArgs);

if nargin < maxArgs || isempty(opts)

    opts = defaults;

else

    fields = fieldnames(defaults);

    for k = 1:numel(fields)

        if ~isfield(opts, fields{k})

            opts.(fields{k}) = defaults.(fields{k});

        end

    end

end

%% Hankel matrix input validation.
nSamples = numel(impulseResponse);

% Defaults.
nRows = opts.nRows;
nColumns = opts.nColumns;

if isempty(svdTolerance)
    svdTolerance = opt.svdTolerance;
end

% If the row argument is empty, maxmimise the coverage of the Hankel
% matrices. 1 offset to account for the shifted Hankel matrix.
if isempty(nRows)
    nRows = floor(nSamples / 2) - 1;
end

if isempty(nColumns)
    nColumns = floor(nSamples / 2) - 1;
end

if opts.validateInputs
    validateERAPipelineInputs(impulseResponse, svdTolerance, fs, ...
        nRows, nColumns, opts);
end

%% Constants.
timeStep = 1 / fs;
firstSample = impulseResponse(1);
hankelIndex = opts.hankelSignalIndex;

%% Main code.

%%% Step 1 - Hankel Matrices
% The Hankel matrix construction starts from the second sample as y(0) = D.
hankelMatrixSignal = impulseResponse(hankelIndex : end);

[hankelMatrix, shiftedHankelMatrix] = ...
    constructHankelMatrices(hankelMatrixSignal, nRows, nColumns);

%%% Step 2 - SVD with truncation.
[Um, Sm, Vm, svdDebug] = svdTruncation(hankelMatrix, svdTolerance, opts);

%%% State space matrix construction
[A, B, C, D, ssmDebug] = stateSpaceMatrixConstruction(Um, Sm, Vm, ...
    shiftedHankelMatrix, nRows, nColumns, firstSample, opts);

%%% Recover frequencies, damping factor and modes.
[frequencies, dampingFactors, amplitudes, phases, modes, decayRates, ...
    At, Bt, Ct, modalDebug] = extractModalParameters(A, B, C, timeStep, ...
    opts);

%%% Single pole states and output signal.
[outputSignal, modalImpulses, impulseDebug] = ...
    buildModalImpulseResponse(At, Bt, Ct, D, nSamples, opts);

%% Optional debug output
if nargout > 4 && opts.returnDebug

    debug = struct();

    % Hankel stage
    debug.Hankel.hankelMatrix            = hankelMatrix;
    debug.Hankel.shiftedHankelMatrix     = shiftedHankelMatrix;
    debug.Hankel.nRows                   = nRows;
    debug.Hankel.nColumns                = nColumns;

    % SVD stage — from svdTruncation diagnostics
    debug.SVD                            = svdDebug;
    debug.SVD.Um                         = Um;
    debug.SVD.Sm                         = Sm;
    debug.SVD.Vm                         = Vm;

    % State space stage — from stateSpaceMatrixConstruction diagnostics
    debug.SSM                            = ssmDebug;
    debug.SSM.A                          = A;
    debug.SSM.B                          = B;
    debug.SSM.C                          = C;
    debug.SSM.D                          = D;

    % Modal stage — from extractModalParameters diagnostics
    debug.Modal                          = modalDebug;
    debug.Modal.At                       = At;
    debug.Modal.Bt                       = Bt;
    debug.Modal.Ct                       = Ct;
    debug.Modal.decayRates               = decayRates;
    debug.Modal.modes                    = modes;

    % Reconstruction stage — quality metrics
    debug.Reconstruction                 = impulseDebug;
    debug.Reconstruction.outputSignal    = outputSignal;
    debug.Reconstruction.modalImpulses   = modalImpulses;
    debug.Reconstruction.impulseNMSE     = sum((impulseResponse(:) - ...
        outputSignal(:)) .^ 2) / sum(impulseResponse(:) .^ 2);
    debug.Reconstruction.relativeError   = ...
        norm(impulseResponse(:) - outputSignal(:)) / ...
        norm(impulseResponse(:));

    varargout{1} = debug;

end

end