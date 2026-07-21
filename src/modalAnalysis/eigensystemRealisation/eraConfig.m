function opts = eraConfig()
% ERACONFIG  Top-level configuration for the ERA pipeline.
%
% DESCRIPTION:
%   Returns a struct containing all configuration options for the
%   Eigensystem Realisation Algorithm pipeline. Edit values here to
%   change pipeline behaviour globally. Pass the returned struct to
%   eigensystemRealisation and any subfunction directly.
%
%   All subfunctions accept this struct and extract only the fields
%   they require, so the full opts struct can be passed anywhere in
%   the pipeline without modification.
%
% USAGE:
%   opts = eraConfig();
%   [sig, f, z, phi] = eigensystemRealisation(ir, nModes, fs, opts);
%
% AUTHOR:   -
% DATE:     June 2026
% VERSION:  1.0
%
% See also: eigensystemRealisation, stateSpaceMatrixConstruction,
%           extractModalParameters, buildImpulseResponse

%% -----------------------------------------------------------------------
%  PIPELINE CONTROL
%  -----------------------------------------------------------------------
% Enable or disable input validation across all pipeline stages.
% Recommended: true during development, may be set false for production
% batch runs where inputs are well-controlled.
opts.validateInputs       = true;

% Return debug struct as optional output from eigensystemRealisation.
% Requires nargout > 4 at the call site.
opts.returnDebug          = false;

%% -----------------------------------------------------------------------
%  HANKEL MATRIX
%  -----------------------------------------------------------------------

% Start index into the impulse response for Hankel matrix construction.
% Index 2 skips y(0), which is assigned to D (direct feedthrough).
% Do not change unless the ERA signal indexing convention changes.
opts.hankelSignalIndex    = 2;

opts.nRows = [];
opts.nColumns = [];

%% -----------------------------------------------------------------------
%  SVD TRUNCATION
%  -----------------------------------------------------------------------
% The tolerance level for the singular values retained in the SVD of the
% Hankel matrix. Setting a threshold value < 1 will retain all normalised
% singular values above the threshold. Setting a value greater than one
% will retain the number of modes specified (rounded up to the next even
% number).
opts.svdTolerance = 0.01;

%% -----------------------------------------------------------------------
%  MODAL PARAMETER EXTRACTION
%  -----------------------------------------------------------------------
% Scale the residues by a factor of the poles.
opts.residueScaling       = true;
% Remove poles on or outside the unit circle (unstable modes).
opts.filterUnstable       = true;

% Remove modes with negative damping ratios (non-physical).
opts.filterNegativeDamp   = true;

% Remove modes with negative frequency (conjugate pair reduction).
% Retains only the positive imaginary half of each conjugate pair.
opts.filterNegativeFreq   = true;

opts.filterOverdamped     = true;

% Upper damping ratio bound for physical mode retention.
% Modes with damping ratio above this threshold are considered
% overdamped and non-resonant. Adjust for systems with high damping.
% Recommended range: 0.99 to 1.0.
opts.dampingThreshold     = 0.99;

%% -----------------------------------------------------------------------
%  VALIDATION
%  -----------------------------------------------------------------------
validateConfig(opts);

end

%% Validation subfunction
function validateConfig(opts)
% VALIDATECONFIG  Checks that all eraConfig fields are self-consistent.

% --- Pipeline control ---
validateattributes(opts.validateInputs, {'logical', 'numeric'}, ...
    {'scalar', 'binary'}, 'eraConfig', 'opts.validateInputs');

validateattributes(opts.returnDebug, {'logical', 'numeric'}, ...
    {'scalar', 'binary'}, 'eraConfig', 'opts.returnDebug');

% --- Hankel matrix ---
validateattributes(opts.hankelSignalIndex, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, ...
    'eraConfig', 'opts.hankelSignalIndex');

% --- SVD truncation ---
validateattributes(opts.svdTolerance, {'numeric'}, ...
    {'scalar', 'positive'}, 'eraConfig', 'opts.svdTolerance');

% --- Modal parameter extraction ---
validateattributes(opts.filterUnstable, {'logical', 'numeric'}, ...
    {'scalar', 'binary'}, 'eraConfig', 'opts.filterUnstable');

validateattributes(opts.filterNegativeDamp, {'logical', 'numeric'}, ...
    {'scalar', 'binary'}, 'eraConfig', 'opts.filterNegativeDamp');

validateattributes(opts.filterNegativeFreq, {'logical', 'numeric'}, ...
    {'scalar', 'binary'}, 'eraConfig', 'opts.filterNegativeFreq');

validateattributes(opts.dampingThreshold, {'double'}, ...
    {'scalar', 'real', 'positive', 'finite', '<=', 1}, ...
    'eraConfig', 'opts.dampingThreshold');

end