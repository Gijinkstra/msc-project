function [outputSignal, signalStates] = buildImpulseResponse(A, B, C, D, ...
    nSamples, opts)
% BUILDIMPULSERESPONSE  Computes the impulse response of an ERA state space
%                       system via direct state recursion.
%
% DESCRIPTION:
%   Simulates the discrete time state space system:
%
%       x(k+1) = A * x(k) + B * u(k)
%       y(k)   = C * x(k) + D * u(k)
%
%   from zero initial conditions with a unit impulse input u(0) = 1,
%   u(k) = 0 for k > 0. This recovers the Markov parameter sequence
%   y(k) = C * A^(k-1) * B for k >= 1, y(0) = D, which is the impulse
%   response the ERA Hankel matrix was constructed from.
%
%   The loop is retained here deliberately — the state recursion is
%   inherently sequential (x(k+1) depends on x(k)) and cannot be
%   parallelised without transforming to modal coordinates. For the
%   modal superposition alternative see modalSuperposition.m.
%
% INPUTS:
%   A        - [double, n x n]  State transition matrix (from ERA).
%   B        - [double, n x 1]  Input matrix (from ERA).
%   C        - [double, 1 x n]  Output matrix, row vector (from ERA).
%   D        - [double, scalar] Direct feedthrough term (from ERA).
%   nSamples - [integer, scalar] Number of impulse response samples
%              to compute.
%   opts     - [struct]          Optional configuration (see below).
%
% OPTIONS (opts fields):
%   validateInputs - Run input validation checks.          [true]
%   verbose        - Print diagnostics to console.         [false]
%   returnStates   - Return full state matrix.             [true]
%                    Set false to save memory for large n or nSamples.
%
% OUTPUTS:
%   outputSignal - [double, 1 x nSamples]     Impulse response sequence
%                  y(0), y(1), ..., y(nSamples-1).
%   signalStates - [double, n x nSamples]     State trajectory matrix.
%                  Column k contains the state vector x(k).
%                  Returns empty matrix if opts.returnStates = false.
%
% NOTES:
%   Zero initial conditions are enforced by preallocation with zeros().
%   The unit impulse is applied at sample k=0 only (inputSignal(1) = 1).
%
%   The final output sample y(nSamples) requires a separate evaluation
%   outside the loop since the loop runs from k=1 to k=nSamples-1,
%   updating states one step ahead. The state at k=nSamples is available
%   after the loop completes.
%
%   For large systems or long simulations, set opts.returnStates = false
%   to avoid allocating the n x nSamples state matrix.
%
% ASSUMPTIONS & LIMITATIONS:
%   - Zero initial conditions (required by ERA formulation).
%   - SISO system (scalar D, 1 x n C, n x 1 B).
%   - Unit impulse input. For arbitrary inputs use modalSuperposition.m.
%   - A must represent a stable system (all eigenvalues inside unit
%     circle) to produce a decaying impulse response.
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
% VERSION:  1.0
%
% See also: eigensystemRealisation, stateSpaceMatrixConstruction,
%           extractModalParameters, modalSuperposition
%% Default options
defaults.validateInputs = true;
defaults.verbose        = false;
defaults.returnStates   = true;

if nargin < 6 || isempty(opts)

    opts = defaults;

else

    fields = fieldnames(defaults);

    for k = 1:numel(fields)

        if ~isfield(opts, fields{k})

            opts.(fields{k}) = defaults.(fields{k});

        end

    end

end

%% Input validation
if opts.validateInputs

    validateImpulseInputs(A, B, C, D, nSamples);

end

%% Preallocation
nStates = size(A, 1);

if opts.returnStates

    signalStates = zeros(nStates, nSamples);

else

    signalStates = [];                           % suppress state storage

end

outputSignal = zeros(1, nSamples);

% Unit impulse function.
inputSignal = zeros(nSamples, 1);
inputSignal(1) = 1;

thisState = zeros(nStates, 1);

%% Main code.

for thisSample = 1 : nSamples - 1

    % Current input.
    thisInput = inputSignal(thisSample);

    % Output at current step.
    outputSignal(thisSample) = C * thisState + D * thisInput;

    % State update — propagate to next step.
    nextState = A * thisState + B * thisInput;

    % Store state if requested.
    if opts.returnStates
        signalStates(:, thisSample + 1) = nextState;
    end

    thisState = nextState;
end

% Last sample of the output signal.
outputSignal(nSamples) = C * thisState + D * inputSignal(nSamples);

%% Verbose output
if opts.verbose
    poles        = eig(A);
    stablePoles  = sum(abs(poles) < 1);
    peakResponse = max(abs(outputSignal));
    finalValue   = abs(outputSignal(end));

    fprintf('--- Build Impulse Response ---\n');
    fprintf('State dimension  : %d\n',   n);
    fprintf('Samples computed : %d\n',   nSamples);
    fprintf('Stable poles     : %d / %d\n', stablePoles, n);
    fprintf('Peak response    : %.6f\n', peakResponse);
    fprintf('Final value      : %.6f\n', finalValue);

    if finalValue > 0.01 * peakResponse
        warning('buildImpulseResponse:insufficientDecay', ...
            ['Final sample magnitude (%.6f) is %.1f%% of peak (%.6f). ' ...
             'Impulse response may not have decayed sufficiently. ' ...
             'Consider increasing nSamples.'], ...
            finalValue, 100 * finalValue / peakResponse, peakResponse);
    end

    fprintf('------------------------------\n');
end

end

%% Validation subfunction
function validateImpulseInputs(A, B, C, D, nSamples)
% VALIDATEIMPULSEINPUTS  Validates inputs to buildImpulseResponse.

% --- A ---
validateattributes(A, {'double'}, ...
    {'2d', 'square', 'nonempty', 'finite'}, ...
    'buildImpulseResponse', 'A');

% --- B ---
validateattributes(B, {'double'}, ...
    {'column', 'nonempty', 'finite'}, ...
    'buildImpulseResponse', 'B');

% --- C ---
validateattributes(C, {'double'}, ...
    {'row', 'nonempty', 'finite'}, ...
    'buildImpulseResponse', 'C');

% --- D ---
validateattributes(D, {'double'}, ...
    {'scalar', 'real', 'finite'}, ...
    'buildImpulseResponse', 'D');

% --- nSamples ---
validateattributes(nSamples, {'numeric'}, ...
    {'scalar', 'integer', 'positive'}, ...
    'buildImpulseResponse', 'nSamples');

%% Cross-parameter checks
n = size(A, 1);

assert(size(B, 1) == n, ...
    'buildImpulseResponse:dimensionMismatch', ...
    'B has %d rows but A is %d x %d. B rows must equal state dimension.', ...
    size(B, 1), n, n);

assert(size(C, 2) == n, ...
    'buildImpulseResponse:dimensionMismatch', ...
    'C has %d columns but A is %d x %d. C columns must equal state dimension.', ...
    size(C, 2), n, n);

% Stability check — warn rather than error, caller may intentionally
% simulate an unstable system over a short window
poles = eig(A);
if any(abs(poles) >= 1)
    warning('buildImpulseResponse:unstableSystem', ...
        ['A has %d pole(s) on or outside the unit circle. ' ...
         'Impulse response will not decay. ' ...
         'Check ERA truncation order and Hankel dimensions.'], ...
        sum(abs(poles) >= 1));
end

end