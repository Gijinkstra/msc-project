function [hankelMatrix, shiftedHankelMatrix] = ...
    constructHankelMatrices(impulseResponse, nRows, nColumns)
% CONSTRUCTHANKELMATRICES  Constructs the Hankel and shifted Hankel
%                          matrices for the ERA pipeline.
%
% DESCRIPTION:
%   Constructs two Hankel matrices from a pre-sliced impulse response
%   signal vector:
%
%       H(0) : built from impulseResponse(1), impulseResponse(2), ...
%       H(1) : built from impulseResponse(2), impulseResponse(3), ...
%
%   H(0) and H(1) are related by a one-sample shift
%
%   The Hankel matrix is constructed using MATLAB's built-in hankel()
%   function, which requires the first column vector and the last row
%   vector as inputs. Refer to MATLAB documentation for detail.
%
% INPUTS:
%   impulseResponse - [double, N x 1]    Impulse response of a SISO system.
%   nRows           - [integer, scalar]  Number of Hankel matrix rows.
%                     Set by opts.nRows in eigensystemRealisation().
%   nColumns        - [integer, scalar]  Number of Hankel matrix columns.
%                     Set by opts.nColumns in eigensystemRealisation().
%   opts            - [struct]           Configuration struct.
%
% OPTIONS (opts fields consumed by this function):
%   None.
%
% OUTPUTS:
%   hankelMatrix        - [double, nRows x nColumns]  H[0] SISO Hankel 
%                         matrix.
%   shiftedHankelMatrix - [double, nRows x nColumns]  H(1). Shifted
%                         Hankel matrix  offset by one sample from H(0).
%% Input check.

assert(numel(impulseResponse) >= nRows + nColumns, ...
    'constructHankelMatrices:insufficientSamples', ...
    'Signal too short — pipeline validation failed to catch this.');

%% Main code.

% Build the generalised Hankel Matrix, Eq. (18) [1].
hankelMatrix = buildHankelMatrix(impulseResponse, nRows, nColumns);

% Build the shifted Hankel Matrix, Eq. (19) [1].
shiftedHankelSignal = impulseResponse(2 : end);
shiftedHankelMatrix = buildHankelMatrix(shiftedHankelSignal, nRows, ...
    nColumns);

end

function hankelMatrix = buildHankelMatrix(impulseResponse, nRows, ...
    nColumns)

% Construct Hankel matrix. Hankel matrix input requires the first column
% and last row. Consult MATLAB documentation for more detail.
hankelColumn = impulseResponse(1 : nRows);
hankelRow = impulseResponse(nRows : nRows + nColumns - 1);
hankelMatrix = hankel(hankelColumn, hankelRow);

end