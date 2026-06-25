function [A, B, C, D, varargout] = stateSpaceMatrixConstruction(U, S, V, ...
    shiftedHankelMatrix, nRows, nColumns, firstSample, opts)
% STATESPACEMATRIXCONSTRUCTION  Constructs ERA state space matrices from
%                               truncated SVD factors.
%
% DESCRIPTION:
%   Recovers the discrete time state space matrices (A, B, C, D) from the
%   truncated SVD of the Hankel matrix and the shifted Hankel matrix.
%
%       A = S^(-1/2) * U' * H(1) * V * S^(-1/2)
%       B = S^(1/2)  * V' * e_c
%       C = S^(1/2)  * U' * e_r
%       D = y(0)
%
%   where:
%       e_c is the first canonical basis vector of length nColumns,
%           selecting the first block column of the controllability matrix.
%       e_r is the first canonical basis vector of length nRows,
%           selecting the first block row of the observability matrix.
%
% INPUTS:
%   U                   - [double, nRows x n]    Truncated left singular
%                         vectors.
%   S                   - [double, n x n]        Truncated diagonal
%                         singular value matrix.
%   V                   - [double, nColumns x n] Truncated right singular
%                         vectors.
%   shiftedHankelMatrix - [double, nRows x nColumns]
%                         Hankel matrix shifted by one time step, H(1).
%   nRows               - [integer, scalar]  Number of Hankel rows.
%                         Sets length of e_r basis vector.
%                         Set by opts.nRows in eigensystemRealisation().
%   nColumns            - [integer, scalar]  Number of Hankel columns.
%                         Sets length of e_c basis vector.
%                         Set by opts.nColumns in eigensystemRealisation().
%   firstSample         - [double, scalar]   First sample of the impulse
%                         response y(0). Assigned to direct feedthrough D.
%                         In MATLAB indexing this is impulseResponse(1).
%   opts                - [struct]           Configuration struct.
%
% OPTIONS (opts fields consumed by this function):
%   returnDiagnostics - [logical] Return diagnostics struct as 5th
%                       output. [false] Controlled by caller for debug
%                       collection.
%
% OUTPUTS:
%   A            - [double, n x n]  State transition matrix.
%   B            - [double, n x 1]  Input matrix (SISO).
%   C            - [double, 1 x n]  Output matrix (row vector, SISO).
%   D            - [double, scalar] Direct feedthrough term.
%   varargout{1} - [struct]         Diagnostics struct (if requested via
%                                   opts.returnDiagnostics and nargout > 4).
%                                   See NOTES for fields.
%
% NOTES:
%   Diagnostics struct fields (when returned):
%     .stateDimension      - Dimension of A (n)
%     .conditionNumberA    - Condition number of identified A
%     .BNorm               - Frobenius norm of B
%     .CNorm               - Frobenius norm of C
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
% VERSION:  1.1
%
% See also: eigensystemRealisation, constructHankelMatrices,
%           svdTruncation, extractModalParameters
%% Input validation
n = size(S, 1);
assert(size(U, 2) == n && size(V, 2) == n, ...
    'stateSpaceMatrixConstruction:truncationMismatch', ...
    'U, V column dimensions must match S size — SVD truncation inconsistent.');

assert(isequal(size(shiftedHankelMatrix), [size(U, 1), size(V, 1)]), ...
    'stateSpaceMatrixConstruction:hankelSizeMismatch', ...
    'shiftedHankelMatrix dimensions must match U row count by V row count.');

%% Main code

% Basis vectors. Eq (22-24).
columnBasisVector = zeros(nColumns, 1);
columnBasisVector(1) = 1;

rowBasisVector = zeros(nRows, 1);
rowBasisVector(1) = 1;

% Half power and inverse matrices. Transform to avoid computing over the
% entire S matrix. Only saves time for large S (assumed).
singularValues  = diag(S);
sHalf = diag(sqrt(singularValues));       
sHalfInv = diag(1 ./ sqrt(singularValues)); 

% Eq. (22-24).
A = sHalfInv * U' * shiftedHankelMatrix * V * sHalfInv;
B = sHalf * V' * columnBasisVector;
C = sHalf * U' * rowBasisVector;
D = firstSample;

% Transpose C to enforce canonical row vector convention (1 x n).
% Note: To take out this step, C matrix form could be written as:
%
% C = rowBasisVector' * U * sHalf;
%
% This has been left as is as it matches the paper implementation.
C = C';

%% Optional diagnostics output
if nargout > 4 && isfield(opts, 'returnDebug') && opts.returnDebug

    diagnostics.stateDimension   = n;
    diagnostics.conditionNumberA = cond(A);
    diagnostics.BNorm            = norm(B);
    diagnostics.CNorm            = norm(C);

    varargout{1} = diagnostics;
end

end

