function [B, A] = buildFilterCoefficients(fullPoles, fullResidues, firstSample)
% Input check the number of residues and poles.


%%
SOS_SECTIONS_WIDTH = 3;

POLE_RESIDUE_SCALING_FACTOR = 2;

POLE_RESIDUE_STRIDE = 2;

%%
fullPoles = fullPoles(:);
fullResidues = fullResidues(:);

% Select only one half of the poles (complex conjugate pairs).
poles = fullPoles(1 : POLE_RESIDUE_STRIDE : end);
conjPoles = fullPoles(2 : POLE_RESIDUE_STRIDE : end);
residues = fullResidues(1 : POLE_RESIDUE_STRIDE : end);

% Check the number of filters.
nFilters = numel(poles);

% b - filter numerators.
% a - filter denominators.
% d - direct feedthrough term (gain).
%
% Note: The gain terms CANNOT be distributed to the cascaded filter
% sections. If there are a greater number of modes specified in the ERA
% decomposition, this will over-attenuate the gain term and cause an offset
% in the first sample, thus, the direct feedthrough term is folded into the
% first filter section.
B = zeros(nFilters, SOS_SECTIONS_WIDTH);
A = ones(nFilters, SOS_SECTIONS_WIDTH);
D = firstSample;

% In the first filter, we fold in the feedthrough term. Scaled gain to
% every section does not work for overspecified mode orders.
B(1, 1) = D;
B(1, 2) = POLE_RESIDUE_SCALING_FACTOR * (real(residues(1)) - D * ...
    real(poles(1)));
B(1, 3) = D * abs(poles(1)) ^ 2 - ...
    POLE_RESIDUE_SCALING_FACTOR * real(residues(1) * conjPoles(1));

% Define every other modes' filter weights.
B(2 : nFilters, 2) = POLE_RESIDUE_SCALING_FACTOR * real(residues(2 : end));
B(2 : nFilters, 3) = -POLE_RESIDUE_SCALING_FACTOR * real(residues(2 : end) .* ...
    conjPoles(2 : end));
 
% A is unchanged across all filters.
A(:, 2) = -POLE_RESIDUE_SCALING_FACTOR * real(poles);
A(:, 3) = abs(poles) .^ 2;

end