% Alexander Hong (997584706)
% November 20, 2013

% modalfrfanalysis.m
% This code approximates the frequency response using a
...modal analysis approach
function [R_freq_response] = ...
    modalfrfanalysis (K, mode_shapes, natural_freq, excite_freq, mesh, ...
    nModes)
% Variables
% nModes: number of modes range
% q: intermediate value needed to calculate displacement response
% u: displacement response for all joints
% R_freq_response: frequency response at point of interest

% Assemble harmonic excitation vector and get index of joint of interest
[f, r_indx] = assembleforce (mesh, -1, K);

% For each mode, calculate q for all frequency ranges
for k = 1:nModes
    for j = 1:size (excite_freq, 2)
        q{k, j} = ((mode_shapes(:, k))' * f) / ...
             ((natural_freq(k) * (2 * pi))^2 - excite_freq(j)^2 + ...
             sqrt(-1) * excite_freq(j) * 10);
    end
end

% Displacement response is the sum of the product of q and mode shapes
for j = 1:size (excite_freq, 2)
    u{j} = zeros(size(K, 1), 1);
    for k = 1:nModes
        u{j} = u{j} + q{k, j} * mode_shapes(:,k);
        R_freq_response{k}(j, 1) = abs(u{j}(r_indx));
    end
end