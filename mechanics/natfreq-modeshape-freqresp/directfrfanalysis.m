% Alexander Hong (997584706)
% November 20, 2013

% directfrfanalysis.m
% This code calculates the frequency response using the full-order ...
% dynamic stiffness matrix
function [R_freq_response] = ...
    directfrfanalysis (M, K, excite_freq, mesh)
% Variables
% C: Damping model given in the problem
% f_indx: force being applied to the specified DOF
% r_indx: joint R (joint of interest for frequency response)
% u: displacement response for all joints
C = 10 * M;
[f, r_indx] = ...
    assembleforce (mesh, -1, K);
% Use equation (6.52) in notes to solve for full-order solution
for j = 1:size (excite_freq, 2)
    u{j} = (K - (excite_freq(j)^2) * M + ...
        sqrt(-1) * excite_freq(j) * C) \ f;
    R_freq_response(j, 1) = abs(u{j}(r_indx));
end