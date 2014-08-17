% Alexander Hong (997584706)
% November 20, 2013

% freevibration.m
% This code computes natural frequencies and mode shapes for a given mesh.
function [natural_freq, mode_shapes, K, M] = freevibration (X, Y, NOD,...
    EA, EI, RHO, dof_active, mesh)
% Variables
% K: global stiffness matrix
% M: global mass matrix
% mode_shape: mode shape of the structure when vibrating at the
... corresponding natural frequency
% natural_freq: natural frequencies computed from K and M
%% Precalculations
[coords, connec, nElements, nDOF, connec_dof] ...
    = precalculations (X, Y, NOD);

%% Assembly of global stiffness and mass matrices
[K, M] = assembly (EA, EI, RHO, nElements, ...
    coords, connec, connec_dof, nDOF);
K = K(dof_active, dof_active);
M = M(dof_active, dof_active);

%% Computing natural frequences and mode shapes
% mode_shapes is a matrix whose columns corresponds to eigenvectors
% D is diagonal of smallest magnitude eigenvalues
if (mesh == 4)
     % For modal analysis, we need a bit more eigenvectors
    [mode_shapes, D] = eigs(K, M, 126, 'SM'); 
else
    % only need first 12 eigenvalues
    [mode_shapes, D] = eigs(K, M, 12, 'SM');
end
eigenvalues = eig(D);
natural_freq = sqrt(abs(eigenvalues)) / (2 * pi);