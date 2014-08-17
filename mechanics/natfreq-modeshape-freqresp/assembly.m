% Alexander Hong (997584706)
% November 20, 2013

% assembly.m
% This code assembles the global stiffness and mass matrices.

function [K, M] = assembly (EA, EI, RHO, nElements, ...
    coords, connec, connec_dof, nDOF)
% Variables
% EA: axial rigidity
% EI: flexural rigidity
% RHO: mass per unit length
% nElements: number of elements
% coods: nodal coordinates
% connec: connectivity matrix
% connec_dof: relationship between local and global DOF numbering
% nDOF: number of DOF


%% Assemble local stiffness and mass matrices
K = zeros(nDOF, nDOF);
M = zeros(nDOF, nDOF);
nDOF_element = size(connec_dof,2);
K_element = {};
M_element = {};
for i = 1:nElements
    X1 = coords(connec(i,1),1);
    Y1 = coords(connec(i,1),2);
    X2 = coords(connec(i,2),1);
    Y2 = coords(connec(i,2),2);
    [K_element{i}] = ElementStiffMat(EA, EI, X1, Y1, X2, Y2);
    [M_element{i}] = ElementMassMat(RHO, X1, Y1, X2, Y2);
end

%% Assemble global stiffness and mass matrices
for i = 1:nElements
    for j = 1:nDOF_element
        for k = 1:nDOF_element
            % local to global DOF mapping info
            ip = connec_dof(i,j);
            jp = connec_dof(i,k);
            % Add element contribution
            K(ip, jp) = K(ip, jp) +...
                K_element{i}(j,k);
            M(ip, jp) = M(ip, jp) +...
                M_element{i}(j,k);
        end
    end
end
