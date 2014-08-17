% Alexander Hong (997584706)
% November 20, 2013

% precalculations.m
% This code computes essential values for assembly of matrices.

function [coords, connec, nElements, nDOF,...
    connec_dof_matrix] = precalculation (X, Y, NOD)
% Variables
% coords: nodal coordinates for each node
% connec: connectivity matrix
% nNodes: number of nodes
% nElements: number of elements
% nDOF: total number of DOF
% connec_dof_matrix: relationship between local and global DOF numbering

% Nodal coordinates for each node
% (x, y) for each node
% Note: each row is an NODE
coords = [X', Y'];

% Connectivity matrices
% (node 1, node 2) for each ELEMENT
% Note: each row is an ELEMENT
connec = NOD;

% Number of nodes and elements in structure
nNodes = size(coords, 1);
nElements = size(connec, 1);

% DOF labelling for each node
% i.e. for dof = 3, node 1 is labelled [1,2,3], node 2 is labelled
... [4,5,6], etc
dof = 3;
nDOF = dof * nNodes;
nodes_dof = zeros (nNodes, dof);
dof_ind = 1;
for i = 1:size(nodes_dof, 1)
    for j = 1:size (nodes_dof, 2)
        nodes_dof(i, j) = dof_ind;
        dof_ind = dof_ind + 1;
    end
end

% Relationship between local and global degree of freedom (DOF) numbering
connec_dof = cell(nElements, 2);
for i = 1:nElements
    for j = 1:2
        node = nodes_dof(connec(i,j),:);
        connec_dof{i, j} = node;
    end
end
connec_dof_matrix = cell2mat(connec_dof);