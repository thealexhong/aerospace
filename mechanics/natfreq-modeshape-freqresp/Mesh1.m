%
%  2     4     6     8
% |>o-----o-----o-----o
%    \    |\    |\    |
%      \  |  \  |  \  |
%        \|    \|    \|
% |>o-----o-----o-----o
%   1     3     5     7
%
%
% Mesh with one element per member
%

% connectivity matrix
NOD = [1,3; 3,4; 2,4; 2,3; 3,5; 5,6; 4,6; 4,5; 5,7; 7,8; 6,8; 6,7];

% vector containing x-coordinates of each node
X = [0,0,1,1,2,2,3,3]; 

% vector containing y-coordinates of each node
Y = [0,1,0,1,0,1,0,1];

% Indices of constrained nodes
constr_node(1) = 1;
constr_node(2) = 2;

% Generate indices of active and inactive DOFs
dof_restr = 1:6;
dof_active = 7:24;

%
% Let K and M denote the global stiffness and mass matrices after assembly.
% The rows and columns corresponding to the zero displacement BCs can be
% deleted as follows:
%
% K = K(dof_active, dof_active)
% M = M(dof_active, dof_active)
%
% Using the mesh information and Fig 1 in Assignment2.pdf you can identify 
% the DOF where the external force is applied and the DOF whose displacement 
% response is of interest. 
%


