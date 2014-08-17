%
% |>o-----o-----o-----o
%    \    |\    |\    |
%      \  |  \  |  \  |
%        \|    \|    \|
% |>o-----o-----o-----o
%
% Mesh with four elements per member
%

% connectivity matrix
NOD = [ 1    2; 2    3; 3    4; 4    5; 5    6; 6    7; 7    8; 8    9; 9   10; 10   11; 11   12; 12   13; 13   14; 14   15; 15   16; 16    5; 5   17; 17   18; 18   19; 19   20; 20   21; 21   22; 22   23; 23   24; 24   25; 25   26; 26   27; 27    9; 9   28; 28   29; 29   30; 30   20; 20   31; 31   32; 32   33; 33   34; 34   35; 35   36; 36   37; 37   38; 38   39; 39   40; 40   41; 41   24; 24   42; 42   43; 43   44; 44   34]; 

% vector containing x-coordinates of each node
X = [   0.00000 0.25000 0.50000 0.75000 1.00000 1.00000 1.00000 1.00000 1.00000 0.75000 0.50000 0.25000 0.00000 0.25000 0.50000 0.75000 1.25000 1.50000 1.75000 2.00000 2.00000 2.00000 2.00000 2.00000 1.75000 1.50000 1.25000 1.25000 1.50000 1.75000 2.25000 2.50000 2.75000 3.00000 3.00000 3.00000 3.00000 3.00000 2.75000 2.50000 2.25000 2.25000 2.50000 2.75000 ];

 
% vector containing y-coordinates of each node
Y = [ 0.00000 0.00000 0.00000 0.00000 0.00000 0.25000 0.50000 0.75000 1.00000 1.00000 1.00000 1.00000 1.00000 0.75000 0.50000 0.25000 0.00000 0.00000 0.00000 0.00000 0.25000 0.50000 0.75000 1.00000 1.00000 1.00000 1.00000 0.75000 0.50000 0.25000 0.00000 0.00000 0.00000 0.00000 0.25000 0.50000 0.75000 1.00000 1.00000 1.00000 1.00000 0.75000 0.50000 0.25000]; 

% Indices of constrained nodes
constr_node(1) = 1;
constr_node(2) = 13;

% Generate indices of active and inactive DOFs
itmp = constr_node(2) - 1;
dof_restr = [1, 2, 3, 3*itmp + 1, 3*itmp + 2, 3*itmp + 3]; % indices of restrained DOF 
dof_active = [4:3*itmp, 3*itmp+4:3*length(X)]; % indices of unrestrained DOF

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



