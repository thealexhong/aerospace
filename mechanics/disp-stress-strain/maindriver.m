%% maindriver.m

%% qn2

% clear all
% clc

numTruss = 10; 

% Parameters
E = 70e9; 
A = 1e-4; 
P = 100;

if numTruss == 3
    num_nodes = 4; 
    num_elements = 3;
    coords = [ -1, 1; 0, 1; 1, 1; 0, 0];
    connect = [1, 4; ...    % element 1
                2, 4; ...    % element 2
                3, 4];       % element 3
    % Force vector
    F = zeros(2*num_nodes,1); 
    ind_knownForce = 4;
    F_known = [P, 0]'; 
    F(2*ind_knownForce-1: 2*ind_knownForce) = F_known; 
    nodes_knownDisp_x = [1, 2, 3];
    nodes_knownDisp_y = [1, 2, 3];
elseif numTruss == 10
    num_nodes = 6; 
    num_elements = 10; 
    coords = [0, 1; 1, 1; 2, 1; 2, 0; 1, 0; 0, 0];
    connect = [1, 2; 2, 3; 1, 5; 2, 6; 2, 5; 2, 4; 3, 5; 3, 4; 4, 5; 5, 6];
    % Force vector
    F = zeros(2*num_nodes,1); 
    ind_knownForce = [4, 5];
    F_known = [0, -P, 0, -P]'; 
    nodes_knownDisp_x = [1, 6];
    nodes_knownDisp_y = [1, 6];
    BC = [0.0, 0.0, 0.0, 0.0];  % known displacements
    F(2*ind_knownForce(1)-1: 2*ind_knownForce(2)) = F_known; 
end
           

% build the global stiffness matrix 

K = assembly(E, A, num_nodes, num_elements, coords, connect);

ind_unknownDisp = [];   % will contain indices for unknown displacements
d = zeros(2*num_nodes,1); 
ind_knownDisp = [];     % will contain indices for known displacements

for i = 1:num_nodes
    if ~ismember(i, nodes_knownDisp_x)
        ind_unknownDisp(end+1) = 2*i-1;
    else
        ind_knownDisp(end+1) = 2*i-1;
    end
    if ~ismember(i, nodes_knownDisp_y)
        ind_unknownDisp(end+1) = 2*i;
    else
        ind_knownDisp(end+1) = 2*i;
    end
end

ind_unknownDisp = ind_unknownDisp(:);
ind_knownDisp = ind_knownDisp(:);

% The K matrix to account for known boundary conditions 
Kdisp = K(ind_unknownDisp, ind_knownDisp);

% reduced system of equation 
redK = K(ind_unknownDisp, ind_unknownDisp); 
redF = F(ind_unknownDisp) - Kdisp*BC';  % subtracted by the force due to displacements
redd = redK\redF;

d(ind_unknownDisp) = redd;
d(ind_knownDisp) = BC';

disp(d)

[sigma, epsilon] = stressstrain(num_elements, coords, connect, d, E)

