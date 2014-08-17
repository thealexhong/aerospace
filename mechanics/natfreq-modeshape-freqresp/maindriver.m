% Alexander Hong (997584706)
% November 20, 2013

% AER501: Advanced Mechanics of Structures
% Assignment 2
% This MATLAB code computes natural frequencies, mode shapes, and
... frequency responses of two-dimensional frame structures.

% maindriver.m
% This code is the main function that controls the flow of the program

function maindriver ()

% clear
clear;
clc;
close all;
fclose all;
warning('off','all');
warning;

%% Compute natural frequencies and mode shapes
% Physical properties of the frame structure
physical_input();

% Compute natural frequencies for each mesh
Mesh1();
[natural_freq{1}, mode_shapes{1}, K{1}, M{1}] = ...
    freevibration (X, Y, NOD, EA, EI, RHO, dof_active, 1);
Mesh2();
[natural_freq{2}, mode_shapes{2}, K{2}, M{2}] = ...
    freevibration (X, Y, NOD, EA, EI, RHO, dof_active, 2);
Mesh3();
[natural_freq{3}, mode_shapes{3}, K{3}, M{3}] = ...
    freevibration (X, Y, NOD, EA, EI, RHO, dof_active, 3);
Mesh4();
[natural_freq{4}, mode_shapes{4}, K{4}, M{4}] = ...
    freevibration (X, Y, NOD, EA, EI, RHO, dof_active, 4);

% Graph Frequency vs. Eigenmode number for all meshes
x = 1:12;
figure(1);
plot (x, natural_freq{1}(1:12), '--+', ...
    x, natural_freq{2}(1:12), '-o', ...
    x, natural_freq{3}(1:12), ':*', ...
    x, natural_freq{4}(1:12), '-..');
legend ('Mesh 1','Mesh 2','Mesh 3','Mesh 4');
xlabel ('Eigenmode Number');
ylabel ('Frequency (Hz)');
title ('Convergence of the First 12 Eigenvalues');

%% Displacement response of node R
% Define excitation frequency (rad/s)
excite_freq = (0:1:250) * (2 * pi);

% Calculate frequency response using FULL-ORDER DYNAMIC STIFFNESS
... MATRIX and graph it on a logarithmic scale for displacement
% For meshes 1 to 4
for i = 1:4
    [R_freq_response{1, i}] = directfrfanalysis (M{i}, K{i}, excite_freq, i);
    figure(2);
    semilogy ((excite_freq / (2 * pi)), R_freq_response{1, i});
    hold all
end
legend ('Mesh 1','Mesh 2','Mesh 3','Mesh 4');
xlabel ('Frequency (Hz)');
ylabel ('Displacement Response');
title (['Total displacement response at joint R as a function of the'...
    ' excitation frequency']);
hold off

% Approximate the frequency response using a modal analysis approach
... and graph it on a logarithmic scale for displacement
% For mesh 4
mesh = 4;
nModes = 50; % number of modes // converges at mode 17
[R_freq_response{2}] = ...
        modalfrfanalysis (K{mesh}, mode_shapes{mesh}, natural_freq{mesh},...
        excite_freq, mesh, nModes);

colour = {'b', 'm', 'c', 'r', 'g', 'b', 'k', [.3 .7 .7], [.8 .2 .6],...
    [.5 .5 .5]};
c = 0;
% Graph the following modes in k
for mode = [1 3 10 15 16 17 18 30]
    c = c + 1;
    figure(3);
    semilogy ((excite_freq / (2 * pi)), R_freq_response{2}{mode},...
        'color', colour{c});
    hold all
end
% Also plot full-order solution
semilogy ((excite_freq / (2 * pi)), R_freq_response{1, 4}, 'color', [.6 .2 .9]);
legend ('Mode 1','Mode 3','Mode 10','Mode 15','Mode 16',...
    'Mode 17','Mode 18','Mode 30',...
    'Full-Order Solution','Location','EastOutside');
xlabel ('Frequency (Hz)');
ylabel ('Displacement Response');
title (['Total displacement response at joint R as a function of the'...
    ' excitation frequency using Modal Analysis']);
hold off

% Error analysis between full-order and modal analysis
for mode = 1:50
    error(mode) = abs (mean(R_freq_response{1, 4} - R_freq_response{2}{mode}));   
end
x = 1:50;
figure(4);
plot (x, error);
xlabel ('Mode Number');
ylabel ('Error Response');
title ('Average error between the full-order and modal analysis');





