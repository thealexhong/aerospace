% Alexander Hong (997584706)
% November 20, 2013

% assembleforce.m
% This code assembles the harmonic excitation vector and returns the DOF
... numbering of joint R
function [f, r_indx] = ...
    assembleforce (mesh, force, K)
% Variables
% f_indx: force being applied to the specified DOF
% r_indx: joint R (joint of interest for frequency response)
f = zeros(size(K, 1), 1);
% Hard code f_indx and r_indx as user input can vary
switch mesh
    case 1
        f_indx = (2);
        r_indx = ((7 * 3 - 1) - 6);
    case 2
        f_indx = (5);
        r_indx = ((16 * 3 - 1) - 6);
    case 3
        f_indx = (8);
        r_indx = ((25 * 3 - 1) - 6);
    case 4
        f_indx = (11);
        r_indx = ((34 * 3 - 1) - 6);
end
f(f_indx) = force; % force being applied