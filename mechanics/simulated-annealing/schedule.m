% schedule.m
% This function returns temperature based on the exponential cooling
... schedule T(t) = Tstart * c^t, where 0 < c < 1. Tstart = 1000.

function [T] = schedule(Tstart, c, t)
% Variables
% T: Temperature based on exponential cooling schedule
% c: cooling schedule parameter
% t: time
T = Tstart * (c^t);
