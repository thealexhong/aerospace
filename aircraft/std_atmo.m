% Generates the standard atmosphere covering the altitude from sea-level up to 
% 100km

% h geopotential height (km)
% h_g geometric height (km)
% t temperature (K)
% p pressure (N/m^2)
% d density (kg/m^3)
clc;
clear;

% constants
R_earth = 6400; %radius
g = 9.81; %gravity
R = 287; %universal gas

% slope in standard atm model
a1 = -0.0065;
a2 = 0.003;
a3 = -0.0045;
a4 = 0.004;

% generating variables
h=0:100;
% column matrices
h_g = 0:100;
t = 0:100;
p = 0:100;
d = 0:100;
% generate
for i = h
    n = i + 1; %index shift
    a = h(n) * 1000; %convert to metres
    
    h_g(n) = R_earth * h(n) / (R_earth - h(n));
    
    
    if i > 90
        t(n) = 165.66 + a4 * (a-90000);
        p(n) = p(91) * (t(n)/165.66)^(-g/(a4*R));
        d(n) = d(91) * (t(n)/165.66)^(-g/(a4*R)-1);
    elseif i > 79
        t(n) = 165.66;
        p(n) = p(80)*exp(-g/(R*t(80))*(a-79000));
        d(n) = d(80)*exp(-g/(R*t(80))*(a-79000));
    elseif i > 53
        t(n) = 282.66 + a3 * (a-53000);
        p(n) = p (54) * (t(n)/282.66)^(-g/(a3*R));
        d(n) = d (54) * (t(n)/282.66)^(-g/(a3*R)-1);
    elseif i > 47
        t(n) = 282.66;
        p(n) = p(48)*exp(-g/(R*t(48))*(a-47000));
        d(n) = d(48)*exp(-g/(R*t(48))*(a-47000));
    elseif i > 25
        t(n) = 216.66 + a2 * (a-25000);
        p(n) = p(26) * (t(n)/216.66)^(-g/(a2*R));
        d(n) = d(26) * (t(n)/216.66)^(-g/(a2*R)-1);
    elseif i > 11
        t(n) = 216.66;
        p(n) = p(12) * exp(-g/(R * t(12)) * (a-11000));
        d(n) = d(12)*exp(-g/(R * t(12)) * (a-11000));
    else
        t(n) = 288.16 + a1 * a;
        p(n) = 101325 * (t(n)/288.16)^(-g/(a1*R));
        d(n) = 1.225 * (t(n)/288.16)^(-g/(a1*R)-1);
    end
    
end

% format & display
fprintf (1, 'h (km)\th_G (km) \t Temperature T (K) \tPressure P (N/m^2) \t Density rho (kg/m^3)\n\n')
for i = h
    n = i + 1;
    fprintf (1, '%d\t\t%.2f\t\t\t%.2f\t\t\t%.5E\t\t\t%.3E\n', h(n), h_g(n),t(n),p(n),d(n));
end

