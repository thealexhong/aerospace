% Plots data of a twin-jet airplane with given initial data

clc;
clear;
close all;

% Data
S = 80; % wing planform [m^2]
b = 23; % wing span [m]
W_0 = 320000; % Gross Weight of full tank of fuel[N]
W_f = 130000; % fuel weight[N]
W_1 = W_0 - W_f; % Weight of plane + empty tank of fuel[N]
c_t = 0.69 / 3600; % specific fuel consumption [1/sec]
h = 10000; % altitude[m]
T_A = 2 * 60000; % Maximum thurst for both engines @ sea level[N]
C_L_max = 2.39;
C_D_0 = 0.015; % parasite drag
K = 0.08;
rho_sea = 1.225; % sea level density[kg/m^3]
temp_sea = 288.16; % sea level temperature[K]
a1 = -0.0065; % slope constant for altitude vs temperature
g = 9.81; % acceleration due to gravity [kg m/s^2]

% Drag Polar C_L  vs C_D
C_L = -C_L_max : 0.0001 : C_L_max; % C_L values (function of AoA)
C_D = C_D_0 + K * C_L .^ 2; % Drag Polar
figure;
plot (C_L, C_D);
xlabel ('C_L');
ylabel ('C_D');
title ('Drag Polar C_D vs. C_L');


% P_R and P_A plot @ [sealevel h]
temp_h = temp_sea + a1 * h; % temperature at h
rho_h = rho_sea * (temp_h / temp_sea) ^ 4.25; % density at h
T_A_h = T_A * (rho_h / rho_sea); % available thurst at h


% Sea Level
% Finding V_P_R_min
C_L_P_R_min = sqrt (3*C_D_0/ K);
V_P_R_min = sqrt (2*(W_0/S)/(rho_sea * C_L_P_R_min));

% Finding max speed given available thrust
q_max = roots ([C_D_0 -T_A/S K*(W_0/S)^2]);
v_max = sqrt(2*max(q_max)/rho_sea);

% Plotting P_R
Vs = V_P_R_min:1:v_max+100; % Velocity Range
P_R = 1/2 * rho_sea * S * C_D_0 * Vs.^3 + ...
    1./ Vs  * ((2 * K * S * (W_0 / S) ^ 2) / rho_sea); % Power Req
P_A = T_A * Vs; % Power available
figure;
plot (Vs, P_R, Vs, P_A, 'b--');
hold on;

% At h = 10000 m
% Finding V_P_R_min at h
V_P_R_min_h = sqrt (2*(W_0/S)/(rho_h * C_L_P_R_min));


% Finding max speed given available thrust at h
q_maxh = roots ([C_D_0 -T_A_h/S K*(W_0/S)^2]);
v_maxh = sqrt(2*max(q_maxh)/rho_h);

% Plotting P_R at h
Vh = V_P_R_min_h:1:v_maxh+100; % Velocity Range
P_R_h = 1/2 * rho_h * S * C_D_0 * Vh.^3 + ...
    1./Vh * ((2 * K * S * (W_0 / S)^2) / rho_h); % Power Req
P_A_h = T_A_h * Vh; % power availabe at h
plot (Vh, P_R_h, 'g', Vh, P_A_h, 'g--'); % green line for altitude
xlabel ('V (m/s)');
ylabel ('Power (W)');
title ('Comparison of altitude effects on Required Power and Power Available');
legend ('P_R @ sea level', 'P_A @ sea level', 'P_R @ h = 10000m', ...
    'P_A @ h = 10000m', 'Location', 'BestOutside');


% Finding V_max at altitude 10000m graphically
df = P_R_h - P_A_h; % finding smallest difference between the two function
iIntersection = find (df == min(abs(df))); % finding index of smallest error - approximately their intersection
disp(sprintf ...
    ('V_max at altitude 10000m found graphically is %.2fm/s.\n',Vh(iIntersection)));


% Hodograph @ sea level
V = eps:1:358; % velocity range
q_sea = 1/2 * rho_sea * V.^ 2;
T_R = q_sea * S .* (C_D_0 + K.*(W_0./(q_sea*S))); % thrust required
climbAngle = asin ((T_A - T_R)/W_0); % climb angle
figure;
dh = V.*sin(climbAngle); % y variable
dx = V.*cos(climbAngle); % x variable
plot (dx, dh);
xlabel ('dx/dt (m/s)');
ylabel ('dh/dt (m/s)');
title ('Hodograph at sea level');
iRCmax = find (dh == max (dh)); % find the maximum value's index
disp(sprintf ...
    ('RC_max and V_RC_max at sea level found graphically are %.2fm/s and %.fm/s\nrespectively.\n'...
    ,dh(iRCmax), V(iRCmax)));


% Estimating the absolute and service ceiling
% At sea level
V_RCmax_sea = sqrt (T_A / (S * 3 * rho_sea * C_D_0) * ...
    (1 + sqrt(1 + 12 * C_D_0 * K / (T_A / W_0)^2)));
RC_max_sea = V_RCmax_sea * (T_A / W_0 - 1/2 * rho_sea * V_RCmax_sea^2 * C_D_0 / (W_0/S)...
    - 2 * K * (W_0 / S) / (rho_sea * V_RCmax_sea^2));

% At altitude alt
alt = 0:1:h+4600; % altitude range
temp_alt = temp_sea + a1 * alt; % temperature at alt
rho_alt = rho_sea * (temp_alt / temp_sea) .^ 4.25; % density at alt
T_A_alt = T_A * (rho_alt / rho_sea); % available thurst at alt
V_RCmax_alt = sqrt (T_A_alt ./ (S * 3 * rho_alt * C_D_0) .* ...
    (1 + sqrt(1 + 12 * C_D_0 * K ./ (T_A_alt / W_0).^2))); % velocity of RC_max_alt
% RC_max at alt
RC_max_alt = V_RCmax_alt .* (T_A_alt / W_0 - 1/2 * rho_alt .* V_RCmax_alt.^2 * C_D_0 / (W_0/S)...
    - 2 * K * (W_0 / S) ./ (rho_alt .* V_RCmax_alt.^2));

% absolute ceiling
iAbsCeil = find (RC_max_alt >= 0); % find index of absolute ceiling
% Service ceiling
iServCeil = find (RC_max_alt >= 0.004 & RC_max_alt <= 0.006); % find index close to 0.005 m/s
disp(sprintf ...
    ('Absolute and service ceiling found graphically are %.fm and %.fm\nrespectively.\n\n'...
    ,alt(iAbsCeil(length (iAbsCeil))), alt(iServCeil)));

figure;
plot (RC_max_alt, alt);
title ('Altitude vs RC_m_a_x');
xlabel ('R/C_m_a_x (m/s)');
ylabel ('Altitude (m)');


% Plot of n vs V @ sea level and n_max vs V under C_L_max
VE = eps:10:400; % velocity range for part E (different velocity range for better view of plots)
q_seaE = 1/2 * rho_sea * VE .^ 2; % dynamic pressure
n = sqrt ( q_seaE/(K*(W_0/S)) .* (T_A/W_0 - q_seaE*C_D_0/(W_0/S))); % calculate load factor

VE2 = eps:10:150; % velocity range for n_max
q_seaE2 = 1/2 * rho_sea * VE2 .^ 2; % dynamic pressure for n_max
n_max = q_seaE2 * S * C_L_max / W_0; % max load factor

% Plot
figure;
plot (VE, n, VE2, n_max, '-.');
title ('Load factor n and n_m_a_x vs V');
xlabel ('V (m/s)');
ylabel ('Load factor');
legend ('Available load factor', 'Max load factor');




% V_max at h
disp(sprintf ...
    ('V_max at altitude 10000m found analytically is %.2fm/s.\n\n'...
    ,v_maxh));

% V_RC max and RC_max @ sea level
disp(sprintf ...
    ('RC_max and V_RC_max at sea level found analytically are %.2fm/s and %.2fm/s\nrespectively.\n\n'...
    ,RC_max_sea, V_RCmax_sea));


% maximum range at h
R_max = 2 * sqrt (2 / (rho_h * S)) * 1/c_t * ... % calculate maximum range
    3/4 * (1 / (3 * K * C_D_0^3)) ^ (1/4) * (sqrt(W_0) - sqrt(W_1));
disp(sprintf ...
    ('Range at altitude 10000m found analytically is %.fm.\n\n', R_max));

% minimum level turning radius and its load factor @ sea level
R_c_min = 4 * K * (W_0 / S) / (rho_sea * g * ...
    sqrt((T_A / W_0)^2 - 4 * K * C_D_0)); % calculate turning radius
n_R_c_min = sqrt (2- 4 * K * C_D_0 / (T_A / W_0)^2); % corresponding load factor
disp(sprintf ...
    ('Minimum level turning radius and its load factor at sea level found analytically\nare %.2fm and %.2f respectively.\n\n'...
    ,R_c_min, n_R_c_min));













