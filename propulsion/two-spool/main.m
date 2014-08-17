% This program minimizes a cost function for a two-spool turbofan
... engine.
    
clear;
clc;
fclose all;
close all;
format short;

%% Important Engine parameters
% Thrust per unit flow rate
thrust_flowrate = 125; % Ns/kg  (Fixed specific thrust)
% Cruise condition
M = 0.85; % Mach
P_a = 24000; % ambient pressure (Pa)
T_a = 220; % ambient temperature (K)

% Two spool Engine with 
% i) high pressure turbine in core compressor
% ii) low pressure turbine in fan

% Efficiencies
n_d = 0.98; % inlet
n_b = 0.99; % fuel combustion
n_n = 0.95; % nozzle (fan & core)
e_c = 0.91; % compressor polytropic
e_f = 0.93; % fan polytropic
e_t = 0.93; % turbine polytropic (high & low pressure)
pi_b = 0.95; % combustor pressure loss
Q_r = 48E6; % fuel heating value (J/kg)
T_04 = 1380; % turbine inlet temperature (K)
% no mechanical losses n_m = 1;
% nozzle expands flow to atmospheric pressure
gamma = 1.4;
c_p = 1000; %(J/kg)

% Design space:
% i) pi_c     %core compression ratio
% ii) pi_f    %fan compression ratio
% iii) beta   %bypass ratio

% Cost Function:
% Minimize C = W + F;
% W = beta_0 / beta;
% F = TSFC/TSFC_0;
beta_0 = 4;
TSFC_0 = 2.5E-5; % (kg/Ns)
% (fixed specific thrust + fixed engine diameter
... = weight increase w/ decrease beta)


%% Calculations
% Fan inlet conditions
T_02 = T_a * (1 + (gamma - 1)/2 * M^2);
P_02 = P_a * (1 + n_d * (T_02/T_a - 1))^(gamma / (gamma - 1));

% Vary core compression ratio, fan compression ratio
pi_f = 1.4:0.1:2.0;
pi_c = 3:1:60;

R = (gamma - 1) * c_p / gamma;
u = M * sqrt (gamma * R * T_a);

for i = 1:size(pi_f, 2)
    % Fan outlet conditions
    P_03f(i) = P_02 * pi_f(i);
    T_03f(i) = T_02 * (pi_f(i))^((gamma - 1)/(e_f * gamma));
        
    % Exit at fan nozzle
    u_ef(i) = sqrt (2 * c_p * n_n * T_03f(i) * ...
        (1 - (P_a/P_03f(i))^((gamma - 1) / gamma)));
    
    for j = 1:size(pi_c, 2)
        % Compressor outlet condition
        P_03(i, j) = P_03f(i) * pi_c(j);
        T_03(i, j) = T_03f(i) * (pi_c(j))^((gamma - 1)/(e_c * gamma));
        
        % f value
        f(i, j) = (T_04 / T_03(i, j) - 1) / ...
                ((n_b * Q_r) / (c_p * T_03(i,j)) - T_04 / T_03(i,j));
        P_04(i,j) = P_03(i, j) * pi_b;
        
        % First Turbine
        T_05(i, j) = T_04 - 1/(1 + f(i,j)) * (T_03(i,j) - T_03f(i));
        P_05(i, j) = P_04(i,j) * ...
                (T_05(i,j) / T_04)^(gamma / (e_t * (gamma - 1)));
        
        % Initial guess at variables:
        % T_05p = x(1)
        % P_05p = x(2)
        % beta = x(3)
        % u_e = x(4)
        x_0 = [600, 40000, 8, 400];
        % Option to display output
        opts = optimset('Display','off','MaxFunEvals',10000,'MaxIter',10000);
        [x, fval] = fsolve(@(x) turbofanEqn(x, T_05(i,j), f(i,j),...
            T_03f(i), T_02, P_05(i,j),gamma, e_t, c_p, n_n, P_a, ...
            u_ef(i), u, thrust_flowrate), x_0, opts);
        T_05p(i,j) = x(1);
        P_05p(i,j) = x(2);
        beta(i,j) = x(3);
        u_e(i,j) = x(4);
        
        % Calculate TSFC, and cost function
        TSFC(i,j) = f(i,j) / (thrust_flowrate * (1 + beta(i,j)));
        W(i,j) = beta_0 / beta(i,j);
        F(i,j) = TSFC(i,j) / TSFC_0;
        C(i,j) = W(i,j) + F(i,j);
        
    end

end

%% Results

% Contour Plots
figure(1);
contourf(pi_f, pi_c, beta', 40);
title('Bypass Ratio');
xlabel('\pi_f');
ylabel('\pi_c');
colorbar;

figure(2);
contourf(pi_f, pi_c, TSFC', 40);
title('TSFC (kg/Ns)');
xlabel('\pi_f');
ylabel('\pi_c');
colorbar;

figure(3);
contourf(pi_f, pi_c, C', 40);
title('Cost Function');
xlabel('\pi_f');
ylabel('\pi_c');
colorbar;

figure(4);
contourf(pi_f, pi_c, f', 40);
title('f');
xlabel('\pi_f');
ylabel('\pi_c');
colorbar;

figure(5);
contourf(pi_f, pi_c, T_03', 40);
title('T_0_3');
xlabel('\pi_f');
ylabel('\pi_c');
colorbar;

% Final Engine Design
[C_min, ind] = min(C(:));
[i,j] = ind2sub(size(C),ind); % indices for optimal values
fprintf('FINAL TURBOFAN ENGINE DESIGN\n====================\n');
fprintf('Design Parameters:\n');
fprintf('pi_c = %.1f\n', pi_c(j));
fprintf('pi_f = %.2f\n', pi_f(i));
fprintf('beta = %.2f\n', beta(i,j));
fprintf('\nResultant Values:\n');
fprintf('TSFC = %d kg/(Ns)\n', TSFC(i,j));
fprintf('W = %.3f\n', W(i,j));
fprintf('F = %.3f\n', F(i,j));
fprintf('C = %.2f\n', C(i,j));
fprintf('\nCalculated Parameters:\n');
fprintf('P_02 = %.2e Pa\n', P_02);
fprintf('T_02 = %3.f K\n', T_02);
fprintf('P_03f = %.2e Pa\n', P_03f(i));
fprintf('T_03f = %3.f K\n', T_03f(i));
fprintf('P_03 = %.2e Pa\n', P_03(i,j));
fprintf('T_03 = %3.f K\n', T_03(i,j));
fprintf('P_04 = %.2e Pa\n', P_04(i,j));
fprintf('T_04 = %3.f K\n', T_04);
fprintf('P_05 = %.2e Pa\n', P_05(i,j));
fprintf('T_05 = %3.f K\n', T_05(i,j));
fprintf('P_05p = %.2e Pa\n', P_05p(i,j));
fprintf('T_05p = %3.f K\n', T_05p(i,j));
fprintf('Specific Thrust = %3.f Ns/kg\n', thrust_flowrate);
fprintf('R = %3.f J/(kgK)\n', R);
fprintf('f = %.4f\n', f(i,j));
fprintf('u = %3.f m/s\n', u);
fprintf('u_ef = %3.f m/s\n', u_ef(i));
fprintf('u_e = %3.f m/s\n', u_e(i,j));

