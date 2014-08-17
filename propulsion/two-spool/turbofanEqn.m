% 4 turbofan engine equations with 4 unknowns

function F = turbofanEqn(x, T_05, f, T_03f, T_02, P_05, gamma, e_t, ...
                         c_p, n_n, P_a, u_ef, u, thrust_flowrate)
T_05p = x(1);
P_05p = x(2);
beta = x(3);
u_e = x(4);

F(1) = T_05 - (1 + beta)/(1 + f) * ...
       (T_03f - T_02) -  T_05p;
F(2) = P_05 * ...
       (T_05p / T_05)^(gamma / (e_t * (gamma - 1))) - P_05p;
F(3) = (2 * c_p * n_n * T_05p * ...
       (1 - (P_a / P_05p)^((gamma - 1) / gamma))) - u_e^2;
F(4) = ((1 + f) / (1 + beta)) * u_e + beta / (1 + beta) * u_ef - u - ...
        thrust_flowrate;

end