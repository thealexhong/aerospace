function [s, u] = pkm_inv(pose) % pose is 4x4 matrix

% PKM parameter
[R_b, R_p, L] = getDim;

% Process inputs
P_C = pose(1:3,4); % translational components

% Rotation matrix
R = pose(1:3, 1:3);

% Spherical joints with respect to platform frame
for i = 1:3
    P_P{i} = [R_p * cos(2 / 3 * pi * (i - 1)),
              R_p * sin(2 / 3 * pi * (i - 1)),
              0];
end

% Home configuration
s_(1:3) = 0;
u_(1:3) = 0;
phi_(1:3) = acos((R_b - R_p)/L);

% Initial guess - home positions
x_0 = [u_, phi_, s_];
%opts = optimset('Display', 'off', 'MaxFunEvals', 10000, 'MaxIter', 10000);
opts = optimoptions('fsolve','Display','off');
[x, ~] = fsolve(@(x) pkm_eqn_inv(x, R, P_P, P_C, L, R_b), x_0, opts);

% Assign to problem parameters
u = x(1:3);
phi = x(4:6);
s = x(7:9);
end

function [R_b, R_p, L] = getDim
% get dimension of PKM
R_b = 150; % mm (Radius of base)
R_p = 34.89; % mm (Radius of platform)
L = 164; % mm (Link length)
end

% Rotate any frame with respect to z-axis
function Rz = rotateRz(rad)
Rz = [cos(rad), -sin(rad), 0;
      sin(rad),  cos(rad), 0;
             0,         0, 1;];
end

% System of non-linear equation to solve for inverse kinematics
function F = pkm_eqn_inv(x, R, P_P, P_C, L, R_b)
F = [
    % SOLVE for:
    % x(1:3) = u_i
    % x(4:6) = phi_i
    % x(7:9) = s_i
    
    R * P_P{1} + P_C - [R_b + x(1) - L * cos(x(4));...
                        x(7);
                        L * sin(x(4))];
    R * P_P{2} + P_C - rotateRz(2 / 3 * pi) * ...
                        [R_b + x(2) - L * cos(x(5));...
                        x(8);
                        L * sin(x(5))];
    R * P_P{3} + P_C - rotateRz(4 / 3 * pi) * ...
                        [R_b + x(3) - L * cos(x(6));...
                        x(9);
                        L * sin(x(6))];
    ];
end