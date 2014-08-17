function [pose, XYZpose, jointPos] = pkm_fwd(s, u)

% PKM parameter
[R_b, R_p, L] = getDim;
l = sqrt(2 * R_p^2 * (1 - cos(2 / 3 * pi)));

% Home configuration
s_(1:3) = 0;
u_ (1:3) = 0;
phi_(1:3) = acos((R_b - R_p)/L);

% Initial guess - close to home position as workspace is very small
phi_0 = phi_;

% opts = optimset('Display', 'off', 'MaxFunEvals', 10000, 'MaxIter', 10000);
opts = optimoptions('fsolve','Display','off');
[phi, ~] = fsolve(@(phi) pkm_eqn(phi, R_b, u, s, L, l), phi_0, opts);

% Calculate positions of spherical joints
for i = 1:3
    P{i} = rotateRz(2 / 3 * pi * (i - 1)) * ...
              [R_b + u(i) - L * cos(phi(i));
               s(i);
               L * sin(phi(i))];
end

jointPos = P;
% Pose of end effector frame
pose = getPose(P, l);
XYZpose = getXYZPose(P, l);
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

% Forward Kinematic Equations
function F = pkm_eqn(phi, R_b, u, s, L, l)
F = [
     (R_b + u(1) - L * cos(phi(1)) - ...
        ((R_b + u(2) - L * cos(phi(2))) * cos(2 / 3 * pi) - ...
            s(2) * sin(2 / 3 * pi)))^2 + ...
     (s(1) - ((R_b + u(2) - L * cos(phi(2))) * sin(2 / 3 * pi) + ...
        s(2) * cos(2 / 3 * pi)))^2 + ...
     (L * sin(phi(1)) - L * sin(phi(2)))^2 - l^2;
     (R_b + u(1) - L * cos(phi(1)) - ...
        ((R_b + u(3) - L * cos(phi(3))) * cos(4 / 3 * pi) - ...
            s(3) * sin(4 / 3 * pi)))^2 + ...
     (s(1) - ((R_b + u(3) - L * cos(phi(3))) * sin(4 / 3 * pi) + ...
        s(3) * cos(4 / 3 * pi)))^2 + ...
     (L * sin(phi(1)) - L * sin(phi(3)))^2 - l^2;
     (((R_b + u(2) - L * cos(phi(2))) * cos(2 / 3 * pi) - ...
            s(2) * sin(2 / 3 * pi)) - ...
        ((R_b + u(3) - L * cos(phi(3))) * cos(4 / 3 * pi) - ...
            s(3) * sin(4 / 3 * pi)))^2 + ...
     (((R_b + u(2) - L * cos(phi(2))) * sin(2 / 3 * pi) + ...
        s(2) * cos(2 / 3 * pi)) - ...
        ((R_b + u(3) - L * cos(phi(3))) * sin(4 / 3 * pi) + ...
        s(3) * cos(4 / 3 * pi)))^2 + ...
     (L * sin(phi(2)) - L * sin(phi(3)))^2 - l^2;
    ];
end

% Get pose from 3 points (P_W)
function [pose, P_C, R] = getPose (P_W, l)
P_C = [1 / 3 * (P_W{1} + P_W{2} + P_W{3})];

% Solve for rotation matrix
A = (P_W{1} - P_C) / (l / sqrt(3));
B = (P_W{2} - P_W{3}) / l;
R = [A, B, cross(A, B)];

% a 4x4 matrix
pose = [R, P_C; 0, 0, 0, 1];
end

% XYZ Euler Angles of pose
function [XYZpose] = getXYZPose (P_W, l)
[~, P_C, R] = getPose(P_W, l);

% [Rx, Ry, Rz] = [alpha, beta, gamma];
alpha = atan2(-R(3,1), sqrt(R(1,1)^2 + R(2,1)^2));
beta = atan2(R(2,1) / cos(alpha), R(1,1) / cos(alpha));
gamma = atan2(R(3,2) / cos(alpha), R(3,3) / cos(alpha));

XYZpose = [P_C', alpha, beta, gamma];
end