function K_element = elementmatrix(A,E,P1,P2)

% This Matlab function generates the 4 x 4 stiffness matrix for a planar
% truss element in global coordinates. The syntax is: 
%           K = elementmatrix(A,E,P1,P2)
% where: A is the cross-sectional area
%        E is the Young's modulus
%        P1 and P2 are vectors containing the (x,y) coordinates of the
%        endpoints.

% compute the length of the element
l = norm(P2-P1);    
theta = atan2(P2(2)-P1(2), P2(1)-P1(1));
% transformation from local to global coordinate
transform = [cos(theta) sin(theta) -cos(theta) -sin(theta)]';
K_element = A*E/l*transform*transform';