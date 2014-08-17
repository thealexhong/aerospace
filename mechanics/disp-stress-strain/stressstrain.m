function [sigma, epsilon] = stressstrain(num_elements, coords, connect, d, E)

epsilon = zeros(num_elements, 1); 
for i = 1:num_elements
    % coordinate of first element
    P1 = coords(connect(i,1),:); 
  
    % coordinate of second element
    P2 = coords(connect(i,2),:); 
   
    l = norm(P2-P1);    % length of the element 
    theta = atan2(P2(2)-P1(2), P2(1)-P1(1));
    
    T = zeros(4,4);
    subT = [cos(theta) sin(theta); -sin(theta) cos(theta)];
    T = [subT, zeros(2,2); zeros(2,2), subT];
    % disp(T)
    newd = d([2*connect(i,1)-1, 2*connect(i,1), 2*connect(i,2)-1, 2*connect(i,2)]);
    % disp(newd)
    locald = T * newd; 
    epsilon(i) = (locald(3) - locald(1))/l; 
end

% compute stress
sigma = zeros(num_elements);    
sigma = E*epsilon;