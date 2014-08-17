function K = assembly(E, A, num_nodes, num_elements, coords, connect)

K = zeros(2*num_nodes, 2*num_nodes);

K_element = {};
for i = 1:num_elements
   % coordinate of first element
   P1 = coords(connect(i,1),:); 
  
   % coordinate of second element
   P2 = coords(connect(i,2),:); 
   K_element{i} = elementmatrix(A,E,P1,P2);
   
   in1 = 2*connect(i,1)-1; 
   in2 = 2*connect(i,2)-1; 
   
   K(in1:in1+1, in1:in1+1) = K(in1:in1+1, in1:in1+1) + K_element{i}(1:2, 1:2); 
   K(in1:in1+1, in2:in2+1) = K(in1:in1+1, in2:in2+1) + K_element{i}(1:2, 3:4);
   K(in2:in2+1, in1:in1+1) = K(in2:in2+1, in1:in1+1) + K_element{i}(3:4, 1:2); 
   K(in2:in2+1, in2:in2+1) = K(in2:in2+1, in2:in2+1) + K_element{i}(3:4, 3:4); 
end

