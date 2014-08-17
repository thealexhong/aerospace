  function pVal = HornerN(c,x,z)
% pVal = HornerN(c,x,z)
% Evaluates the Newton interpolant on z where
% c and x are n-vectors and z is an m-vector.
%
% pVal is a vector the same size as z with the property that if
%
%         p(x) = c(1) +  c(2)(x-x(1))+ ... + c(n)(x-x(1))...(x-x(n-1))
% then 
%         pval(i) = p(z(i)) , i=1:m.

n = length(c); 
pVal = c(n)*ones(size(z));
for k=n-1:-1:1
   pVal = (z-x(k)).*pVal + c(k);
end