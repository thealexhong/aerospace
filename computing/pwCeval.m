  function Cvals = pwCEval(a,b,c,d,x,zVals)
% Cvals = pwCEval(a,b,c,d,x,zVals)
%
% Evaluates the piecewise cubic polynomial defined by the column (n-1)-vectors a,b,c, and
% d and the column n-vector x. It is assumed that x(1) < ... < x(n).
% zVals  is a column m-vector with each component in [x(1),x(n)].
%
% CVals is a column m-vector with the property that CVals(j) = C(zVals(j)) 
% for j=1:m where on the interval [x(i),x(i+1)]
%
%   C(z)= a(i) + b(i)(z-x(i)) + c(i)(z-x(i))^2 + d(i)(z-x(i))^2(z-x(i+1))
 
m = length(zVals); 
Cvals = zeros(m,1); 
g=1;
for j=1:m
   i = Locate(x,zVals(j),g);
   Cvals(j) = d(i)*(zVals(j)-x(i+1)) + c(i);
   Cvals(j) = Cvals(j)*(zVals(j)-x(i)) + b(i);
   Cvals(j) = Cvals(j)*(zVals(j)-x(i)) + a(i);
   g = i;
end