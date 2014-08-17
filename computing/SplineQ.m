  function numI = SplineQ(x,y)
% numI = SplineQ(x,y)
%
% Integrates the spline interpolant of the data specified by the
% column n-vectors x and y. It is a assumed that x(1) < ... < x(n)
% and that the spline is produced by SPLINE. The integral is from
% x(1) to x(n).

S = spline(x,y);
[x,rho,L,k] = unmkpp(S);
sum = 0;
for i=1:L
   % Add in the integral from x(i) to x(i+1).
   h = x(i+1)-x(i);
   subI = h*(((rho(i,1)*h/4 + rho(i,2)/3)*h + rho(i,3)/2)*h + rho(i,4));
   sum = sum + subI;
end
numI = sum;