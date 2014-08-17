  function numI = CompQNC(fname,a,b,m,n)
% numI = CompQNC(fname,a,b,m)
%
% Integrates a function of the form f(x) named by the string fname from a to b. 
% f must be defined on [a,b] and it must return a column vector if x is a column vector.
% m is an integer that satisfies 2 <= m <= 11.
% numI is the composite m-point Newton-Cotes approximation of the integral of f 
% from a to b with n equal length subintervals.
 
Delta = (b-a)/n;
h = Delta/(m-1);
x = a+h*(0:(n*(m-1)))';  
w = NCWeights(m);
x = linspace(a,b,n*(m-1)+1)';
f = feval(fname,x);
numI = 0; 
first = 1; 
last = m;
for i=1:n
   %Add in the inner product for the i-th subintegral.
   numI = numI + w'*f(first:last);
   first = last;
   last = last+m-1;
end 
numI = Delta*numI;