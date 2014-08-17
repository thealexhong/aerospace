  function numI = QGL(fname,a,b,m)
% numI = QGL(fname,a,b,m,tol)
%
% Integrates a function from a to b 
% fname is a string that names an available function of the form f(x) that 
% is defined on [a,b]. f should return a column vector if x is a column vector.
% a,b are real scalars, m is an integer that satisfies 2 <= m <= 6.
%
% numI is the m-point Gauss-Legendre approximation of the 
% integral of f(x) from a to b. 

[w,x] = GLWeights(m);
fvals = feval(fname,((b-a)/2)*x + ((a+b)/2)*ones(m,1));
numI = ((b-a)/2)*w'*fvals;