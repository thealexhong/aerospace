  function [x,fx,nEvals,aF,bF] = GlobalNewton(fName,fpName,a,b,tolx,tolf,nEvalsMax)
% [ x,fx,nEvals,aF,bF] = GlobalNewton(fName,fpName,a,b,tolx,tolf,nEvalsMax)
% fName       string that names a function f(x).
% fpName      string that names the derivative function f'(x).
% a,b         A root of f(x) is sought in the interval [a,b]
%               and f(a)*f(b)<=0.
% tolx,tolf   Nonnegative termination criteria.
% nEvalsMax   Maximum number of derivative evaluations.
%
% x          An approximate zero of f.
% fx         The value of f at x.
% nEvals     The number of derivative evaluations required. 
% aF,bF      The final bracketing interval is [aF,bF].
%
% The iteration terminates as soon as x is within tolx of a true zero or
% if |f(x)|<= tolf or after nEvalMax f-evaluations.

fa  = feval(fName,a);
fb  = feval(fName,b);
if fa*fb>0
   disp('Initial interval not bracketing.')
   return
end
x   = a;
fx  = feval(fName,x);
fpx = feval(fpName,x);
disp(sprintf('%20.15f  %20.15f    %20.15f',a,x,b))

nEvals = 1;
while (abs(a-b) > tolx ) & (abs(fx) > tolf) & ((nEvals<nEvalsMax) | (nEvals==1))
   %[a,b] brackets a root and x = a or x = b.
   if StepIsIn(x,fx,fpx,a,b)
      %Take Newton Step
      disp('Newton')
      x   = x-fx/fpx;
   else
      %Take a Bisection Step:
      disp('Bisection')
      x = (a+b)/2;
   end
   fx  = feval(fName,x);
   fpx = feval(fpName,x); 
   nEvals = nEvals+1;  
   if fa*fx<=0
      % There is a root in [a,x]. Bring in right endpoint.
      b  = x; 
      fb = fx;
   else
      % There is a root in [x,b]. Bring in left endpoint.
      a  = x; 
      fa = fx;
   end 
   disp(sprintf('%20.15f  %20.15f    %20.15f',a,x,b))
end
aF = a;
bF = b;