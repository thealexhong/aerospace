  function ok = StepIsIn(x,fx,fpx,a,b)
% ok = StepIsIn(x,fx,fpx,a,b)
% Yields 1 if the next Newton iterate is in [a,b] and 0 otherwise.  
% x is the current iterate, fx is the value of f  at x, and fpx is 
% the value of f' at x.

if fpx > 0
   ok = ((a-x)*fpx <= -fx) & (-fx <= (b-x)*fpx);
elseif fpx < 0
   ok = ((a-x)*fpx >= -fx) & (-fx >= (b-x)*fpx);
else
   ok = 0;
end