  function i = Locate(x,z,g)
% i = Locate(x,z,g)
% Locates z in a partition x.
%  
% x is column n-vector with x(1) < x(2) <...<x(n) and
% z is a scalar with x(1) <= z <= x(n).
% g (1<=g<=n-1) is an optional input parameter
%
% i is an integer such that x(i) <= z <= x(i+1). Before the general 
% search for i begins, the value i=g is tried.

if nargin==3 
   % Try the initial guess.
   if (x(g)<=z) & (z<=x(g+1))
      i = g;
      return
   end
end
n = length(x);
if z==x(n)
   i = n-1;
else
   % Binary Search
   Left = 1; 
   Right = n;
   while Right > Left+1
      % x(Left) <= z <= x(Right)
      mid = floor((Left+Right)/2);
      if z < x(mid) 
         Right = mid;
      else 
         Left = mid;
      end
   end
   i = Left;
end