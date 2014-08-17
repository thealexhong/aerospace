  function [a,b,c,d] = CubicSpline(x,y,derivative,muL,muR)
% [a,b,c,d] = CubicSpline(x,y,derivative,muL,muR)
% Cubic spline interpolation with prescribed end conditions.
% 
% x,y are column n-vectors. It is assumed that n >= 4 and x(1) < ... x(n).
% derivative is an integer (1 or 2) that specifies the order of the endpoint derivatives.
% muL and muR are the endpoint values of this derivative.
%
% a,b,c, and d are column (n-1)-vectors that define the spline S(z). On [x(i),x(i+1)], 
%  
%          S(z) =  a(i) + b(i)(z-x(i)) + c(i)(z-x(i))^2 + d(i)(z-x(i))^2(z-x(i+1).
%
% Usage:
%   [a,b,c,d] = CubicSpline(x,y,1,muL,muR)   S'(x(1))  = muL, S'(x(n))  = muR
%   [a,b,c,d] = CubicSpline(x,y,2,muL,muR)   S''(x(1)) = muL, S''(x(n)) = muR
%   [a,b,c,d] = CubicSpline(x,y)             S'''(z) continuous at x(2) and x(n-1)
       

% First, set up all but the first and last equations that
% define the vector of interior knot slopes s(2:n-1).

n = length(x); 
Dx = diff(x);
yp = diff(y) ./ Dx;
T = zeros(n-2,n-2);
r = zeros(n-2,1);
for i=2:n-3
   T(i,i)   = 2*(Dx(i) + Dx(i+1));
   T(i,i-1) = Dx(i+1);
   T(i,i+1) = Dx(i);
   r(i)     = 3*(Dx(i+1)*yp(i) + Dx(i)*yp(i+1));
end

% For each of the 3 cases, finish setting up the linear system,
% solve the system, and set s(1:n) to be the vector of slopes.

if nargin==5
   %Derivative information available.
   if derivative==1
      % End values for S'(z) specified.         
      T(1,1) = 2*(Dx(1) + Dx(2));
      T(1,2) = Dx(1);
      r(1) = 3*(Dx(2)*yp(1)+Dx(1)*yp(2)) - Dx(2)*muL;
      T(n-2,n-2) = 2*(Dx(n-2)+Dx(n-1));
      T(n-2,n-3) = Dx(n-1);
      r(n-2) = 3*(Dx(n-1)*yp(n-2) + Dx(n-2)*yp(n-1)) -Dx(n-2)*muR;
      s = [muL; T\r; muR];
   end
   
   if derivative==2
      % End values for S''(z) specified.  
      T(1,1) = 2*Dx(1) + 1.5*Dx(2);
      T(1,2) = Dx(1);
      r(1) = 1.5*Dx(2)*yp(1) + 3*Dx(1)*yp(2) + Dx(1)*Dx(2)*muL/4;
      T(n-2,n-2) = 1.5*Dx(n-2)+2*Dx(n-1);
      T(n-2,n-3) = Dx(n-1);
      r(n-2) = 3*Dx(n-1)*yp(n-2) + 1.5*Dx(n-2)*yp(n-1)-Dx(n-1)*Dx(n-2)*muR/4;
      stilde = T\r;
      s1 = (3*yp(1) - stilde(1) - muL*Dx(1)/2)/2;
      sn = (3*yp(n-1) - stilde(n-2) + muR*Dx(n-1)/2)/2;
      s = [s1;stilde;sn];
   end;
else
   % No derivative information. Compute the not-a-knot spline.
   q = Dx(1)*Dx(1)/Dx(2);
   T(1,1) = 2*Dx(1) +Dx(2) + q;
   T(1,2) = Dx(1) + q;
   r(1) = Dx(2)*yp(1) + Dx(1)*yp(2)+2*yp(2)*(q+Dx(1));
   q = Dx(n-1)*Dx(n-1)/Dx(n-2);
   T(n-2,n-2) = 2*Dx(n-1) + Dx(n-2)+q;
   T(n-2,n-3) = Dx(n-1)+q;
   r(n-2) = Dx(n-1)*yp(n-2) + Dx(n-2)*yp(n-1) +2*yp(n-2)*(Dx(n-1)+q);
   stilde = T\r;
   s1 = -stilde(1)+2*yp(1);
   s1 = s1 + ((Dx(1)/Dx(2))^2)*(stilde(1)+stilde(2)-2*yp(2));
   sn = -stilde(n-2) +2*yp(n-1);
   sn = sn+((Dx(n-1)/Dx(n-2))^2)*(stilde(n-3)+stilde(n-2)-2*yp(n-2));
   s = [s1;stilde;sn];
end

% Compute the a,b,c,d vectors.
   
a = y(1:n-1);
b = s(1:n-1);
c = (yp - s(1:n-1)) ./ Dx;
d = (s(2:n) + s(1:n-1) - 2*yp) ./ (Dx.* Dx);