  function tmin = Golden(fname,a,b)
% tmin = Golden(fname,a,b)
% Golden Section Search
%
% fname   string that names function  f(t) of a single variable.
% a,b   define an interval [a,b] upon which f is unimodal.
%
% tmin   approximate global minimizer of f on [a,b].

r = (3 - sqrt(5))/2;
c = a + r*(b-a);     fc = feval(fname,c);
d = a + (1-r)*(b-a); fd = feval(fname,d);
while (d-c) > sqrt(eps)*max(abs(c),abs(d))
   if fc >= fd
      z = c + (1-r)*(b-c); 
      % [a c d b ] <--- [c d z b]
      a = c; 
      c = d; fc = fd;
      d = z; fd = feval(fname,z);
   else
      z = a + r*(d-a); 
      % [a c d b ] <--- [a z c d]
      b = d; 
      d = c; fd = fc;
      c = z; fc = feval(fname,z);
   end
end
tmin = (c+d)/2;