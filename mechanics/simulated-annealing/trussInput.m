% User-defined inputs for 10 bar structure
Tstart = 1000;
lb = 1E-4 * ones(1,10);
ub = 5E-3 * ones(1,10);
maxiter = 5000;

% best selection of c and epsilon
c = 0.99;
epsilon = 0.3;