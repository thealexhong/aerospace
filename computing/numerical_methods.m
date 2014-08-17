% Comparing different numerical methods to solve spring/mass and beam problems.

clc;
clear;
close all;

% mass spring system: solving for displacement using Thomas algorithm

kconst = 10000;
Fconst = 200;
k = [kconst kconst kconst kconst kconst]; % Stiffness (same springs)
F = [Fconst Fconst 100*Fconst Fconst]; % external forces (b matrix)

% Tridiagonals
a = [0 -k(2) -k(3) -k(4)];
d = [k(1)+k(2) k(2)+k(3) k(3)+k(4) k(4)+k(5)];
c = [-k(2) -k(3) -k(4)];

n = length (d); % number of rows
l = zeros (1, n); % lower triangular matrix
u = zeros (1, n); % upper triangular matrix

% Thomas Algorithm
u(1) = d(1);
for i = 2:n
    l(i) = a(i) / u(i-1);
    u(i) = d(i) - l(i) * c(i-1);
end

% Forward substitution
z = zeros (1, n);
z(1) = F(1);
for i = 2:n
    z(i) = F(i) - l(i) * z(i-1);
end

% Back substitution
x = zeros (n, 1); % spring displacements
x(n) = z(n) / u(n);
for i = n-1:-1:1
    x(i) = (z(i)-c(i)*x(i+1))/u(i);
end

% Confirming solution (uncomment to check)
% A = [d(1) c(1) 0 0; a(2) d(2) c(2) 0; 0 a(3) d(3) c(3); 0 0 a(4) d(4)];
% A*x

disp (sprintf ('The displacements, x_1, x_2, x_3, x_4, are'));
disp (x);


% cantilevered beam: solving for deflection using modified Thomas algorithm

clear;
% different values for n
for n = [5 10 15]

    % Setting up the A matrix
    a = ones (1, n - 2); % first upper/lower diagonal (Symmetric)
    c = -4 .* ones (1, n - 1); % second upper/lower diagonal (Symmetric)
    c (n-1) = -2;
    d = 6 .* ones (1, n); % diagonal
    d(1) = 9;
    d(n-1) = 5;
    d(n) = 1;

    % Construct A matrix
    A = diag (d);
    for i = 1:1:n-1
        A (i, i+1) = c(i);
        A (i+1, i) = c(i);
    end
    for i = 1:1:n-2
        A(i,i+2)= a(i);
        A(i+2,i)= a(i);
    end
    
    % Construct b, x matrices
    b = ones (n, 1); % load on the bar (b matrix)
    x = zeros (n, 1); % deflection (x matrix) Ax = b
        
    x = pentsolver (A,b); % solves the pentadiagonal linear system
      
    % Results
    disp (sprintf...
        ('For n = %.f:\nThe deflection of the bar for %.fx%.f A matrix is (solved with \nmodified Thomas algorithm)',n,n,n));
    disp (x);
    
    X = A\b; % blackslash command
    disp ('The same A matrix solved with backslash command gives');
    disp (X);
    
    
    disp ('The same A matrix solved with lu function gives')
    [L U] = lu (A);
    Z = L\b;
    Xlu = U\Z;
    disp (Xlu);
    
    fprintf(1,'Comparing the two methods with the written function, \nthe error for each of the deflection component is\n');
    % error between the two approaches
    fprintf (1, '\t\tBackslash \t\t\t     lu function\n');
    fprintf (1, '\t\t%d \t\t\t %d\n', abs(X-x), abs(Xlu-x));
    % Answers agree with each other well as the error is relatively small
    % compared to its value.
        
    disp (sprintf ('\nThe condition number of A for n = %.f is %.2f.\n',n,cond(A)));
end
disp ('From calculations, the lu function has a larger error than the backslash method.');
disp('One can see that the condition number increases as n gets large.')
disp ('This means the error increases as n increases.');



% calculating the processor time required to solve beam deflections as a
% function of system size, n.

clear;

S = 10; % initial n value
I = 10; % incremental value
N = 1000; % max n value

% Elapsed time for each operation
Tpenta = zeros(1,(N-S)/I);
Tbacks = zeros(1,(N-S)/I);
Tlu = zeros(1,(N-S)/I);
nvar = S:I:N; % on the x-axis of the graph

for n = nvar
    % Setting up the A matrix
    a = ones (1, n - 2); % first upper/lower diagonal (Symmetric)
    c = -4 .* ones (1, n - 1); % second upper/lower diagonal (Symmetric)
    c (n-1) = -2;
    d = 6 .* ones (1, n); % diagonal
    d(1) = 9;
    d(n-1) = 5;
    d(n) = 1;

    % Construct A matrix
    A = diag (d);
    for i = 1:1:n-1
        A (i, i+1) = c(i);
        A (i+1, i) = c(i);
    end
    for i = 1:1:n-2
        A(i,i+2)= a(i);
        A(i+2,i)= a(i);
    end
    
    % Construct b, x matrices
    b = ones (n, 1); % load on the bar (b matrix)
    x = zeros (n, 1); % deflection (x matrix) Ax = b
    
    % timing pentadiagonal solver
    tstart_penta = tic;
    x = pentsolver (A,b); % solves the pentadiagonal linear system
    Tpenta((n-S)/I+1) = toc(tstart_penta);
    
    % timing backslash solving
    tstart_backs = tic;
    X = A\b;
    Tbacks((n-S)/I+1) = toc(tstart_backs);
    
    % timing lu command solving
    tstart_lu = tic;
    [L U] = lu (A);
    Z = L\b;
    Xlu = U\Z;
    Tlu((n-S)/I+1) = toc(tstart_lu);
end

% Graph
figure;
loglog(nvar, Tpenta, nvar, Tbacks,':', nvar, Tlu,'--');
xlabel ('system size: n');
ylabel ('CPU time (s)')
title ('Processing times for different approaches to solve pentadiagonal system -- Alexander Hong 997584706');
legend ('Pentadiagonal Solver','Backslash','lu command');
fprintf (1, 'Looking at the graph, one can see the most efficient (fastest) solver is the \npentadiagonal solver.\n');
disp ('One can also see that backslash command is slightly faster than the lu command.');
disp ('It is recommended to use the pentadiagonal solver for large n.');

