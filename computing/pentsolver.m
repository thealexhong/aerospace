function x = pentsolver(A, b)
% Solves a general pentadiagonal linear system of equations (modified
% Thomas algorithm)

    % isolate diagonals of A
    d = diag(A);
    c = diag(A,1);
    a = diag(A,2);
    n = length(d);

    % Factor A into 3 matrices: u, l, q
    u = zeros (n, 1);
    l = zeros (n-1, 1);
    q = zeros (n-2, 1);
    
    % Calculating u, l, q
    u(1) = d(1);
    l(1) = c(1) / u(1);
    q(1) = a(1) / u(1);
    
    u(2) = d(2) - c(1) * l(1);
    l(2) = (c(2) - a(1) * l(1)) / u(2);
    q(2) = a(2) / u(2);
    
    for i = 3:n-2
        u(i) = d(i) - a(i-2) * q(i-2) - u(i-1) * l(i-1)^2;
        l(i) = (c(i) - a(i-1) * l(i-1)) / u(i);
        q(i) = a(i) / u(i);
    end
    
    u(n-1) = d(n-1) - a(n-3) * q(n-3) - u(n-2) * l(n-2)^2;
    l(n-1) = (c(n-1) - a(n-2) * l(n-2)) / u(n-1);
    u(n) = d(n) - a(n-2) * q(n-2) - u(n-1) * l(n-1)^2;
    
    % Forward Substitution Lz = b
    z = zeros (n, 1);
    z(1) = b(1);
    z(2) = b(2) - l(1) * z(1);
    for i=3:n
        z(i) = b(i) - l(i-1) * z(i-1) - q(i-2) * z(i-2);
    end
    
    
    % Back substitution Ux = z
    t = z./u;
    x = zeros(n,1);
    x(n) = t(n);
    x(n-1) = t(n-1) - l(n-1) * x(n);
    for i = n-2:-1:1
        x(i) = t(i) - l(i) * x(i+1) - q(i) * x(i+2);
    end

end

