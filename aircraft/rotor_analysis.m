% Rotor Analysis

clc;
clear;
close all
% --- Read Data ---

%Lower surface
[x,y,p,M,Cp] = textread ('NASA-Rotor37.dat','%n %n %n %n %n',97,'headerlines',4);
%upper surface
[x2,y2,p2,M2,Cp2] = textread ('NASA-Rotor37.dat','%n %n %n %n %n',97,'headerlines',104);


% -- Process & Display Results --

% Interpolate (x,y) coordinates using spline function: cubic spline
% interpolant

x_int = 0:0.004:1; %251 equally-spaced points x =0:1

% Lower surface - Display Cubic Spline
ySpline = spline (x,y, x_int);
figure
plot (x_int, ySpline)
xlabel ('x')
ylabel ('y')
title ('Lower Surface Cubic Spline Interpolant')

% Upper surface - Display Cubic Spline
ySpline2 = spline (x2,y2,x_int);
figure
plot (x_int, ySpline2)
xlabel ('x')
ylabel ('y')
title ('Upper Surface Cubic Spline Interpolant')

% Improved predictions of the cubic spline by taking into account path
% length along the rotor surface


% Calculating estimated path length
l1 = zeros(size(x)); %lower surface
l2 = zeros(size(x)); %upper surface
for i = 2:length(x)
    % lower surface
    l1(i) = l1(i - 1) + sqrt ((x(i) - x(i-1))^2 + (y(i) - y(i-1))^2);
    %upper surface
    l2(i) = l2(i - 1) + sqrt ((x2(i) - x2(i-1))^2 + (y2(i) - y2(i-1))^2);
end

%calculating piecewise cubic spline interpolant
%lower surface
%x=s(l); y=s(l)
sx1 = spline (l1,x);
sy1S = spline (l1,y);

%coefficients of equation in form of a polynomial
coeffx1 = sx1.coefs;
coeffy1 = sy1S.coefs;
sy1 = zeros(size(x_int)); %spline of y
sy1(1)= y(1);

for i = 2:length(x_int) %calculate spline y
    n = Locate (x, x_int(i)); % determine which subinterval n
    eqn = coeffx1(n,:);
    eqn(4) = eqn(4) - x_int(i); % constant coefficient accounts for x_int
    rt = roots (eqn);
    sy1(i) = polyval(coeffy1(n,:), rt(3)); %evaluate with 3rd root
end
figure
plot (ppval (sx1,x_int),sy1) %plot interpolant
xlabel ('x')
ylabel ('y')
title ('Lower Surface Cubic Spline Improved Interpolant')

% same procedure for upper surface
sx2 = spline (l2,x2);
sy2S = spline (l2,y2);

%coefficients of equation in form of a polynomial
coeffx2 = sx2.coefs;
coeffy2 = sy2S.coefs;
sy2 = zeros(size(x_int));
sy2(1)= y2(1);

for i = 2:length(x_int)
    n = Locate (x2, x_int(i)); % determine which subinterval n
    eqn = coeffx2(n,:);
    eqn(4) = eqn(4) - x_int(i); % constant coefficient accounts for x_int
    rt = roots (eqn);
    sy2(i) = polyval(coeffy2(n,:), rt(3)); %evaluate
end
figure
plot (ppval (sx2,x_int),sy2)
xlabel ('x')
ylabel ('y')
title ('Upper Surface Cubic Spline Improved Interpolant')

% Spline interpolants for Mach number, pressure coefficient

% MACH NUMBER INTERPOLANT
%lower surface

sm1S = spline (l1,M);

%coefficients of equation in form of a polynomial
coeffm1 = sm1S.coefs;
sm1 = zeros(size(x_int));
sm1(1)= M(1);

for i = 2:length(x_int)
    n = Locate (x, x_int(i)); % determine which subinterval n
    eqn = coeffx1(n,:);
    eqn(4) = eqn(4) - x_int(i); % constant coefficient accounts for x_int
    rt = roots (eqn);
    sm1(i) = polyval(coeffm1(n,:), rt(3)); %evaluate
end
figure
plot (ppval (sx1,x_int),sm1)
xlabel ('x')
ylabel ('M')
title ('Lower Surface Cubic Spline Improved Interpolant (Mach Number)')

%upper surface

sm2S = spline (l2,M2);

%coefficients of equation in form of a polynomial
coeffm2 = sm2S.coefs;
sm2 = zeros(size(x_int));
sm2(1)= M2(1);

for i = 2:length(x_int)
    n = Locate (x2, x_int(i)); % determine which subinterval n
    eqn = coeffx2(n,:);
    eqn(4) = eqn(4) - x_int(i); % constant coefficient accounts for x_int
    rt = roots (eqn);
    sm2(i) = polyval(coeffm2(n,:), rt(3)); %evaluate
end
figure
plot (ppval (sx2,x_int),sm2)
xlabel ('x')
ylabel ('M')
title ('Upper Surface Cubic Spline Improved Interpolant (Mach Number)')




% PRESSURE CONSTANT INTERPOLANT
%lower surface

sc1S = spline (l1,Cp);

%coefficients of equation in form of a polynomial
coeffc1 = sc1S.coefs;
sc1 = zeros(size(x_int));
sc1(1)= Cp(1);

for i = 2:length(x_int)
    n = Locate (x, x_int(i)); % determine which subinterval n
    eqn = coeffx1(n,:);
    eqn(4) = eqn(4) - x_int(i); % constant coefficient accounts for x_int
    rt = roots (eqn);
    sc1(i) = polyval(coeffc1(n,:), rt(3)); %evaluate
end
figure
plot (ppval (sx1,x_int),sc1)
xlabel ('x')
ylabel ('Cp')
title ('Lower Surface Cubic Spline Improved Interpolant (Pressure Constant)')

%upper surface

sc2S = spline (l2,Cp2);

%coefficients of equation in form of a polynomial
coeffc2 = sc2S.coefs;
sc2 = zeros(size(x_int));
sc2(1)= Cp2(1);

for i = 2:length(x_int)
    n = Locate (x2, x_int(i)); % determine which subinterval n
    eqn = coeffx2(n,:);
    eqn(4) = eqn(4) - x_int(i); % constant coefficient accounts for x_int
    rt = roots (eqn);
    sc2(i) = polyval(coeffc2(n,:), rt(3)); %evaluate
end
figure
plot (ppval (sx2,x_int),sc2)
xlabel ('x')
ylabel ('Cp')
title ('Upper Surface Cubic Spline Improved Interpolant (Pressure Constant)')

