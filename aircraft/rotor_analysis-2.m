clc;
clear;
close all

% File input using NASA Rotor 37 surface data at 70% span
[xl,yl,pl,Ml,Cpl] = ... % lower surface data
    textread ('NASA-Rotor37.dat','%n %n %n %n %n',97,'headerlines',4);
[xu,yu,pu,Mu,Cpu] = ... % upper surface data
    textread ('NASA-Rotor37.dat','%n %n %n %n %n',97,'headerlines',104);


% Integrates the pressure coefficient on the upper and lower surfaces in
% order to determine the lift coefficient for the lift forces based on
% piecewise cubic spline interpolant.

% Generate pressure coefficient spline interpolants for lower and upper
% surfaces
x_int_small = 0:0.001:1; % small equally-spaced points x = 0:1
Cpls = spline (xl, Cpl, x_int_small);
Cpus = spline (xu, Cpu, x_int_small);

% Integrate interpolants using spline quadrature to find lift coefficient
ICpls = SplineQ (x_int_small, Cpls);
ICpus = SplineQ (x_int_small, Cpus);
% Lift Coefficient = integral(0->1) [Cpl-Cpu] d(x/c)
liftCoeff = ICpls - ICpus;

fprintf ('The lift coefficient using spline quadrature is %f.\n\n', liftCoeff); % display result


% Determines the pressure drag coefficient

% Calculate y(x)
yls = spline (xl, yl);
yus = spline (xu, yu);
% Copy of interpolant for future modification (i.e. change to interpolant
% dy/dx)
Dyls = yls;
Dyus = yus;

% Find the coefficients of dy/dx for both lower and upper surfaces
dyls = zeros (length (yls.coefs),4);
dyus = zeros (length (yls.coefs),4);
for j = 1:length (yls.coefs)
    dyls(j,2:4) = polyder (yls.coefs(j,:)); % Coefficients are stored
    dyus(j,2:4) = polyder (yus.coefs(j,:)); % in a 1x4 matrix
end

% A new interpolant with dy/dx coefficients
Dyls.coefs = dyls;
Dyus.coefs = dyus;

% Pressure Drag Coefficient = integral (0->1) [(Cpl-Cpu)*(dy/dx)] dx
% Integrate using spline quadrature
Ipdls = SplineQ (x_int_small, Cpls .* ppval(Dyls, x_int_small));
Ipdus = SplineQ (x_int_small, Cpus .* ppval(Dyus, x_int_small));
pDragCoeff = Ipdls - Ipdus;

% display result
fprintf ('The pressure drag coefficient using spline quadrature is %f.\n\n', pDragCoeff); 



% Uses cuplic spline interpolant for surface coordinates, pressure
% coefficient as function of path length to integrate the pressure forces
% acting on the rotor in order to determine the lift coefficient, pressure
% drag coefficient, and pitching moment coefficient.

% Calculate path length
ll = zeros(size(xl)); % lower surface
lu = zeros(size(xu)); % upper surface
for i = 2:length(xl)
    % lower surface path length
    ll(i) = ll(i - 1) + sqrt ((xl(i) - xl(i-1))^2 + (yl(i) - yl(i-1))^2);
    % upper surface path length
    lu(i) = lu(i - 1) + sqrt ((xu(i) - xu(i-1))^2 + (yu(i) - yu(i-1))^2);
end

% Cubic spline interpolant for surface coordinates, pressure coefficient as
% function of path length

% lower surface
xlsl = spline (ll, xl);
ylsl = spline (ll, yl);
Cplsl = spline (ll, Cpl);

% upper surface
xusl = spline (lu, xu);
yusl = spline (lu, yu);
Cpusl = spline (lu, Cpu);

% y-coordinates and pressure coefficient as function of x-coordinates for
% both lower and upper surfaces
yxl = zeros(size(x_int_small));
yxu = zeros(size(x_int_small));
Cpxl = zeros(size(x_int_small));
Cpxu = zeros(size(x_int_small));
yxl(1) = yl(1);
yxu(1) = yu(1);
Cpxl(1) = Cpl(1);
Cpxu(1) = Cpu(1);

% Calculation of y-coordinates and pressure coefficient
for i = 2:length(x_int_small)
    % Lower surface y-coordinates and pressure coefficient
    j = Locate (xl, x_int_small(i)); % determine subinterval j
    eqn = xlsl.coefs(j,:); % coefficients of cubic spline polynomial
    eqn(4) = eqn(4) - x_int_small(i); % constant coefficient accounts for x_int_small (i)
    root = roots (eqn); % solve for roots
    yxl(i) = polyval(ylsl.coefs(j,:),root(3)); % evaluate Y-coordinate
    Cpxl(i) = polyval(Cplsl.coefs(j,:),root(3)); % evaluate Pressure Coefficient
    
    % Upper surface y-coordinates and pressure coefficient
    j = Locate (xu, x_int_small(i)); % determine subinterval j
    eqn = xusl.coefs(j,:); % coefficients of cubic spline polynomial
    eqn(4) = eqn(4) - x_int_small(i); % constant coefficient accounts for x_int_small (i)
    root = roots (eqn); % solve for roots
    yxu(i) = polyval(yusl.coefs(j,:),root(3)); % evaluate Y-coordinate
    Cpxu(i) = polyval(Cpusl.coefs(j,:),root(3)); % evaluate Pressure Coefficient
end

% Integration interval
x_int = 0:0.01:1;

% Using spline function to map x_int_small with the calculated interpolant
% Purpose: e.g. maps yxl(i) and x_int_small(i) together into yxl(x_int_small)
mapyxl = spline (x_int_small, yxl);
mapyxu = spline (x_int_small, yxu);
mapCpxl = spline (x_int_small, Cpxl);
mapCpxu = spline (x_int_small, Cpxu);

% Copy of interpolant for future modification (i.e. change to interpolant
% dy/dx)
Dmapyxl = mapyxl;
Dmapyxu = mapyxu;

% Find the coefficients of dy/dx for both lower and upper surfaces
dmapyxl = zeros (length (mapyxl.coefs),4);
dmapyxu = zeros (length (mapyxu.coefs),4);
for j = 1:length (mapyxl.coefs)
    dmapyxl(j,2:4) = polyder (mapyxl.coefs(j,:)); % Coefficients are stored
    dmapyxu(j,2:4) = polyder (mapyxu.coefs(j,:)); % in a 1x4 matrix
end

% A new interpolant with dy/dx coefficients
Dmapyxl.coefs = dmapyxl;
Dmapyxu.coefs = dmapyxu;


% Create function handles (i.e. integrands)
CpxlF = @(x) ppval (mapCpxl,x); % integrand for lift coef
CpxuF = @(x) ppval (mapCpxu,x); % "
pdlF = @(x) ppval (mapCpxl,x) .* ppval (Dmapyxl,x); % integrand for ...
pduF = @(x) ppval (mapCpxu,x) .* ppval (Dmapyxu,x); % pressure drag coef
CpxlFmom = @(x) x .* ppval (mapCpxl,x); % integrand for ...
CpxuFmom = @(x) x .* ppval (mapCpxu,x); % pitching moment coef

for m = [2 3] % m-point rule used for Newton-Cotes
   for n = [length(x_int) length(x_int_small)] % subintervals
      % Calculate Newton-Cotes integration for both lower and upper
      % surfaces
      ICpl = CompQNC (CpxlF, x_int(1), x_int(length(x_int)), m, n); 
      ICpu = CompQNC (CpxuF, x_int(1), x_int(length(x_int)), m, n);
      Ipdl = CompQNC (pdlF, x_int(1), x_int(length(x_int)), m, n);
      Ipdu = CompQNC (pduF, x_int(1), x_int(length(x_int)), m, n);
      ICpl25 = CompQNC (CpxlFmom, x_int(1), 0.25 * x_int(length(x_int)), m, n); % 25% of the rotor
      ICpl75 = CompQNC (CpxlFmom, 0.25 * x_int(length(x_int)), x_int(length(x_int)), m, n); % last 75% of the rotor
      ICpu25 = CompQNC (CpxuFmom, x_int(1), 0.25 * x_int(length(x_int)), m, n); % 25% of the rotor
      ICpu75 = CompQNC (CpxuFmom, 0.25 * x_int(length(x_int)), x_int(length(x_int)), m, n); % last 75% of the rotor
      
      % Lift Coefficient = integral(0->1) [Cpl-Cpu]d(x/c)
      liftCoeffQnc = ICpl - ICpu;
      % Pressure Drag Coefficient = integral(0->1) [(Cpl-Cpu) * (dy/dx)] dx
      pDragCoeffQnc = Ipdl - Ipdu;
      % Pitching moment = integral(0->0.25) [(x/c)(Cpl-Cpu)d(x/c)] - integral(0.25->0.75) [(x/c)(Cpl-Cpu)d(x/c)]
      pitchMomQnc = (ICpl25 - ICpu25) - (ICpl75 - ICpu75);
      
      % Display Results
      if m == 2
         fprintf ('The lift coefficient using trapezoidal quadrature and %.0f subintervals is %f.\n', n - 1, liftCoeffQnc);
         fprintf ('The pressure drag coefficient using trapezoidal quadrature and %.0f subintervals \nis %f.\n', n - 1, pDragCoeffQnc);
         fprintf ('The pitching moment coefficient about the quarter-chord using trapezoidal quadrature\nand %.0f subintervals is %f.\n', n - 1, pitchMomQnc);
      else
         fprintf ('The lift coefficient using Simpsons quadrature and %.0f subintervals is %f.\n', n - 1, liftCoeffQnc);
         fprintf ('The pressure drag coefficient using Simpsons quadrature and %.0f subintervals \nis %f.\n', n - 1, pDragCoeffQnc);
         fprintf ('The pitching moment coefficient about the quarter-chord using Simpsons quadrature\nand %.0f subintervals is %f.\n', n - 1, pitchMomQnc);
      end
      % Compare estimates with Part 1 & 2
      fprintf ('[ESTIMATE COMPARISON]\n')
      fprintf ('The error of lift coefficient compared to Part 1 is %f.\n', abs (liftCoeffQnc - liftCoeff))
      fprintf ('Error of lower surface pressure coefficient (compared to P1): %f\n', abs (ICpl - ICpls))
      fprintf ('Error of upper surface pressure coefficient (compared to P1): %f\n', abs (ICpu - ICpus))
      fprintf ('The error of pressure drag coefficient compared to Part 2 is %f.\n\n', abs (pDragCoeffQnc - pDragCoeff))
   end
end
