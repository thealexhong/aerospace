/**
 * Using Newton's method to solve a given
 * nonlinear system of equations
 * (1) 3*x1 - cos(x2*x3) - 0.5 = 0
 * (2) x1^2 - 81*(x2 + 0.1)^2 + sin(x3) + 1.06 = 0
 * (3) exp (-x1*x2) + 20*x3 + (10*PI - 3) / 3
 */

#include <cmath>
#include <iostream>
#define PI 4*atan(1)

using namespace std;

int main()
{
  const double tol = 1.0e-6; // tolerance
  double x[3] = {0.10, 0.10, -0.10}; // initial guess
  double f[3]; // nonlinear system of equations
  double fp [3][3]; // Jacobian
  double detfp; // determinant of Jacobian
  double detv1, detv2, detv3; // determinants for Cramer's Rule
  double v[3]; // changes to initial guess
  double xtemp[3]; // comparison for tolerance

  do
  {
    // Calculate nonlinear systems of equations; Matrix F
    f[0] = 3 * x[0] - cos(x[1] * x[2]) - 0.5;
    f[1] = pow(x[0], 2) - 81 * pow(x[1] + 0.10, 2)
            + sin(x[2]) + 1.06;
    f[2] = exp(-x[0] * x[1]) + 20 * x[2]
            + (10 * PI - 3) / 3;
    // Calculate Jacobian elements; Matrix J
    fp [0][0] = 3;
    fp [0][1] = x[2] * sin(x[1] * x[2]);
    fp [0][2] = x[1] * sin(x[1] * x[2]);
    fp [1][0] = 2 * x[0];
    fp [1][1] = -162 * x[1] - 16.2;
    fp [1][2] = cos(x[2]);
    fp [2][0] = -x[1] * exp(-x[0] * x[1]);
    fp [2][1] = -x[0] * exp(-x[0] * x[1]);
    fp [2][2] = 20;
    // Calculate determinant of Jacobian for Cramer's Rule
    detfp = fp[0][0] * fp[1][1] * fp[2][2] +
            fp[0][1] * fp[1][2] * fp[2][0] +
            fp[0][2] * fp[1][0] * fp[2][1] -
            fp[2][0] * fp[1][1] * fp[0][2] -
            fp[2][1] * fp[1][2] * fp[0][0] -
            fp[2][2] * fp[1][0] * fp[0][1];
    // Calculate determinants for Cramer's Rule
    detv1 = -f[0] * fp[1][1] * fp[2][2] +
            fp[0][1] * fp[1][2] * -f[2] +
            fp[0][2] * -f[1] * fp[2][1] -
            -f[2] * fp[1][1] * fp[0][2] -
            fp[2][1] * fp[1][2] * -f[0] -
            fp[2][2] * -f[1] * fp[0][1];
    detv2 = fp[0][0] * -f[1] * fp[2][2] +
            -f[0] * fp[1][2] * fp[2][0] +
            fp[0][2] * fp[1][0] * -f[2] -
            fp[2][0] * -f[1] * fp[0][2] -
            -f[2] * fp[1][2] * fp[0][0] -
            fp[2][2] * fp[1][0] * -f[0];
    detv3 = fp[0][0] * fp[1][1] * -f[2] +
            fp[0][1] * -f[1] * fp[2][0] +
            -f[0] * fp[1][0] * fp[2][1] -
            fp[2][0] * fp[1][1] * -f[0] -
            fp[2][1] * -f[1] * fp[0][0] -
            -f[2] * fp[1][0] * fp[0][1];

    // Cramer's Rule solving for v in Jv = -F
    v[0] = 1 / detfp * detv1;
    v[1] = 1 / detfp * detv2;
    v[2] = 1 / detfp * detv3;

    // Remember for comparison
    xtemp [0] = x[0];
    xtemp [1] = x[1];
    xtemp [2] = x[2];

    // Update: x(k+1) = x(k) + v(k)
    x[0] += v[0];
    x[1] += v[1];
    x[2] += v[2];
  }
  // iterate until tolerance is satisfied
  while ( abs(xtemp[0] - x[0]) >= tol ||
          abs(xtemp[1] - x[1]) >= tol ||
          abs(xtemp[2] - x[2]) >= tol );

  // Output approximate solution to nonlinear equation
  cout << "Solution:" << endl << "x_1 = " << x[0] << endl
       << "x_2 = " << x[1] << endl << "x_3 = " << x[2] << endl;
}
