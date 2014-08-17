/**
 * Using Newton's method to find local minimum of a given objective
 * function:
 * f(x,y,z) = x^2 + 2*y^2 + z^2 - 2*x*y + 2*x - 5/2*y - z + 2*exp(x*y*z)
 */

#include <cmath>
#include <iostream>

using namespace std;

int main ()
{
  const double tol = 1.0e-6; // tolerance
  double x[3] = {0, 0, 0}; // initial guess (x,y,z)
  double gradf[3]; // gradient of f
  double h[3][3]; // Jacobian (Hessian Matrix for f)
  double deth; // determinant of Jacobian for Cramer's Rule
  double detv1, detv2, detv3; // determinants for Cramer's Rule
  double v[3]; // changes to initial/previous guess
  double xtemp [3]; // comparison for tolerance

  do
  {
    // Calculate gradient of f
    gradf[0] = 2 * x[0] - 2 * x[1] + 2 + 2 * x[1] * x[2] * exp(x[0] * x[1] * x[2]);
    gradf[1] = 4 * x[1] - 2 * x[0] - 2.5 + 2 * x[0] * x[2] * exp(x[0] * x[1] * x[2]);
    gradf[2] = 2 * x[2] - 1 + 2 * x[0] * x[1] * exp(x[0] * x[1] * x[2]);

    // Calculate Hessian/Jacobian matrix
    h[0][0] = 2 + 2 * pow(x[1], 2) * pow(x[2], 2) * exp(x[0] * x[1] * x[2]);
    h[0][1] = -2 + 2 * x[0] * x[1] * pow(x[2], 2) * exp(x[0] * x[1] * x[2]) +
              2 * x[2] * exp(x[0] * x[1] * x[2]);
    h[0][2] = 2 * x[1] * exp(x[0] * x[1] * x[2]) + 2 * x[0] * pow(x[1], 2) *
              x[2] * exp(x[0] * x[1] * x[2]);
    h[1][0] = h[0][1];
    h[1][1] = 4 + 2 * pow(x[0], 2) * pow(x[2], 2) * exp(x[0] * x[1] * x[2]);
    h[1][2] = 2 * x[0] * exp(x[0] * x[1] * x[2]) + 2 * pow(x[0], 2) * x[1]
              * x[2] * exp(x[0] * x[1] * x[2]);
    h[2][0] = h[0][2];
    h[2][1] = h[1][2];
    h[2][2] = 2 + 2 * pow(x[0], 2) * pow(x[1], 2) * exp(x[0] * x[1] * x[2]);

    // Calculate determinant of Jacobian
    deth = h[0][0] * h[1][1] * h[2][2] +
           h[0][1] * h[1][2] * h[2][0] +
           h[0][2] * h[1][0] * h[2][1] -
           h[2][0] * h[1][1] * h[0][2] -
           h[2][1] * h[1][2] * h[0][0] -
           h[2][2] * h[1][0] * h[0][1];

    // Calculate determinants for Cramer's Rule
    detv1 = -gradf[0] * h[1][1] * h[2][2] +
            h[0][1] * h[1][2] * -gradf[2] +
            h[0][2] * -gradf[1] * h[2][1] -
            -gradf[2] * h[1][1] * h[0][2] -
            h[2][1] * h[1][2] * -gradf[0] -
            h[2][2] * -gradf[1] * h[0][1];
    detv2 = h[0][0] * -gradf[1] * h[2][2] +
            -gradf[0] * h[1][2] * h[2][0] +
            h[0][2] * h[1][0] * -gradf[2] -
            h[2][0] * -gradf[1] * h[0][2] -
            -gradf[2] * h[1][2] * h[0][0] -
            h[2][2] * h[1][0] * -gradf[0];
    detv3 = h[0][0] * h[1][1] * -gradf[2] +
            h[0][1] * -gradf[1] * h[2][0] +
            -gradf[0] * h[1][0] * h[2][1] -
            h[2][0] * h[1][1] * -gradf[0] -
            h[2][1] * -gradf[1] * h[0][0] -
            -gradf[2] * h[1][0] * h[0][1];

    // Cramer's Rule solving for v in Jv = -gradf
    v[0] = 1 / deth * detv1;
    v[1] = 1 / deth * detv2;
    v[2] = 1 / deth * detv3;

    // Remember for comparison
    xtemp [0] = x[0];
    xtemp [1] = x[1];
    xtemp [2] = x[2];

    // Update: x(k+1) = x(k) + v(k)
    x[0] += v[0];
    x[1] += v[1];
    x[2] += v[2];
  }
  // Iterate until tolerance is satisfied
  while ( abs(xtemp[0] - x[0]) >= tol ||
          abs(xtemp[1] - x[1]) >= tol ||
          abs(xtemp[2] - x[2]) >= tol );

  // Output approximate solution to nonlinear equation
  cout << "Solution:" << endl << "x = " << x[0] << endl
       << "y = " << x[1] << endl << "z = " << x[2] << endl;
}
