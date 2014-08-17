/**
 * RK4 to solve ODE: du/dt = -<lamda>*u+a*sin(wt)
 */


#include <cmath>
#include <iostream>

using namespace std;

// constants
const double lam = 3.0 / 2.0;
const double ome = 10;
const double a = 20;
const double u0 = 10;

// functions
double f (double t, double u); // evaluate derivative
double rk4 (double t, double u, double dt); // evaluate runge kutta 4


/** f (t, U)
 * Evaluate given ODE with parameters t and u
 */
double f (double t, double u)
{
  return -lam * u + a * sin(ome * t);
}

/** rk4 (t, u dt)
 * Fourth order Runge Kutta with parameters t and u and timestep dt
 */
double rk4 (double t, double u, double dt)
{
  double f1 = f(t, u);
  double f2 = f(t + dt / 2.0, u + dt * f1 / 2.0);
  double f3 = f(t + dt / 2.0, u + dt * f2 / 2.0);
  double f4 = f(t + dt, u + dt * f3);
  return u + dt / 6.0 * (f1 + 2.0 * f2 + 2.0 * f3 + f4);
}

// Main Program
int main()
{
  double dt = 5.0 / 1000.; // time step
  double u = u0; // initial conditions
  double t = 0.0; // t0 = 0

  // RK4 until t = 5 with dt defined
  do
  {
    u = rk4 (t, u, dt);
    t = t + dt;
  }
  while (t < 5);

  /** Exact Solution
   * u(t) = 10/409 *(489*exp(-3*t/2)+12*sin(10*t)-80*cos(10*t))
   */
  double soln = 10. / 409. * (489. * exp (-3. * 5. / 2.) + 12. * sin (10. * 5.)
                - 80. * cos (10. * 5.));
  double error = abs (soln - u);

  cout << "Question 1\n\nThe solution at t = 5 is " << u
  << endl << "Error (1000 steps): " << error << endl << endl <<
  "Alexander Hong 997584706" << endl;
}
