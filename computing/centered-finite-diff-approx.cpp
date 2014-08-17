/**
 * Using Second-order centred finite-difference approximation
 * to solve for PDE: du/dt = d2u/dx2 - 6x
 */
 
 #include <cmath>
 #include <iostream>
 
using namespace std;

// Global variables
int neqn;
double err;
double xl, xu, dx, dxs;

void intpar(); // integration parameters
void par (); // problem parameters
double f (double uf, double u, double ub, double dxs, double g); // evaluate derivative
double rk4 (double uf, double u, double ub, double dxs, double g,
              double dt); // evaluate runge kutta 4

/**
 * Set parameters to control integration of PDE system
 */
void intpar()
{
  // Number of first order ODEs
  neqn = 41;
  err = 1.0e-05;
}

/**
 * Set variables for PDE system
 */
void par()
{
  xl = 0.0;
  xu = 1.0;
  dx = (xu - xl)/(neqn - 1);
  dxs = dx * dx;
}

/** f (t, U)
 * Evaluate ODE with parameters t and u
 */
double f (double uf, double u, double ub, double dxs, double g)
{
  return (uf - 2.0 * u + ub) / dxs - g;
}

/** rk4 (t, u dt)
 * Fourth order Runge Kutta with parameters t and u and timestep dt
 */
double rk4 (double uf, double u, double ub, double dxs, double g, double dt)
{
  double f1 = f(uf, u, ub, dxs, g);
  double f2 = f(uf, u + dt * f1 / 2.0, ub, dxs, g);
  double f3 = f(uf, u + dt * f2 / 2.0, ub, dxs, g);
  double f4 = f(uf, u + dt * f3, ub, dxs, g);
  return u + dt / 6.0 * (f1 + 2.0 * f2 + 2.0 * f3 + f4);
}


// Main Program
int main ()
{
  intpar(); // integration parameters
  par(); // PDE parameters
  double u[neqn]; // initial conditions
  double ut[neqn]; // system of ODE
  
  
  // Set initial conditions for PDE system
  double x;
  for (int i = 1; i <= neqn; i++)
  {
    x = xl + (float)(i - 1) / (neqn - 1) * (xu - xl);
    u[i - 1] = x;
  }
  
  // Boundary Conditions
  ut[1 - 1] = 0.0; // BC @ x = xl
  ut[neqn - 1] = 1.0; // BC @ x = xu
  
  // Evaluate derivative for PDE system
  for (int i = 2; i <= (neqn - 1); i++)
  {
    x = xl + (float)(i - 1) / (neqn - 1) * (xu - xl);
    // ODE approximation for PDE
    ut[i - 1] = (u[i + 1 - 1] - 2.0 * u[i - 1] + u[i - 1 - 1]) / dxs
                - 6 * x;
  }
  
  double t = 0.0; // t0
  double dt = 1.0/1000.0; // time step
  double diff [neqn - 1]; // difference for error tolerance
  
  for (int i = 2; i < neqn; i++)
  {
    diff[i - 2] = err + 1; // arbitrary num greater than err
  }
  
  bool done = false; // true if PDE system reach steady state
  int iter = 0;
  // 
  while (true)
  {
    for (int i = 2;i<neqn;i++)
    {
      if (diff[i-2] > err)
        break;
      done = true;
    }
    if (done)
      break;
      
    for (int i = 2; i <= (neqn - 1); i++)
    {
      double before = ut[i-1];
      x = xl + (float)(i - 1) / (neqn - 1) * (xu - xl);
      ut[i-1]=rk4(u[i + 1-1],u[i-1],u[i - 1-1],dxs, -6*x, dt);
      t = t + dt;
      double after = ut[i-1];
      diff[i-2] = abs (before - after);
    }
    iter++;
  }
  
  for (int i = 0; i<neqn; i++)
  {
    cout<<ut[i]<<endl;
  }
  cout<<iter<<endl;
  
  
  
}
 