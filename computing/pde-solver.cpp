/**
 * Determines numerical approximation to steady state solution of PDE:
 * du/dt = d^2u/dx^2 - 6x
 * Domain: 0 <= x <= 1
 * u(x=0,t)=0
 * u(x=1,t)=1
 * u(x,t=0)=x
 * Using Point Jacobi, Gauss Seidel Relaxation, and SOR Methods
 */

#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

// Global variables
int neqn;
double err;
double xl, xu, dx, dxs;

/**
 * Set parameters to control integration of PDE system
 */
void intpar()
{
    // Number of first order ODEs
    neqn = 41;
    // error tolerance
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

// Main Program
int main ()
{
    intpar(); // integration parameters
    par();   // PDE parameters
    double u [neqn]; // initial conditions
    double ut [neqn]; //system of ODEs
    double j [neqn]; // point-jacobi solution
    double gs [neqn]; // Gauss Seidel Solution
    double sor [neqn]; // successive over relaxation method

    // Set initial conditions for PDE system
    double x;
    for (int i = 1; i <= neqn; i++)
    {
        x = xl + (float)(i - 1) / (neqn - 1) * (xu - xl);
        u[i - 1] = x;
    }

    // Boundary Conditions (minus 1 as arrays start at 0)
    ut[1 - 1] = 0.0; // BC @ x = xl
    ut[neqn - 1] = 1.0; // BC @ x = xu
    j[0]= 0.0;
    j[neqn-1]=1.0;
    gs[0]=0.0;
    gs[neqn-1]=1.0;
    sor[0]=0.0;
    sor[neqn-1]=1.0;


    // Evaluate derivative for PDE system
    for (int i = 2; i <= (neqn - 1); i++)
    {
        x = xl + (float)(i - 1) / (neqn - 1) * (xu - xl);
        // ODE approximation for PDE (second order centred FDA)
        ut[i - 1] = (u[i + 1 - 1] - 2.0 * u[i - 1] + u[i - 1 - 1]) / dxs
                    - 6 * x;
        // Point Jacobi Method
        j[i-1]=1/2*(j[i-2]+j[i]-ut[i-1]);

        // Gauss Seidel Method
        gs[i-1] = 1/2*(gs[i-2]+gs[i]-ut[i-1]);

        // Successive Over Relaxation Method
        double w = 1.85; // optimal relaxation parameter for M = 39
        double utilda = 1/2 * (ut[i-2]+ut[i]-ut[i-1]);
        sor[i-1] = sor[i-1]+w*(utilda-sor[i-1]);

    }


    // Output
    ofstream myfile;
    myfile.open ("out.txt");

    myfile<<"Point Jacobi Method"<<endl;
    for (int i = 0; i < neqn; i++)
    {
        myfile << j[i] << endl;
    }
    myfile<<"\nGuass seidel"<<endl;
    for (int i = 0; i < neqn; i++)
    {
        myfile << gs[i] << endl;
    }
    myfile<<"\nSuccesive over relaxation"<<endl;
    for (int i = 0; i < neqn; i++)
    {
        myfile << sor[i] << endl;
    }
    myfile.close();
}
