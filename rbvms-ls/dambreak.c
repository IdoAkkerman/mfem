//--------------------------------------------------------------
// Solution function for the dambreak problem
//
// To compile run, for example:
// mpicc -shared -o rbvms-ls.so -fPIC dambreak.c
// gcc   -shared -o rbvms-ls.so -fPIC dambreak.c
//--------------------------------------------------------------

#include <math.h>
#include <stdio.h>


void sol_u(double *coord, int dim, double time, double *u, int vdim)
{
   u[0] = 0.0;
   u[1] = 0.0;
}
double sol_p(double *coord, int dim, double time)
{
   return 0.0;
}


double sol_phi(double *coord, int dim, double time)
{
   double x = coord[0] - 0.3;
   double y = coord[1] - 0.7;
   double r = sqrt(x*x + y*y);
   //return fmin(-x,-z);
   return ((x > 0) && (y >0)) ? -r : fmin(-x,-y);
}

void force(double *coord, int dim, double time, double *force, int vdim)
{
   force[0] = 0.0;
   force[1] = -9.81;
}

double mu(double *coord, int dim, double time)
{
   return 0.01;
}

