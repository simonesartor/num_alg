#include <iostream>
#include <cmath>
#include <math.h>
#include <iomanip>
#include <ctime>
#include <fstream>

using namespace std;

//root-finder methods functions prototypes
double bisection(double (*func)(double), double, double, double, double, double &, int &);
double falseposition(double (*func)(double), double, double, double, double, double &, int &);
double secant(double (*func)(double), double, double, double, double, int &);
double newton(double (*func)(double), double (*dfunc)(double), double, double, double, double, int &);
void bracket(double (*func)(double), double, double, int, double*, double*, int &);


// ode solvers prototypes
void EulerStep (double t, double *Y, void (*RHS_Func)(double, double *,
                double *), double dt, int neq);
void RK2Step(double t, double *Y, void (*RHS_Func)(double, double *, double *),
             double h, int neq);
void RK4Step(double t, double *Y, void (*RHS_Func)(double, double *, double *),
             double h, int neq);
void VerletPosition(double *xx, double *vv, int neq, double dt,
                    void (*Acceleration)(double *, double *));
void VerletVelocity(double *xx, double *vv, int neq, double dt,
                    void (*Acceleration)(double *, double *));
void ForestRuth (double *x, double *v, int neq, double h,
                    void (*Acceleration)(double *, double *))
