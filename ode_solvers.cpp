#include "prototypes.h"


//////////////////////////////////////////////////////////////// EulerStep
void EulerStep (double t, double *Y, void (*RHS_Func)(double, double *, double *),
                double dt, int neq){
  int k;
  double rhs[neq];
  RHS_Func (t, Y, rhs);
  for (k = 0; k < neq; k++){
    Y[k] += dt*rhs[k];
  }
}


//////////////////////////////////////////////////////////////// Runge-Kutta2
void RK2Step(double t, double *Y, void (*RHS_Func)(double, double *, double *),
             double h, int neq){
  int n;
  double Y1[neq], k1[neq], k2[neq];
  RHS_Func(t,Y,k1);
  for (n=0; n < neq; n++){
    Y1[n] = Y[n]+0.5*h*k1[n];
  }
  RHS_Func(t+0.5*h,Y1,k2);
  for (n=0; n< neq; n++){
    Y[n] += h*k2[n];
  }
}

//////////////////////////////////////////////////////////////// Runge-Kutta4
void RK4Step(double t, double *Y, void (*RHS_Func)(double, double *, double *),
             double h, int neq){
  int n;
  double Y1[neq], k1[neq], k2[neq], k3[neq], k4[neq];
  RHS_Func(t,Y,k1);
  for (n=0; n < neq; n++){
    Y1[n] = Y[n]+0.5*h*k1[n];
  }
  RHS_Func(t+0.5*h,Y1,k2);
  for (n=0; n < neq; n++){
    Y1[n] = Y[n]+0.5*h*k2[n];
  }
  RHS_Func(t+0.5*h,Y1,k3);
  for (n=0; n< neq; n++){
    Y1[n] = Y[n]+h*k3[n];
  }
  RHS_Func(t+h,Y1,k4);
  for (n=0; n< neq; n++){
    Y[n] = Y[n]+ (h/6.0)*(k1[n] + 2*k2[n] + 2*k3[n] + k4[n]);
  }
}


//////////////////////////////////////////////////////////////// Position-Verlet
void VerletPosition(double *x, double *v, int neq, double h,
                    void (*Acceleration)(double *, double *)){
  int n;
  double a[neq];
                                //steps explained
  for (n=0; n < neq; n++){      //evolve position by half a step [drift]
    x[n] = x[n]+0.5*h*v[n];
  }
  Acceleration(x, a);           //compute acceleration at t = t(n+1/2)
  for (n=0; n<neq; n++){        //evolve velocities by a full step
    v[n] = v[n] +h*a[n];
  }
  for (n=0; n<neq; n++){        //evolve position by half a step
    x[n] = x[n] + 0.5*h*v[n];
  }
}

//////////////////////////////////////////////////////////////// velocity-Verlet
void VerletVelocity(double *x, double *v, int neq, double h,
                    void (*Acceleration)(double *, double *)){
  int n;
  double a[neq];
  Acceleration(x, a);      //compute acceleration at t = tn

  for (n=0; n < neq; n++){ //evolve velocities to half time step
    v[n] = v[n]+0.5*h*a[n];
  }
  for (n=0; n<neq; n++){   //evolve coordinates by a full step
    x[n] = x[n] +h*v[n];
  }
  Acceleration(x, a);
  for (n=0; n<neq; n++){   //evolve velocities to half time step
    v[n] = v[n] + 0.5*h*a[n];
  }
}


//////////////////////////////////////////////////////////////////////////////
// Note that this method requires three evaluations of the force per time step,
// as opposed to just one for leapfrog. Note too that the steps are symmetric
// about the middle one (this ensures time reversal invariance).
// The algorithm gives an error of order h5 for one interval (and hence of
// order h4 when integrated over n = T /h time steps for a fixed time
// increment T )
//////////////////////////////////////////////////////////////// Forest-Ruth
void ForestRuth (double *x, double *v, int neq, double h,
                    void (*Acceleration)(double *, double *))
{
  double gamma = 1.0/(2.0 - pow(2.0,1.0/3.0));
  double m = 1.0;
  int n;
  double a[neq];

  for (n=0; n < neq; n++){
    x[n] = x[n]+ gamma*0.5*h*v[n];
  }
  Acceleration(x, a);
  for (n=0; n<neq; n++){
    v[n] = v[n] + gamma*h*m*a[n];
  }
  for (n=0; n<neq; n++){
    x[n] = x[n] + (1.0 - gamma)*0.5*h*v[n];
  }
  Acceleration(x, a);
  for (n=0; n<neq; n++){
    v[n] = v[n] + (1.0 - 2.0*gamma)*h*m*a[n];
  }
  for (n=0; n<neq; n++){
    x[n] = x[n] + (1.0 - gamma)*0.5*h*v[n];
  }
  Acceleration(x, a);
  for (n=0; n<neq; n++){
    v[n] = v[n] + gamma*h*m*a[n];
  }
  for (n=0; n < neq; n++){
    x[n] = x[n]+ gamma*0.5*h*v[n];
  }

}
