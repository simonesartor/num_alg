#include "prototypes.h"


///////////////////////////////////////////////////////////////////////////
//RectangularRule
double RectangularRule (double (*func)(double), double a, double b, int N){

  double sum = 0.0;
  double dx = (b - a)/(double)N; //lenght of intervals
  int i;

  for (i=0; i<N; i++){
    sum += func(a + i * dx);
  }

  return sum*dx;
}

///////////////////////////////////////////////////////////////////////////
//MidpointRule
double MidpointRule (double (*func)(double), double a, double b, int N){
  double sum = 0.0;
  double dx = (b - a)/(double)N; //lenght of intervals
  int i;

  for (i=0; i<N; i++){
    sum += func(a + i * dx + 0.5*dx);
  }

  return sum*dx;
}

///////////////////////////////////////////////////////////////////////////
// TrapezoidalRule
double TrapezoidalRule (double (*func)(double), double a, double b, int N){
  double sum = 0.0;
  double dx = (b - a)/(double)N; //lenght of intervals
  int i;

  for (i=0; i<N; i++){
    sum += 0.5*dx*(func(a + i * dx)+func(a + (i+1)*dx));
  }

  return sum;
}

///////////////////////////////////////////////////////////////////////////
// SimpsonRule
//
// ! the total number of intervals must be even
///////////////////////////////////////////////////////////////////////////
double SimpsonRule (double (*func)(double), double a, double b, int N){
  double sum;
  double dx = (b - a)/(double)N; //lenght of intervals
  int i;

  sum = func(a) + func(b); //initialize sum (function evaluated at endpoints)
  for (i=1; i<N; i++){
    if (i%2==0){
      sum += 2.0*func(a + i * dx);
    } else{
      sum += 4.0*func(a + i * dx);
    }
  }
  sum = sum * dx/3.0;

  return sum;
}

///////////////////////////////////////////////////////////////////////////
// Gaussian Quadrature
//
// ! here defined up to gn = 8
///////////////////////////////////////////////////////////////////////////
double GaussianQuad (double (*func)(double), double a, double b, int N, int gn){
  double sum, sum_sub;
  double mid;
  double dx = (b - a)/(double)N; // step size
  int i,k;
  double ga[8];
  double gw[8];

  // define gaussian weights gw and abscissae ga depending on gn (gauss number)
  if (gn == 2){
    ga[0] = -1.0/sqrt(3.0); gw[0] = 1.0;
    ga[1] =  1.0/sqrt(3.0); gw[1] = 1.0;
  } else if (gn == 3){
    ga[0] = -sqrt(3.0/5.0); gw[0] = 5.0/9.0;
    ga[1] = 0.0;            gw[1] = 8.0/9.0;
    ga[2] = sqrt(3.0/5.0);  gw[2] = 5.0/9.0;
  } else if (gn == 4){
    double ga4 = 2.0/7.0*sqrt(6.0/5.0);
    double gw4 = (18.0 + sqrt(30.0));
    ga[0] = -sqrt(3.0/7.0 - ga4); gw[0] = gw4/36.0;
    ga[1] =  sqrt(3.0/7.0 - ga4); gw[1] = gw4/36.0;
    ga[2] = -sqrt(3.0/7.0 + ga4); gw[2] = gw4/36.0;
    ga[3] =  sqrt(3.0/7.0 + ga4); gw[3] = gw4/36.0;
  } else if (gn == 5){
    ga[0] = 0.0;                gw[0] = 0.5688888888888889;
    ga[1] = 0.5384693101056831; gw[1] = 0.4786286704993665;
    ga[2] = 0.9061798459386640; gw[2] = 0.2369268850561891;
    ga[3] = -ga[1];             gw[3] = gw[1];
    ga[4] = -gw[2];             gw[4] = gw[2];
  } else if (gn == 6){
    ga[0] = 0.2386191860831969; gw[0] = 0.4679139345726910;
    ga[1] = 0.6612093864662645; gw[1] = 0.3607615730481386;
    ga[2] = 0.9324695142031521; gw[2] = 0.1713244923791704;
    ga[3] = -ga[0];             gw[3] = gw[0];
    ga[4] = -ga[1];             gw[4] = gw[1];
    ga[5] = -ga[2];             gw[5] = gw[2];
  } else if (gn == 8){
    ga[0] = 0.1834346424956498; gw[0] = 0.3626837833783620;
    ga[1] = 0.5255324099163290; gw[1] = 0.3137066458778873;
    ga[2] = 0.7966664774136267; gw[2] = 0.2223810344533745;
    ga[3] = 0.9602898564975363; gw[3] = 0.1012285362903763;
    ga[4] = -ga[0];             gw[4] = gw[0];
    ga[5] = -ga[1];             gw[5] = gw[1];
    ga[6] = -ga[2];             gw[6] = gw[2];
    ga[7] = -ga[3];             gw[7] = gw[3];
  }else{
    std::cout << "! GaussianQuad() not defined for gn > 8 " << std::endl;
    exit(1);
  }


  sum = 0.0;                // initialize total sum
  for (i=0; i < N; i++) {
    mid = a + (i + 0.5)*dx; // midpoint of the interval
    sum_sub = 0.0;          // initialize sum k-subinterval
    for (k=0; k < gn; k++){
      sum_sub +=  gw[k] *func(0.5*dx * ga[k] + mid); // gaussian to subintervals
    }
    sum_sub *= 0.5*dx;
    sum += sum_sub;         // here the total sum is computed
  }

  return sum;
}
