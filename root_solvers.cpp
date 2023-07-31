#include "prototypes.h"


/////////////////////////////////////////////////Bisection method
double bisection(double (*func)(double), double a, double b, double xtol, double ftol, double &zero, int &attempts)
{
  int maxattempt = 128; //setting iteration max value
  double fa = func(a), fb = func(b);
  double mid, fmid;
  double del = b-a; // distance between endpoints
  double dx = fabs(del);

  if (fa * fb < 0){
    for (int k=1; k<maxattempt; k++){

      mid = (a+b)*0.5; // find midpoint
      fmid = func(mid);
      del = b - a;
      dx = fabs(del);

      cout << "bisection():  " << " # = " << k
           << "      [a,b] =  [ " << a << ", " << b << " ]; "
           << "      mid = " << mid
           << "      fmid = " << fmid
           << "      |dx| = " << dx<< endl;

      if (fabs(del) < xtol || fabs(fmid) < ftol){
        attempts = k;
        zero = mid;
        return 0;
      }

      //updating interval boundaries
      if (fmid*fa < 0.0){
        b = mid;
      }else{
        a = mid;
      }
    }
    cout << "! bisection(): max attempts reached" << endl;
    return 1;
  }else{
    cout << "! bisection(): no root in the initial interval" << endl;
    return -1;
  }
}

/////////////////////////////////////////////////Falseposition method
double falseposition(double (*func)(double), double a, double b, double xtol, double ftol, double &zero, int &attempts)
{
  int maxattempt = 128; //setting iteration max value
  double fa = func(a), fb = func(b);
  double mid, fmid;
  double del = b-a; // distance between endpoints
  double dx = fabs(del);

  if (fa * fb < 0){
    for (int k=1; k<maxattempt; k++){

      mid = (a*fb - b*fa)/(fb - fa); // find midpoint
      fmid = func(mid);
      dx = fabs(del);

      cout << "falseposition():  " << " # = " << k
           << "      [a,b] =  [ " << a << ", " << b << " ]; "
           << "      mid = " << mid
           << "      fmid = " << fmid
           << "      |dx| = " << fabs(del)<< endl;

      //updating interval boundarys
      if (fmid*fa < 0.0){
        del = mid -b;
        b = mid;
        fb = fmid;
      }else{
        del = mid -a;
        a = mid;
        fa = fmid;
      }

      if (dx < xtol || fabs(fmid) < ftol){
        attempts = k;
        zero = mid;
        return 0;
      }
    }
    cout << "! falseposition(): max attempts reached" << endl;
    return 1;
  }else{
    cout << "! falseposition(): no root in the initial interval" << endl;
    return -1;
  }
}

/////////////////////////////////////////////////Secant method
double secant(double (*func)(double), double a, double b, double xtol, double ftol, int &attempts)
{
  int maxattempt = 63;    //setting iteration max value
  double fa = func(a), fb = func(b);
  double delx = b-a;      // distance between endpoints

  if (fa * fb < 0){
    for (int k=1; k<=maxattempt; k++){

      delx = fb*(b - a)/(fb - fa); // computing increment
      a = b;
      fa = fb;
      b = b - delx;
      fb = func(b);

      cout << "secant():  " << " # = " << k
           << "      [a,b] =  [ " << a << ", " << b << " ]; "
           << "      dx = " << delx << endl;

      if (fabs(delx) < xtol || fabs(b) < ftol){
        attempts = k;
        return b;
      }
    }
    cout << "! secant(): max attempts reached" << endl;
    return 1;
  }else{
    cout << "! secant(): no root in the initial interval" << endl;
    return -1;
  }
}

/////////////////////////////////////////////////Newton method
double newton(double (*func)(double), double (*dfunc)(double), double a, double b, double xtol, double ftol, int &attempts)
{
  int maxattempt = 63; //setting iteration max value
  double fa = func(a), fb = func(b);
  double dx, fmid;

  double mid = (a+b)*0.5;
  if (fa * fb < 0){
    int k=0;
    dx = 1.0;
    while ( fabs(dx) > xtol){
      k++;
      if (k == maxattempt){
        cout << "! newton(): max attempts reached" << endl;
        return 1;
      }
      fmid = func(mid);
      dx = fmid/dfunc(mid);
      mid -= dx;

      cout << "newton():  " << " # = " << k
           << "  mid = " << mid << endl;
    }
    attempts = k;
    return mid;
  }else{
    cout << "! newton(): no root in the initial interval" << endl;
    return -1;
  }
}

/////////////////////////////////////////////////Bracket method
void bracket(double (*func)(double), double a, double b, int n, double *xL, double *xR, int &nroots)
{
  int count = 0, i;
  double x,fL,fR,dx;

  dx = (b-a)/n; // determining subintervals number
  fL =  func(x=a);
  for (i = 0; i < n; i++){
    fR = func(x += dx);

    if(fL*fR <= 0.0){ // check a sign change
      xL[count] = x - dx;
      xR[count] = x;
      count++;        // count sign change
      if(count == nroots) {
        return;
      }
      fL = fR;
    }
    nroots = count;
  }


}
