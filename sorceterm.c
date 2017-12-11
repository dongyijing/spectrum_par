#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "solve.h"

double complex sourceterm(double x, double y, double z, int flag_p){
 
  double complex fvalue;
  if (flag_p == 0){
    fvalue = sin(x)*sin(y)*sin(z) + cos(x)*cos(y)*cos(z)*I;
  }
  if (flag_p == 1){
    fvalue = ((cos(x)*cos(y)*cos(z)) + (sin(x)*sin(y)*sin(z))*I)*(3*pow((cos(x)*cos(y)*cos(z)),2) + 3*pow(sin(x)*sin(y)*sin(z),2) + 1);
  }
  return fvalue;
}
 
  
