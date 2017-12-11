#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "solve.h"

void error(fftw_complex ***phi_1, fftw_complex ***phi_0, int n){
  int i, j, k;
  double error = 0.0;
	for (i = 0; i < n; i++){
	  for (j = 0; j < n; j++){
	    for (k = 0; k < n; k++){
		  /**res = max(res, cabs(phi_1[i][j][k] - phi_0[i][j][k]));*/
		  error = 0.0;
        }
	  }
	}
}
