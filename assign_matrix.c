#include <stdio.h>
#include <complex.h>
#include <fftw3.h>
#include "solve.h"

void assign_matrix(fftw_complex ***A, fftw_complex ***B, int s, int n){

  int i, j, k;
  for (i = 0; i < n; i++){
    for (j = 0; j < n; j++){
	  for (k = 0; k < n; k++){
		  A[i][j][k] = B[i][j][k]/s;
	  }
	} 
  }
}
