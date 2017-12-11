#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "solve.h"

int Relaxation_Iteration(fftw_complex ***phi_1_hat, fftw_complex ***phi_0_hat, 
fftw_complex ***g_hat, int n, double lamb, int myid){

  int i, j, k;
  int m; /**save the coefficient of phi_1_hat*/
  int lab[3];
  
  Ini_label(lab, myid);

  for (i = 0; i < N/2; i++){
    for (j = 0; j < N/2; j++){
	  for (k = 0; k < N/2; k++){
	    m = pow(get_mode(i*2 + lab[0], N), 2) + pow(get_mode(j*2 + lab[1], N), 2) +pow(get_mode(k*2 + lab[2], N), 2);
		phi_1_hat[i][j][k] = (lamb*phi_0_hat[i][j][k] + g_hat[i][j][k])/((double)m + lamb);
	  }
	}
  }
}

int get_mode( int i, int n){
  if (i < N/2)
    return i;
  else
    return i-N;
}
