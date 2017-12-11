#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include "solve.h"
/**Aussue we want to transfer data0 to data1*/
int FFT_3D_complex(fftw_complex ***input, fftw_complex ***output, int n, int sign){

  int i, j, k;
  fftw_complex *in, *out;
  //fftw_complex in[N], out[N];
  fftw_plan p;

  /** copy the data */
  for (i = 0; i < n; i++){
    for (j = 0; j < n; j++){
      for (k = 0; k < n; k++){
        output[i][j][k] = input[i][j][k];
      }
    }
  }
  
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n);

  for (i = 0; i < n; i++){
    for (j = 0; j < n; j++){
      /** copy the data into 'in'*/
      for (k = 0; k < n; k++){
        in[k] = output[i][j][k];
      }
        
      /**do FFTs when sign = 1, or do inverse FFTs when sign = -1*/
      p = fftw_plan_dft_1d(n, in, out, sign, FFTW_ESTIMATE);
      fftw_execute(p);
        
      for (k = 0; k < n; k++){
        if(sign == 1){
          output[i][j][k] = out[k];
        }
        else if(sign == -1){
          output[i][j][k] = out[k]/(double)n;
        }
      }
    }
  }
 
  /**The second step: do FFTs in y-direction*/
  for (i = 0; i < n; i++){
    for (k = 0; k < n; k++){
      /** copy the data into 'in'*/
      for (j = 0; j < n; j++){
        in[j] = output[i][j][k];
      }
      /**do FFTs when sign = 1, or do inverse FFTs when sign = -1*/
      p = fftw_plan_dft_1d(n, in, out, sign, FFTW_ESTIMATE);
      fftw_execute(p);
        
      for (j = 0; j < n; j++){
        if(sign == 1){
          output[i][j][k] = out[j];
        }
        else if(sign == -1){
          output[i][j][k] = out[j]/(double)n;
        }
      }
    }
  }

  /**The third step: do FFTs in x-direction*/
  for (j = 0; j < n; j++){
    for (k = 0; k < n; k++){
      /** copy the data into 'in'*/
      for (i = 0; i < n; i++){
        in[i] = output[i][j][k];
      }
        
      /**do FFTs when sign = 1, or do inverse FFTs when sign = -1*/
      p = fftw_plan_dft_1d(n, in, out, sign, FFTW_ESTIMATE);
      fftw_execute(p);
        
      for (i = 0; i < n; i++){
        if(sign == 1){
          output[i][j][k] = out[i];
        }
        else if(sign == -1){
          output[i][j][k] = out[i]/(double)n;
        }
      }
    }
  }

  fftw_destroy_plan(p);
  fftw_free(in);
  fftw_free(out);
  
  return 0;
}
  

  
