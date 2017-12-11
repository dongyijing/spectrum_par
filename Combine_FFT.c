
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "solve.h"

void Combine_FFT (fftw_complex ***combine, fftw_complex ***data, fftw_complex *buff, 
                  fftw_complex *basic, int myid, int dim, int j, int k, int tag){
  /**Do combination in x-direction*/
  int i;
  int lab[3];
  Ini_label(lab, myid);
    if (tag == 0){/**even number*/
	  for(i = 0; i < N/4; i++){  
	  if(dim == 0){
	      combine[i][j][k]       = data[2*i][j][k] + basic[2*i] * (buff[2*i]);
	      combine[i + N/4][j][k] = data[2*i][j][k] - basic[2*i] * (buff[2*i]);
	  }
	  if(dim == 1){
	      combine[j][i][k]       = data[j][2*i][k] + basic[2*i] * buff[2*i];
	      combine[j][i + N/4][k] = data[j][2*i][k] - basic[2*i] * buff[2*i];
	  }
	  if(dim == 2){
	      combine[j][k][i]       = data[j][k][2*i] + basic[2*i] * buff[2*i];
	      combine[j][k][i + N/4] = data[j][k][2*i] - basic[2*i] * buff[2*i];
	  }
	  }
	}
    if (tag == 1){/**odd number*/
	for (i = 0; i < N/4; i++){
	  if(dim == 0){
	      combine[i][j][k]       = buff[2*i + 1] + basic[2*i + 1] * data[2*i + 1][j][k];
	      combine[i + N/4][j][k] = buff[2*i + 1] - basic[2*i + 1] * data[2*i + 1][j][k];
	  }
	  if(dim == 1){
	      combine[j][i][k]       = buff[2*i + 1] + basic[2*i + 1] * data[j][2*i + 1][k];
	      combine[j][i + N/4][k] = buff[2*i + 1] - basic[2*i + 1] * data[j][2*i + 1][k];
	  }
	  if(dim == 2){
	      combine[j][k][i]       = buff[2*i + 1] + basic[2*i + 1] * data[j][k][2*i + 1];
	      combine[j][k][i + N/4] = buff[2*i + 1] - basic[2*i + 1] * data[j][k][2*i + 1];
	  }
	}
	}
}




	  
    
	
