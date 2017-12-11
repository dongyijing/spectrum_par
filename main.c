#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <math.h>
#include <fftw3.h>
#include <time.h>
#include <unistd.h>
#include <float.h>
#include <string.h>
#include "solve.h"
#include "mpi.h"

int main(int argc, char *argv[]){

  int i, j, k;
  int iter_count = 1;
  int dim;
  int corr_rank;
  int lab[3];
  double tol = 1.0e-8;
  double h;
  double m;
  double err;/*error in each process*/
  double totalerr = 1.0;/**the maxerror*/
  FILE *answer, *residual;
  complex source;

  double *xlabel;
  double *ylabel;
  double *zlabel;

  int pos = 1;
  int flag_p = 0;

  int size = 8;
  int rank;

  double t0, t1;

  fftw_complex *send_buff;
  fftw_complex *recv_buff;

  fftw_complex *basic;
  fftw_complex *ibasic;

  /**< data in position place*/
  fftw_complex ***phi_1;
  fftw_complex ***phi_0;
  fftw_complex ***g;
  fftw_complex ***comb_result_p;/**store the value combined*/
  fftw_complex ***comb_result_g;

  /**<data in reciporal space*/
  fftw_complex ***phi_1_hat;
  fftw_complex ***phi_0_hat;
  fftw_complex ***g_hat;

  /**MPI initialize*/
  MPI_Status status;
  MPI_Datatype fftwcplx;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /*while (pos < argc){
    if (!strcmp(argv[pos], "-problem")){
	  flag_p = (int)atoi(argv[++pos]);
	  pos++;
	  if (flag_p == 0)
	    fprintf(stderr, "The true problem\n");
	  else if(flag_p == 1)
	    fprintf(stderr, "The check problem\n");
	  else{
	    fprintf(stderr, "hello, little pig\n");
		flag_p = 0;
		return -1;
	  }
	}
	else
	  break;
  }*/
/**  
  switch(argc){

    case 0:
	  break;
	
	default:
	  fprintf(stderr, "Hello, clever boy! Here's my code! You can run it now!\n");
	  return -1;
  }
*/
  h = 2*PI/N;

  
  /**We need 8 process to parallel, so we divided points into 8 parts,for each process,
     we just store (N/2)^(3) points, so we define the matrix as below*/
  phi_1         = fftw_malloc(sizeof(fftw_complex) * N/2);
  phi_0         = fftw_malloc(sizeof(fftw_complex) * N/2);
  g             = fftw_malloc(sizeof(fftw_complex) * N/2);
  phi_0_hat     = fftw_malloc(sizeof(fftw_complex) * N/2);
  phi_1_hat     = fftw_malloc(sizeof(fftw_complex) * N/2);
  g_hat         = fftw_malloc(sizeof(fftw_complex) * N/2);
  comb_result_p = fftw_malloc(sizeof(fftw_complex) * N/2);
  comb_result_g = fftw_malloc(sizeof(fftw_complex) * N/2);

  for (i = 0; i < N/2; i++){
    phi_0[i]          = fftw_malloc(sizeof(fftw_complex) * N/2);
    phi_1[i]          = fftw_malloc(sizeof(fftw_complex) * N/2);
    g[i]              = fftw_malloc(sizeof(fftw_complex) * N/2);
    phi_0_hat[i]      = fftw_malloc(sizeof(fftw_complex) * N/2);
    phi_1_hat[i]      = fftw_malloc(sizeof(fftw_complex) * N/2);
    g_hat[i]          = fftw_malloc(sizeof(fftw_complex) * N/2);
    comb_result_g [i] = fftw_malloc(sizeof(fftw_complex) * N/2);
    comb_result_p [i] = fftw_malloc(sizeof(fftw_complex) * N/2);
      
	  for (j = 0; j < N/2; j++){
	    phi_0[i][j]         = fftw_malloc(sizeof(fftw_complex) * N/2);
	    phi_1[i][j]         = fftw_malloc(sizeof(fftw_complex) * N/2);
	    g[i][j]             = fftw_malloc(sizeof(fftw_complex) * N/2);
	    phi_0_hat[i][j]     = fftw_malloc(sizeof(fftw_complex) * N/2);
	    phi_1_hat[i][j]     = fftw_malloc(sizeof(fftw_complex) * N/2);
	    g_hat[i][j]         = fftw_malloc(sizeof(fftw_complex) * N/2);
	    comb_result_p[i][j] = fftw_malloc(sizeof(fftw_complex) * N/2);
	    comb_result_g[i][j] = fftw_malloc(sizeof(fftw_complex) * N/2);
	  }
  }
  /******************************************/
  /**< Initialize the value in every process*/
  /******************************************/
  for (i = 0; i < N/2; i++){
    for (j = 0; j < N/2; j++){
	  for (k = 0; k < N/2; k++){
	    phi_0[i][j][k] = 0.25;
	  }
	}
  }
  /*******************************************************/ 
  /** Correspond the label in every process to true label*/
  /*******************************************************/
  xlabel = malloc(sizeof(double) * N/2);
  ylabel = malloc(sizeof(double) * N/2);
  zlabel = malloc(sizeof(double) * N/2);
  label_corr(xlabel, ylabel, zlabel, rank, N/2);
  

  /**prepare buffer to store the value need to communicate*/
  send_buff = malloc(sizeof(fftw_complex) * N/2);
  recv_buff= malloc(sizeof(fftw_complex) * N/2);

  MPI_Type_contiguous(2, MPI_DOUBLE, &fftwcplx);
  MPI_Type_commit(&fftwcplx);

  /**************************/
  /**Define the basic of FFT*/
  /**************************/
  basic  = malloc(sizeof(fftw_complex) * N/2);
  ibasic = malloc(sizeof(fftw_complex) * N/2);
  for (i = 0; i < N/2; i++){
     basic[i] = cos(2*PI*i/N) - I*sin(2*PI*i/N);
    ibasic[i] = cos(2*PI*i/N) + I*sin(2*PI*i/N);
  }

  /* residual = fopen("../spectrum_sol/residual.txt","w");*/

  /** Start to calculate time*/
  t0 = MPI_Wtime();

  /***********************/
  /** Start the iteration*/
  /***********************/
  while (totalerr > tol){
    if (iter_count > iter){
	  break;
    }
    /**calculate g */
    for (i = 0; i < N/2; i++){
      for (j = 0; j < N/2; j++){
	    for (k = 0; k < N/2; k++){
	      source =   sourceterm(xlabel[i], ylabel[j], zlabel[k], flag_p);
		  g[i][j][k] =- conj(phi_0[i][j][k])*phi_0[i][j][k]*phi_0[i][j][k] 
		               + source;
        }
	  }
    }

    /** Do FFT in each process*/
    /*===========================================*/
    /*===========================================*/  	
	/*
    if (rank == 0){
	  for (i = 0; i < N/2; i++){
	    for (j = 0; j < N/2; j++){
	      for (k = 0; k < N/2; k++){
		    fprintf(stderr, "%d :%f + %fi\n",rank, creal(g[i][j][k]), cimag(g[i][j][k]));
		  }
        }
	  }
	}
	*/
    
    /******************************/
    /**Combine FFT of each process*/
    /******************************/

	/*==================================================================================*/
    /**Combine it of x-direction*/
	
    /** Run 2^(s-1) fast Fourier transform in each processor*/
    FFT_3D_complex(phi_0, phi_0_hat, N/2, -1);
    FFT_3D_complex(g, g_hat, N/2, -1);
	/*if (rank == 0){
	  for (i = 0; i < N/2; i++){
	    for (j = 0; j < N/2; j++){
		  for (k = 0; k < N/2; k++){
		    printf("g[%d][%d][%d] = %f + %fi\n",i,j,k, creal(g[i][j][k]), cimag(g[i][j][k]));
		  }
		}
      }
	}*/

   /** Data exchange and merge the 2^(s-1) FFTs into 2^s FFTs in z-direction*/

   Ini_label(lab,rank);
   /** The first step: merge FFTs in x-direction*/
   /** u0_hat */
   dim = 0;
   for (i = 0; i < N/2; i++){
     for (j = 0; j < N/2; j++){

       /** Data exchange */
       for (k = 0; k < N/2; k++){
         send_buff[k] = phi_0_hat[k][i][j];
       }

       corr_rank = corr_procs(rank, dim);
 
       MPI_Send(send_buff, N/2, fftwcplx, corr_rank, 99, MPI_COMM_WORLD);
       MPI_Recv(recv_buff, N/2, fftwcplx, corr_rank, 99, MPI_COMM_WORLD, &status);
	 /* if (rank == 0){
	       printf("send:send_buff[0] = %f + %fi, send_buff[1] = %f + %fi\n",creal(send_buff[0]),cimag(send_buff[0]),creal(send_buff[1]),cimag(send_buff[1]));
	       printf("recv:recv_buff[0] = %f + %fi, recv_buff[1] = %f + %fi\n",creal(recv_buff[0]),cimag(recv_buff[0]),creal(send_buff[1]),cimag(send_buff[1]));
	  }
	   if (rank == 4){
	       printf("send:send_buff[0] = %f + %fi, send_buff[1] = %f + %fi\n",creal(send_buff[0]),cimag(send_buff[0]),creal(send_buff[1]),cimag(send_buff[1]));
	       printf("recv:recv_buff[0] = %f + %fi, recv_buff[1] = %f + %fi\n",creal(recv_buff[0]),cimag(recv_buff[0]),creal(send_buff[1]),cimag(send_buff[1]));
	   }
	   sleep(60);*/
       
       /** Merge the FFTs */
       Combine_FFT(comb_result_p, phi_0_hat, recv_buff, basic, rank, dim, i, j, lab[dim]);
     }
   }
   assign_matrix(phi_0_hat, comb_result_p, 2.0, N/2);
 /* if (rank == 0){
     for (i = 0; i < N/2; i++){
	   for(j = 0; j < N/2; j++){
	     for (k = 0; k < N/2; k++){
		   printf("g_hat[%d][%d][%d] = %f + %fi\n", 2*i+lab[0], 2*j+lab[1], 2*k, creal(phi_0_hat[i][j][k]), cimag(phi_0_hat[i][j][k]));
		 }
	   }
	 }
   }
   sleep(60);
 */
   /** g_hat */
   dim = 0;
   for (i = 0; i < N/2; i++){
     for (j = 0; j < N/2; j++){

       /** Data exchange */
       for (k = 0; k < N/2; k++){
         send_buff[k] = g_hat[k][i][j];
       }

       corr_rank = corr_procs(rank, dim);

       MPI_Send(send_buff, N/2, fftwcplx, corr_rank, 100, MPI_COMM_WORLD);
       MPI_Recv(recv_buff, N/2, fftwcplx, corr_rank, 100, MPI_COMM_WORLD, &status);

       /** Merge the FFTs */
       Combine_FFT(comb_result_g, g_hat, recv_buff, basic, rank, dim, i, j, lab[dim]);
     }
   }
   assign_matrix(g_hat, comb_result_g, 2.0, N/2);
  /*if (rank == 0){
     for (i = 0; i < N/2; i++){
	   for(j = 0; j < N/2; j++){
	     for (k = 0; k < N/2; k++){
		   printf("g_hat[%d][%d][%d] = %f + %fi\n", 2*i+lab[0], 2*j+lab[1], 2*k, creal(phi_0_hat[i][j][k]), cimag(phi_0_hat[i][j][k]));
		 }
	   }
	 }
   }
   sleep(60);*/
   /** The second step: merge FFTs in y-direction*/
   dim = 1;
   for (i = 0; i < N/2; i++){
     for (j = 0; j < N/2; j++){

       /** Data exchange */
       for (k = 0; k < N/2; k++){
         send_buff[k] = phi_0_hat[i][k][j];
       }

       corr_rank = corr_procs(rank, dim);
 
       MPI_Send(send_buff, N/2, fftwcplx, corr_rank, 101, MPI_COMM_WORLD);
       MPI_Recv(recv_buff, N/2, fftwcplx, corr_rank, 101, MPI_COMM_WORLD, &status);
	  /* if (rank == 5){
	       printf("send:send_buff[0] = %f + %fi, send_buff[1] = %f + %fi\n",creal(send_buff[0]),cimag(send_buff[0]),creal(send_buff[1]),cimag(send_buff[1]));
	       printf("recv:recv_buff[0] = %f + %fi, recv_buff[1] = %f + %fi\n",creal(recv_buff[0]),cimag(recv_buff[0]),creal(send_buff[1]),cimag(send_buff[1]));
	   }
	   if (rank == 7){
	       printf("send:send_buff[0] = %f + %fi, send_buff[1] = %f + %fi\n",creal(send_buff[0]),cimag(send_buff[0]),creal(send_buff[1]),cimag(send_buff[1]));
	       printf("recv:recv_buff[0] = %f + %fi, recv_buff[1] = %f + %fi\n",creal(recv_buff[0]),cimag(recv_buff[0]),creal(send_buff[1]),cimag(send_buff[1]));
	   }
	   sleep(60);*/
       Combine_FFT(comb_result_p, phi_0_hat, recv_buff, basic, rank, dim, i, j,lab[dim]);
     }
   }
   assign_matrix(phi_0_hat, comb_result_p, 2.0, N/2);
  /* if (rank == 0){
     for (i = 0; i < N/2; i++){
	   for(j = 0; j < N/2; j++){
	     for (k = 0; k < N/2; k++){
		   printf("g_hat[%d][%d][%d] = %f + %fi\n", 2*i+lab[0], 2*j+lab[1], 2*k, creal(phi_0_hat[i][j][k]), cimag(phi_0_hat[i][j][k]));
		 }
	   }
	 }
   }
   sleep(60);*/
   dim = 1;
   for (i = 0; i < N/2; i++){
     for (j = 0; j < N/2; j++){

       /** Data exchange */
       for (k = 0; k < N/2; k++){
         send_buff[k] = g_hat[i][k][j];
       }

       corr_rank = corr_procs(rank, dim);

       MPI_Send(send_buff, N/2, fftwcplx, corr_rank, 102, MPI_COMM_WORLD);
       MPI_Recv(recv_buff, N/2, fftwcplx, corr_rank, 102, MPI_COMM_WORLD, &status);
	  /* if (rank == 0){
	       printf("send:send_buff[0] = %f + %fi, send_buff[1] = %f + %fi\n",creal(send_buff[0]),cimag(send_buff[0]),creal(send_buff[1]),cimag(send_buff[1]));
	       printf("recv:recv_buff[0] = %f + %fi, recv_buff[1] = %f + %fi\n",creal(recv_buff[0]),cimag(recv_buff[0]),creal(send_buff[1]),cimag(send_buff[1]));
	   }
	   if (rank == 1){
	       printf("send:send_buff[0] = %f + %fi, send_buff[1] = %f + %fi\n",creal(send_buff[0]),cimag(send_buff[0]),creal(send_buff[1]),cimag(send_buff[1]));
	       printf("recv:recv_buff[0] = %f + %fi, recv_buff[1] = %f + %fi\n",creal(recv_buff[0]),cimag(recv_buff[0]),creal(send_buff[1]),cimag(send_buff[1]));
	   }
	   sleep(60);*/

       /** Merge the FFTs */
       Combine_FFT(comb_result_g, g_hat, recv_buff, basic, rank, dim, i, j,lab[dim]);
     }
   }
   assign_matrix(g_hat, comb_result_g, 2.0, N/2);
  /* if (rank == 0){
     for (i = 0; i < N/2; i++){
	   for(j = 0; j < N/2; j++){
	     for (k = 0; k < N/2; k++){
		   printf("g_hat[%d][%d][%d] = %f + %fi\n", 2*i+lab[0], 2*j+lab[1], 2*k, creal(phi_0_hat[i][j][k]), cimag(phi_0_hat[i][j][k]));
		 }
	   }
	 }
   }
   sleep(60);*/
   dim = 2;
   for (i = 0; i < N/2; i++){
     for (j = 0; j < N/2; j++){

       /** Data exchange */
       for (k = 0; k < N/2; k++){
         send_buff[k] = phi_0_hat[i][j][k];
       }

       corr_rank = corr_procs(rank, dim);
 
       MPI_Send(send_buff, N/2, fftwcplx, corr_rank, 103, MPI_COMM_WORLD);
       MPI_Recv(recv_buff, N/2, fftwcplx, corr_rank, 103, MPI_COMM_WORLD, &status);
	   /*if (rank == 0){
	       printf("send:send_buff[0] = %f + %fi, send_buff[1] = %f + %fi\n",creal(send_buff[0]),cimag(send_buff[0]),creal(send_buff[1]),cimag(send_buff[1]));
	       printf("recv:recv_buff[0] = %f + %fi, recv_buff[1] = %f + %fi\n",creal(recv_buff[0]),cimag(recv_buff[0]),creal(send_buff[1]),cimag(send_buff[1]));
	   }
	   if (rank == 1){
	       printf("send:send_buff[0] = %f + %fi, send_buff[1] = %f + %fi\n",creal(send_buff[0]),cimag(send_buff[0]),creal(send_buff[1]),cimag(send_buff[1]));
	       printf("recv:recv_buff[0] = %f + %fi, recv_buff[1] = %f + %fi\n",creal(recv_buff[0]),cimag(recv_buff[0]),creal(send_buff[1]),cimag(send_buff[1]));
	   }
	   sleep(60);*/

       /** Merge the FFTs */
       Combine_FFT(comb_result_p, phi_0_hat, recv_buff, basic, rank, dim, i, j,lab[dim]);
     }
   }
   assign_matrix(phi_0_hat, comb_result_p, 2.0, N/2);
   dim = 2;
   for (i = 0; i < N/2; i++){
     for (j = 0; j < N/2; j++){
       for (k = 0; k < N/2; k++){
         send_buff[k] = g_hat[i][j][k];
       }

       corr_rank = corr_procs(rank, dim);

       MPI_Send(send_buff, N/2, fftwcplx, corr_rank, 104, MPI_COMM_WORLD);
       MPI_Recv(recv_buff, N/2, fftwcplx, corr_rank, 104, MPI_COMM_WORLD, &status);
       Combine_FFT(comb_result_g, g_hat, recv_buff,  basic, rank, dim, i, j, lab[dim]);
     }
   }
   assign_matrix(g_hat, comb_result_g, 2.0, N/2);
  /* if (rank == 0){
	 printf("g_hat[0][0][0] = %f + %fi\n",  creal(phi_0_hat[0][0][0]), cimag(phi_0_hat[0][0][0]));
	 printf("g_hat[2][6][8] = %f + %fi\n",  creal(phi_0_hat[1][3][4]), cimag(phi_0_hat[1][3][4]));
   }
   sleep(5);*/
   /*if (rank == 7){
     for (i = 0; i < N/2; i++){
	   for(j = 0; j < N/2; j++){
	     for (k = 0; k < N/2; k++){
		   printf("g_hat[%d][%d][%d] = %f + %fi\n", 2*i+1, 2*j+1, 2*k+1, creal(g_hat[i][j][k]), cimag(g_hat[i][j][k]));
		 }
	   }
	 }
   }
   sleep(30);*/
    /*************************************************/
    /**iterate phi_0_hat to phi_1_hat in each process*/
    /*************************************************/
   Relaxation_Iteration(phi_1_hat, phi_0_hat, g_hat, N/2, lambda, rank);
   /*for(i = 0; i < N/2; i++){
     for (j = 0; j < N/2; j++){
	   for (k = 0; k < N/2; k++){
	     phi_1_hat[i][j][k] = phi_0_hat[i][j][k];
       }
	 }
   }*/

    /***************************************************************/
    /**transfer the output value phi_1_hat to phi_1 in each process*/
    /***************************************************************/
   FFT_3D_complex (phi_1_hat, phi_1, N/2, 1);
   /*if(rank == 1){
   for(i = 0; i < N/2; i++){
     for (j = 0; j < N/2; j++){
	   for (k = 0; k < N/2; k++){
	    printf("%f + %fi\n",creal(phi_1[i][j][k]), cimag(phi_1[i][j][k]));
       }
	 }
   }
   }
   sleep(60);*/

   /**************************************************/
   /**Do combination to combine phi_1 in each process*/
   /**************************************************/

   /*==================================================================================*/
   /**Combine it of x-direction*/
   dim = 0;
   for (j = 0; j < N/2; j++){
     for (k = 0; k < N/2; k++){
	   /**Copy the values need to communicate with other process to send_buff*/
	   for (i = 0; i < N/2; i++){
	     send_buff[i] = phi_1[i][j][k];
	   }
	   corr_rank = corr_procs(rank, dim);

       MPI_Send(send_buff, N/2, fftwcplx, corr_rank, 105, MPI_COMM_WORLD);
	   MPI_Recv(recv_buff, N/2, fftwcplx, corr_rank, 105, MPI_COMM_WORLD, &status);
	  /* 
	   if (rank == 0){
	       printf("send:send_buff[0] = %f + %fi, send_buff[1] = %f + %fi\n",creal(send_buff[0]),cimag(send_buff[0]),creal(send_buff[1]),cimag(send_buff[1]));
	       printf("recv:recv_buff[0] = %f + %fi, recv_buff[1] = %f + %fi\n",creal(recv_buff[0]),cimag(recv_buff[0]),creal(send_buff[1]),cimag(send_buff[1]));
	   }
	   if (rank == 4){
	       printf("send:send_buff[0] = %f + %fi, send_buff[1] = %f + %fi\n",creal(send_buff[0]),cimag(send_buff[0]),creal(send_buff[1]),cimag(send_buff[1]));
	       printf("recv:recv_buff[0] = %f + %fi, recv_buff[1] = %f + %fi\n",creal(recv_buff[0]),cimag(recv_buff[0]),creal(send_buff[1]),cimag(send_buff[1]));
	   }
	   sleep(60);*/
	   Combine_FFT(comb_result_p, phi_1, recv_buff, ibasic, rank, dim, j, k, lab[dim]);
	 }
   }
   assign_matrix(phi_1, comb_result_p, 1, N/2);
   /*==================================================================================*/

   /*==================================================================================*/
   /**Combine it of y-direction*/
   dim = 1;
   for (j = 0; j < N/2; j++){
     for (k = 0; k < N/2; k++){
	   /**Copy the values need to communicate with other process to send_buff*/
	   for (i = 0; i < N/2; i++){
	     send_buff[i] = phi_1[j][i][k];
	   }
	   corr_rank = corr_procs(rank, dim);

       MPI_Send(send_buff, N/2, fftwcplx, corr_rank, 106, MPI_COMM_WORLD);
	   MPI_Recv(recv_buff, N/2, fftwcplx, corr_rank, 106, MPI_COMM_WORLD, &status);
	   Combine_FFT(comb_result_p, phi_1, recv_buff, ibasic, rank, dim, j, k, lab[dim]);
	 }
   }
   assign_matrix(phi_1, comb_result_p, 1, N/2);
   /*==================================================================================*/

   /*==================================================================================*/
   /**Combine it of z-direction*/
   dim = 2;
   for (j = 0; j < N/2; j++){
     for (k = 0; k < N/2; k++){
	   /**Copy the values need to communicate with other process to send_buff*/
	   for (i = 0; i < N/2; i++){
	     send_buff[i] = phi_1[j][k][i];
	   }
	   corr_rank = corr_procs(rank, dim);

       MPI_Send(send_buff, N/2, fftwcplx, corr_rank, 107, MPI_COMM_WORLD);
	   MPI_Recv(recv_buff, N/2, fftwcplx, corr_rank, 107, MPI_COMM_WORLD, &status);
	   Combine_FFT(comb_result_p, phi_1, recv_buff, ibasic, rank, dim, j, k,lab[dim]);
	 }
   }
   assign_matrix(phi_1, comb_result_p, 1, N/2);
   /*==================================================================================*/
   /*if (rank == 6){
     for (i = 0; i < N/2; i++){
	   for (j = 0; j < N/2; j++){
	     for (k = 0; k < N/2; k++){
		   printf("%f+%fi\n", creal(phi_1[i][j][k]), cimag(phi_1[i][j][k]));
		 }
	   }
	 }
   }*/

   
   /**********************************/
   /**Calculate error of each process*/
   /**********************************/
   err = 0;
   for (i = 0;i < N/2; i++){
     for (j = 0; j < N/2; j++){
	   for (k = 0; k < N/2; k++){
	     m = cabs(phi_1[i][j][k] - phi_0[i][j][k]);
		 if (m > err){
		   err = m;
		 }
	   }
	 }
   }
   MPI_Reduce(&err, &totalerr, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
   MPI_Bcast(&totalerr, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
     /*fprintf(residual,"%f\n",totalerr);*/
   if (rank == 0){
     printf("residual  %f\n", totalerr);
   }
   assign_matrix(phi_0, phi_1, 1, N/2);
   iter_count++;
 }
 /*fclose(residual);*/

 MPI_Finalize();
/*if (flag_p = 1){
   error = 0;
   for (i = 0; i < N; i++){
     for (j = 0; j < N; j++){
	   for (k = 0; k < N; k++){
	       x = -PI + i*h;
	       y = -PI + j*h;
               z = -PI + k*h;
	       m = cabs(phi_0[i][j][k] - (cos(x)*cos(y)*cos(z) + sin(x)*sin(y)*sin(z)*I));
		 if (m > error){
		   error = m;
		 }
	   }
	 }
   }
 }
 printf("error = %f\n",error);
 */

 answer = fopen("../spectrum_sol/par_answer.txt","w");
   for (i = 0; i < N/2; i++){
     for (j = 0; j < N/2; j++){
       for (k = 0; k < N/2; k++){
	     fprintf(answer, "%f + %f*i\n",creal(phi_0[i][j][k]), cimag(phi_0[i][j][k]));
	   }
     }
   }
 
 fclose(answer);
 
 if (rank == 0){
   printf("phi_0[0][0][0] = %f + %fi\n",creal(phi_0[0][0][0]), cimag(phi_0[0][0][0]));
 }

 fftw_free(phi_0);
 fftw_free(phi_1);
 fftw_free(g);
 fftw_free(phi_0_hat);
 fftw_free(phi_1_hat);
 fftw_free(g_hat);
 /*fftw_free(comb_result_p);
 fftw_free(comb_result_g);*/
 if (rank == 0){
   printf("I have printed %d times\n", iter_count);
 }
}

   
 



