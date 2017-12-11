/**
* @file solve.h
* @Brief the host  
* @author Dong Yijing,1701110028@pku.edu.cn
* @version 1
* @date 2017-12-04
*/

#include <stdio.h>
#include <fftw3.h>
#include <complex.h>
#include <math.h>
#define N 16
#define lambda 0.5
#define iter 1e6
#define PI M_PI

/**
* @Brief calculate the value of the right hand side term  
*
* @Param x the value of point along x
* @Param y the value of point along y
* @Param z the value of point along z
*
* @Returns sourceterm value of a point  
*/
double complex sourceterm(double x, double y, double z, int flag_p);


/**
* @Brief Do 3D FFT by using 1D FFT
*
* @Param input the value need to do FFT of iFFT
* @Param output the output value
* @Param n the number of value input 
* @Param sign -1-FFT, 1-iFFT
*
* @Returns the output value 
*/
int FFT_3D_complex (fftw_complex ***input, fftw_complex ***output, int n, int sign);

/**
* @Brief 
*
* @Param i the coefficiant of the number
* @Param N The number of total points
*
* @Returns i or N - i
*/
int get_mode(int i, int n);

/**Do a iteration*/
int Relaxation_Iteration(fftw_complex ***phi_1_hat, fftw_complex ***phi_0_hat,       
  fftw_complex ***g_hat, int n, double lamb, int myid);


/**
* @Brief   send the value of B to A
*
* @Param A the matrix to Recv value
* @Param B the matrix to send value
* @Param s the coefficient need to add to suit the FFT
* @Param n the number of input in each direction
*/
void assign_matrix(fftw_complex ***A, fftw_complex ***B, int s, int n);

/**
* @Brief Calculate maximum error between two calculation  
*
* @Param phi_1 the new value of iteration
* @Param phi_0 the old value of iteration
* @Param n the total number of inputfor each direction
* 
* @Return err the maximum err
*/
void error(fftw_complex ***phi_1, fftw_complex ***phi_0, int n);

/**
* @Brief Initial the starting coordinate of each process in true coordinate
*
* @Param lab vector to store the coordinate 
* @Param myid the rank of each process
*/
void Ini_label(int *lab, int myid);

/**
* @Brief return the coresponding coordinate of each proccess
*
* @Param xlabel x-coordinate
* @Param ylabel y-coordinate
* @Param zlabel z-coordinate
* @Param rank myid of each process
* @Param n total points of each direction
*/
void label_corr(double *xlabel, double *ylabel, double *zlabel, int rank, int n);

/** 
* @Brief define corresponding rule in each process and each direction   
*
* @Param myid the rank of your process
* @Param direction x or y or z
*
* @Returns the corresponding process  
*/
int corr_procs(int myid, int direction);


/**
* @Brief Combine the FFT of even and odd numbers together  
*
* @Param combine the combined result,this result should be equal to 3D_FFT result
* @Param data the data stored in this process
* @Param buff the data send to this process
* @Param basic the basic of this FFT
* @Param myid the rank of this process
* @Param dim the direction of this combination
* @Param j the coordinate of one of another two direction
* @Param k the coordinate of one of another two direction
*/
void Combine_FFT (fftw_complex ***combine, fftw_complex ***data, fftw_complex *buff,
  fftw_complex *basic, int myid, int dim, int j, int k, int tag);
