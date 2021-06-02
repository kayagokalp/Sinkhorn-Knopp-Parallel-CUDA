#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <omp.h>
#include <math.h>       /* fabsf */
#include <string.h>
#include <stdlib.h>
#include <time.h>

#define DEBUG 1
#define NO_BLOCKS 256
#define NO_THREADS 32
//Error check-----
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
  if (code != cudaSuccess) 
    {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
    }
}
//Error check-----
//This is a very good idea to wrap your calls with that function.. Otherwise you will not be able to see what is the error.
//Moreover, you may also want to look at how to use cuda-memcheck and cuda-gdb for debugging.

__global__ void rowSum(int *adj, int *xadj, int* tadj, int* txadj,double *rv, double *cv,int *nov,double *errorOut)
{
	int step = 32*1024;
	int  index = blockDim.x * blockIdx.x + threadIdx.x;
        for(int i = index; i<(*nov)-1; i+=step){
                double  rsum = 0;
                int  row_start  = xadj[i];
                int  row_end    = xadj[i +1];
                for (int jj = row_start; jj < row_end; jj++)
                        rsum += cv[adj[jj]];
                if(rsum != 0)
                {
                        rv[i] = 1.0 / rsum;
                }

        }

}

__global__ void colSum(int *tadj, int *txadj, double *rv, double *cv,int *nov)
{
	int step = 32*1024;
	int  index = blockDim.x * blockIdx.x + threadIdx.x;
	for(int i = index; i<(*nov)-1; i+=step){
                double  csum = 0;
                int  col_start  = txadj[i];
                int  col_end    = txadj[i +1];
                for (int jj = col_start; jj < col_end; jj++)
                        csum += rv[tadj[jj]];
                if(csum != 0)
                {
                        cv[i] = 1.0 / csum;
                }

        }
        __syncthreads();

}

__global__ void errorCheck(int *adj, int *xadj, double *rv, double *cv, int *nov, double *errorOut)
{
	int step = 32*1024;
	__shared__ double sharedMem[1024];
	int  index = blockDim.x * blockIdx.x + threadIdx.x;
	for(int i = index; i<(*nov)-1; i+=step)
        {
                double errorSum = 0;
                int row_start = xadj[i];
                for(int jj = row_start; jj<xadj[i+1]; jj++)
                        errorSum += (cv[adj[jj]] * rv[i]);
                double errorVal = fabsf(1.0 - errorSum);
                sharedMem[threadIdx.x] = sharedMem[threadIdx.x] < errorVal ? errorVal :sharedMem[threadIdx.x];
        }
		__syncthreads();
	int ib = blockDim.x / 2;
	while(ib != 0)
	{
		if(threadIdx.x < ib && sharedMem[threadIdx.x + ib] > sharedMem[threadIdx.x])
		{
			sharedMem[threadIdx.x] = sharedMem[threadIdx.x + ib];
		}
		__syncthreads();
		ib /=2;
	}
	if(threadIdx.x == 0)
		errorOut[blockIdx.x] = sharedMem[0];
	sharedMem[threadIdx.x] = 0;	

}


void wrapper(int* adj, int* xadj, int* tadj, int*txadj, double* rv, double* cv, int* nov, int* nnz, int siter){
  
  printf("Wrapper here! \n");
  
  int *d_adj, *d_xadj, *d_tadj, *d_txadj, *d_nov, *d_nnz, *d_siter;
  double *d_rv,*d_cv, *d_errorOut, *d_errorSortedOut, *errorOut;

  double errorValue;  
  
  errorOut = (double*)malloc(512 * sizeof(double));

  cudaMalloc(&d_adj, *nnz * sizeof(int));
  cudaMalloc(&d_tadj, *nnz * sizeof(int));
  cudaMalloc(&d_xadj, *nov * sizeof(int)); 
  cudaMalloc(&d_txadj, *nov * sizeof(int));

  cudaMalloc(&d_nov, sizeof(int));
  cudaMalloc(&d_nnz, sizeof(int));
  cudaMalloc(&d_siter, sizeof(int));

  cudaMalloc(&d_rv, *nov * sizeof(double));
  cudaMalloc(&d_cv, *nov * sizeof(double));	

  cudaMalloc(&d_errorOut, *nov * sizeof(double));
  cudaMalloc(&d_errorSortedOut, *nov * sizeof(double));

  cudaMemcpy(d_adj, adj, *nnz*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_tadj, tadj, *nnz*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_xadj, xadj, *nov*sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_txadj, txadj, *nov*sizeof(int), cudaMemcpyHostToDevice);
   
  cudaMemcpy(d_nov, nov, sizeof(int), cudaMemcpyHostToDevice);
  cudaMemcpy(d_nnz, nnz, sizeof(int), cudaMemcpyHostToDevice);
  //cudaMemcpy(d_siter, &siter, sizeof(int), cudaMemcpyHostToDevice);

  cudaMemcpy(d_rv, rv, *nov*sizeof(double), cudaMemcpyHostToDevice);
  cudaMemcpy(d_cv, cv, *nov*sizeof(double), cudaMemcpyHostToDevice);
 
  int size = (*nov * sizeof(double));
  errorOut = (double*)malloc(size);
  cudaEvent_t start,stop;
  float elapsedTime;
  cudaEventCreate(&start);
  cudaEventRecord(start, 0);
 
  for(int i = 0; i<siter; i++)
  {
	rowSum<<<32,1024>>>(d_adj, d_xadj,d_tadj,d_txadj, d_rv, d_cv, d_nov,d_errorOut);
  	gpuErrchk( cudaDeviceSynchronize() );

  	colSum<<<32,1024>>>(d_tadj,d_txadj, d_rv, d_cv, d_nov);
  	gpuErrchk( cudaDeviceSynchronize() );

  	errorCheck<<<32,1024>>>(d_adj, d_xadj, d_rv, d_cv, d_nov, d_errorOut);
  	gpuErrchk( cudaDeviceSynchronize() );
	
  	cudaMemcpy(errorOut, d_errorOut, 32*sizeof(double), cudaMemcpyDeviceToHost);
	errorValue = 0;
	for(int k = 0; k<32; k++)
	{
		if(errorValue < errorOut[k])
		{
			errorValue = errorOut[k];
		}	
	}

  	printf("iter-%i error : %f \n\n",i ,errorValue);
	cudaMemset(d_errorOut, 0, *nov*sizeof(double));
  } 
  cudaEventCreate(&stop);
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  
  cudaEventElapsedTime(&elapsedTime, start, stop);
  printf("GPU scale took: %f s\n", elapsedTime/1000);
  free(errorOut);
  cudaFree(d_adj);
  cudaFree(d_xadj);
  cudaFree(d_tadj);
  cudaFree(d_txadj);
  cudaFree(d_errorOut);
  cudaFree(d_errorSortedOut);
  cudaFree(d_rv);
  cudaFree(d_cv);
  cudaFree(d_nov);
  cudaFree(d_nnz);
    
}


