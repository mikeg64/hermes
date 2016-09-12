// includes, system
#include <stdio.h>
#include <assert.h>
 
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

// Simple utility function to check for CUDA runtime errors
void checkCUDAError(const char* msg);
 
// Part3: implement the kernel




__global__ void max_parallel(double *cmax, double *temp,int ndimp, double maxac)
{

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int tid = threadIdx.x;
  extern __shared__ double partialResult[];

  int i;
   partialResult[tid]=0.0;
   if(iindex<ndimp)
              partialResult[tid]=temp[iindex];
  __syncthreads();

      //   if(temp[iindex]==maxac)
      //      printf("max here %d %d \n",tid,blockIdx.x);



     //if(tid==0)
    //     printf("sero %d\n",blockIdx.x);
for(unsigned int s=1; s < blockDim.x; s *= 2) {


        if ((tid % (2*s)) == 0) {
            if(partialResult[tid+s]>partialResult[tid])
                 partialResult[tid]=partialResult[tid + s];
        }
        __syncthreads();
    }

    __syncthreads();
    if(tid==0)
    {
      cmax[blockIdx.x]=partialResult[0];
      //temp[blockIdx.x]=partialResult[0];
     }
    __syncthreads();

}




/////////////////////////////////////////////////////////////////////
// Program main
/////////////////////////////////////////////////////////////////////
int main( int argc, char** argv) 
{
    // pointer for host memory and size
    int *h_a;
    double *h_c, *h_temp;
    double maxc=-1.0;
    double tmax=-1.0;
    int dimA = 256 * 1024; // 256K elements (1MB total)
    dimA=256*2048;
    dimA=2097152;
    // pointer for device memory
    int *d_b, *d_a;

    double *d_c, *d_temp;
 
    // define grid and block size
    int numThreadsPerBlock = 128;
 
    // Part 1: compute number of blocks needed based on 
    // array size and desired block size
    int numBlocks = dimA / numThreadsPerBlock;  
 
    srand (time(NULL));

    // allocate host and device memory
    size_t memSize = numBlocks * numThreadsPerBlock * sizeof(int);
    h_a = (int *) malloc(memSize);
    cudaMalloc( (void **) &d_a, memSize );
    cudaMalloc( (void **) &d_b, memSize );

     int smemSize = numThreadsPerBlock * sizeof(double);

    size_t dmemSize = numBlocks * numThreadsPerBlock * sizeof(double);
    h_c = (double *) malloc(dmemSize);
    h_temp = (double *) malloc(dmemSize);

    cudaMalloc( (void **) &d_c, dmemSize );
    cudaMalloc( (void **) &d_temp, dmemSize );

    int imax;

    int ccount=0;
    int j=0;
   // for( j=0; j<1000; j++)
   // {
   // tmax=-1;

    // Initialize input array on host
    for (int i = 0; i < dimA; ++i)
    {
        h_a[i] = i;
        h_c[i]=(rand()%100000000);

        if(h_c[i]>tmax)
        {
            tmax=h_c[i];
            imax=i;
        }
        //printf(" %g ",h_c[i]);
    }
    printf("\n\n\n %d %f %d\n", dimA, tmax, imax);
 
    // Copy host array to device array
    cudaMemcpy( d_a, h_a, memSize, cudaMemcpyHostToDevice );
    cudaMemcpy( d_c, h_c, dmemSize, cudaMemcpyHostToDevice );
    cudaMemcpy( d_temp, h_c, dmemSize, cudaMemcpyHostToDevice );
 
 
    // device to host copy
    cudaMemcpy( h_a, d_b, memSize, cudaMemcpyDeviceToHost );
 
    // Check for any CUDA errors
    checkCUDAError("memcpy");
 


	  for(int i=0;i<numBlocks;i++)
		       h_temp[i]=0;
	  cudaMemcpy(d_temp, h_temp, numBlocks*sizeof(double), cudaMemcpyHostToDevice);

	  max_parallel<<<numBlocks,numThreadsPerBlock,smemSize>>>(d_temp,d_c,dimA,tmax);
	  cudaThreadSynchronize();
	  cudaMemcpy(h_temp, d_temp, numBlocks*sizeof(double), cudaMemcpyDeviceToHost);

         for(int i=0;i<numBlocks;i++)
         {          		
                if(h_temp[i]>maxc) maxc=h_temp[i]; 
                //printf(" %f ",h_temp[i]);

         }

       if(maxc==tmax) ccount++;

         printf("\n\n\nnumblocks %d %d max=%f %f %d\n",j, numBlocks, maxc, tmax, ccount);



 //     }

  


       
 
    // check if kernel execution generated an error
    // Check for any CUDA errors
    checkCUDAError("kernel invocation");

















 
    // free device memory
    cudaFree(d_a);
    cudaFree(d_b);
 
    // free host memory
    free(h_a);
 
    // If the program makes it this far, then the results are 
    // correct and there are no run-time errors.  Good work!
    printf("Correct!\n");
 
    return 0;
}
void checkCUDAError(const char *msg)
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err) 
    {
        fprintf(stderr, "Cuda error: %s: %s.\n", msg, 
                                  cudaGetErrorString( err) );
        exit(EXIT_FAILURE);
    }                         
}

