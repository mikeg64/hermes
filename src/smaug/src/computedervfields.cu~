//#define MODID pre


#include "../include/cudapars.h"
#include "../include/paramssteeringtest1.h"
#include "../include/iobparams.h"
/////////////////////////////////////
// standard imports
/////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include "../include/smaugcukernels.h"
#include "../include/gradops_cdf.cuh"
#include "../include/dervfields_cdf.cuh"
/////////////////////////////////////
// kernel function (CUDA device)
/////////////////////////////////////
__global__ void computevels_parallel(struct params *p,   real *wmod, real *wd, int order, int dir)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
//  real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     









  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))



     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
	{		               

                        switch(dir)
                        {
                         case 0:
                          #ifdef USE_SAC_3D
       				if(ii[0]<p->n[0] && ii[1]>1 && ii[1]<(p->n[1]-2) && ii[2]>1 && ii[2]<(p->n[2]-2))
     			  #else
       				if(ii[0]<p->n[0] && ii[1]>1 && ii[1]<(p->n[1]-2))
     			  #endif
                         //if(i<(ni)  && j >1 &&  j<(nj-1))
                                           computevel3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
                         break;
                         case 1:
                          #ifdef USE_SAC_3D
       				if(ii[1]<p->n[1] && ii[0]>1 && ii[0]<(p->n[0]-2) && ii[2]>1 && ii[2]<(p->n[2]-2))
     			  #else
       				if(ii[1]<p->n[1] && ii[0]>1 && ii[0]<(p->n[0]-2))
     			  #endif
                         //if(i>1 &&  i<(ni-1) && j<(nj))
                                           computevel3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
                         break;
                          #ifdef USE_SAC_3D
                         case 2:

       				if(ii[2]<p->n[2] && ii[0]>1 && ii[0]<(p->n[0]-2) && ii[1]>1 && ii[1]<(p->n[1]-2))

                         //if(i>1 &&  i<(ni-1) && j<(nj))
                                           computevel3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
                         break;
                         #endif
                        }


         }


              __syncthreads();











  
}


__global__ void computept_parallel(struct params *p,   real *wmod, real *wd, int order, int dir)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
//  real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     




     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
	{		               

                        switch(dir)
                        {
                         case 0:
                          #ifdef USE_SAC_3D
       				if(ii[0]<p->n[0] && ii[1]>1 && ii[1]<(p->n[1]-2) && ii[2]>1 && ii[2]<(p->n[2]-2))
     			  #else
       				if(ii[0]<p->n[0] && ii[1]>1 && ii[1]<(p->n[1]-2))
     			  #endif
                         //if(i<(ni)  && j >1 &&  j<(nj-1))
                                           
                           computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
                         break;
                         case 1:
                          #ifdef USE_SAC_3D
       				if(ii[1]<p->n[1] && ii[0]>1 && ii[0]<(p->n[0]-2) && ii[2]>1 && ii[2]<(p->n[2]-2))
     			  #else
       				if(ii[1]<p->n[1] && ii[0]>1 && ii[0]<(p->n[0]-2))
     			  #endif
                         //if(i>1 &&  i<(ni-1) && j<(nj))
                                           
                            computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
                         break;
                          #ifdef USE_SAC_3D
                         case 2:

       				if(ii[2]<p->n[2] && ii[0]>1 && ii[0]<(p->n[0]-2) && ii[1]>1 && ii[1]<(p->n[1]-2))

                         //if(i>1 &&  i<(ni-1) && j<(nj))
                                          
                                computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
                         break;
                         #endif
                        }


         }




  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))



  /*   ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
	{		               

	     #ifdef ADIABHYDRO
	       
	       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
	     #else
	       
	       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
	     #endif */        
              /* switch(dir)
                        {
                         case 0:
                          #ifdef USE_SAC_3D
       				if(ii[0]<p->n[0] && ii[1]>1 && ii[1]<(p->n[1]-2) && ii[2]>1 && ii[2]<(p->n[2]-2))
     			  #else
       				if(ii[0]<p->n[0] && ii[1]>1 && ii[1]<(p->n[1]-2))
     			  #endif
				     {
				     #ifdef ADIABHYDRO
				       computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				     #else
				       computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				     #endif
				     }
                         break;
                         case 1:
                          #ifdef USE_SAC_3D
       				if(ii[1]<p->n[1] && ii[0]>1 && ii[0]<(p->n[0]-2) && ii[2]>1 && ii[2]<(p->n[2]-2))
     			  #else
       				if(ii[1]<p->n[1] && ii[0]>1 && ii[0]<(p->n[0]-2))
     			  #endif
				     {
				     #ifdef ADIABHYDRO
				       computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				     #else
				       computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				     #endif
				     }
                         break;
                          #ifdef USE_SAC_3D
                         case 2:

       				if(ii[2]<p->n[2] && ii[0]>1 && ii[0]<(p->n[0]-2) && ii[1]>1 && ii[1]<(p->n[1]-2))
				     {

				       computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);

				     }
                         break;
                         #endif
                        }*/


        /* }*/


              __syncthreads();











  
}


__global__ void computeptzero_parallel(struct params *p,   real *wmod, real *wd, int order, int dir)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
//  real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     




     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
	{		               

                        wd[fencode3_cdf(p,ii,pressuret)]=0.0;
                        


         }




  


              __syncthreads();











  
}



__global__ void computepk_parallel(struct params *p,   real *wmod, real *wd, int order, int dir)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
//  real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     









  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))



     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
	{		               

	     #ifdef ADIABHYDRO
	       
	       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
	     #else
	       
	       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
	     #endif         
              /* switch(dir)
                        {
                         case 0:
                          #ifdef USE_SAC_3D
       				if(ii[0]<p->n[0] && ii[1]>1 && ii[1]<(p->n[1]-2) && ii[2]>1 && ii[2]<(p->n[2]-2))
     			  #else
       				if(ii[0]<p->n[0] && ii[1]>1 && ii[1]<(p->n[1]-2))
     			  #endif
				     {
				     #ifdef ADIABHYDRO
				       computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				     #else
				       computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				     #endif
				     }
                         break;
                         case 1:
                          #ifdef USE_SAC_3D
       				if(ii[1]<p->n[1] && ii[0]>1 && ii[0]<(p->n[0]-2) && ii[2]>1 && ii[2]<(p->n[2]-2))
     			  #else
       				if(ii[1]<p->n[1] && ii[0]>1 && ii[0]<(p->n[0]-2))
     			  #endif
				     {
				     #ifdef ADIABHYDRO
				       computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				     #else
				       computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				     #endif
				     }
                         break;
                          #ifdef USE_SAC_3D
                         case 2:

       				if(ii[2]<p->n[2] && ii[0]>1 && ii[0]<(p->n[0]-2) && ii[1]>1 && ii[1]<(p->n[1]-2))
				     {

				       computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);

				     }
                         break;
                         #endif
                        }*/


         }


              __syncthreads();











  
}


__global__ void computepbg_parallel(struct params *p,   real *wmod, real *wd, int order, int dir)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
//  real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

  int ii[NDIM];
  int dimp;
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#else
   dimp=((p->n[0]))*((p->n[1]));
#endif  
   int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     









  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))



     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
	{		               

	     #ifdef ADIABHYDRO
	       computepbg3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
	       
	     #else
	       computepbg3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
	       
	     #endif         
              /* switch(dir)
                        {
                         case 0:
                          #ifdef USE_SAC_3D
       				if(ii[0]<p->n[0] && ii[1]>1 && ii[1]<(p->n[1]-2) && ii[2]>1 && ii[2]<(p->n[2]-2))
     			  #else
       				if(ii[0]<p->n[0] && ii[1]>1 && ii[1]<(p->n[1]-2))
     			  #endif
				     {
				     #ifdef ADIABHYDRO
				       computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				     #else
				       computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				     #endif
				     }
                         break;
                         case 1:
                          #ifdef USE_SAC_3D
       				if(ii[1]<p->n[1] && ii[0]>1 && ii[0]<(p->n[0]-2) && ii[2]>1 && ii[2]<(p->n[2]-2))
     			  #else
       				if(ii[1]<p->n[1] && ii[0]>1 && ii[0]<(p->n[0]-2))
     			  #endif
				     {
				     #ifdef ADIABHYDRO
				       computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				     #else
				       computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				     #endif
				     }
                         break;
                          #ifdef USE_SAC_3D
                         case 2:

       				if(ii[2]<p->n[2] && ii[0]>1 && ii[0]<(p->n[0]-2) && ii[1]>1 && ii[1]<(p->n[1]-2))
				     {

				       computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
				       computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);

				     }
                         break;
                         #endif
                        }*/


         }


              __syncthreads();











  
}


__global__ void computemaxc_parallel(struct params *p,   real *wmod, real *wd, int order, int dir)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
//  real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     




   /*for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   #ifdef USE_SAC_3D
     for(kpg=0;kpg<(p->npgp[2]);kpg++)
   #endif
   {

     ii[0]=ip*(p->npgp[0])+ipg;
     ii[1]=jp*(p->npgp[1])+jpg;
     #ifdef USE_SAC_3D
	   ii[2]=kp*(p->npgp[2])+kpg;
     #endif

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
	{
 //determin cmax
               computec3_cdf(wmod+(order*dimp*NVAR),wd,p,ii,dir);
               //p->cmax=0.0;
        }

}*/
              __syncthreads();



if(iindex==0)
{
   
 //  for(ipg=0;ipg<(p->npgp[0]);ipg++)
 //  for(jpg=0;jpg<(p->npgp[1]);jpg++)
  // {

  //   i=ip*(p->npgp[0])+ipg;
 //    j=jp*(p->npgp[1])+jpg;
   //if( i<((p->n[0])) && j<((p->n[1])))
  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
    p->cmax=0.0;
    for(ii[0]=0;ii[0]<((p->n[0]));ii[0]++)
      for(ii[1]=0;ii[1]<((p->n[1]));ii[1]++)
     #ifdef USE_SAC_3D
        for(ii[2]=0;ii[2]<((p->n[2]));ii[2]++)
     #endif
	{ 
               computecmax3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);




	}

 //  }
}
 __syncthreads(); 

//p->cmax=1.0;



  
}

//from http://www.nvidia.com/object/cuda_sample_data-parallel.html#reduction
/*
    This version uses n/2 threads --
    it performs the first level of reduction when reading from global memory
*/
__global__ void fastcomputemaxc_parallel(struct params *p,   real *wmod, real *wd, int order, int dir)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int tid = threadIdx.x;
  int i,j;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp,ipg,jpg;
        extern __shared__ real sdata[];
  // __shared__ float sdata[];
 //real sdata[dimp];
  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     



    // perform first level of reduction,
    // reading from global memory, writing to shared memory
   sdata[tid]=0.0;


     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
	{
 //determin cmax
               //computec3_cdf(wmod+(order*dimp*NVAR),wd,p,ii,dir);
               //p->cmax=0.0;
               
               if(wd[fencode3_cdf(p,ii,cfast)]>sdata[tid])
                    sdata[tid]=wd[fencode3_cdf(p,ii,cfast)];
        }


              __syncthreads();

    // do reduction in shared mem
    for(unsigned int s=blockDim.x/2; s>0; s>>=1) {
        if (tid < s) {
            //sdata[tid] += sdata[tid + s];
            //if(sdata[tid]>sdata[0])
             //   sdata[0]=sdata[tid];
            if(sdata[tid+s]>sdata[tid])
                sdata[tid]=sdata[tid+s];
        }
        __syncthreads();
    }


    if (tid == 0) p->cmax = sdata[0];

  /* for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   #ifdef USE_SAC_3D
     for(kpg=0;kpg<(p->npgp[2]);kpg++)
   #endif
   {

     ii[0]=ip*(p->npgp[0])+ipg;
     ii[1]=jp*(p->npgp[1])+jpg;
     #ifdef USE_SAC_3D
	   ii[2]=kp*(p->npgp[2])+kpg;
     #endif

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
	{
 //determin cmax
               //computec3_cdf(wmod+(order*dimp*NVAR),wd,p,ii,dir);
               //p->cmax=0.0;
               
               if(wd[fencode3_cdf(p,ii,cfast)]>(p->cmax))
                    sdata[tid]=wd[fencode3_cdf(p,ii,cfast)];
        }

}
              __syncthreads();*/

    // do reduction in shared mem
    /*for(unsigned int s=blockDim.x/2; s>0; s>>=1) {
        if (tid < s) {
            //sdata[tid] += sdata[tid + s];
            //if(sdata[tid]>sdata[0])
             //   sdata[0]=sdata[tid];
            if(sdata[tid+s]>sdata[tid])
                sdata[tid]=sdata[tid+s];
        }
        __syncthreads();
    }


    if (tid == 0) p->cmax = sdata[0];*/

  
}



/*
    Parallel sum reduction using shared memory
    - takes log(n) steps for n input elements
    - uses n threads
    - only works for power-of-2 arrays
*/

/* This reduction interleaves which threads are active by using the modulo
   operator.  This operator is very expensive on GPUs, and the interleaved 
   inactivity means that no whole warps are active, which is also very 
   inefficient */


//from http://www.nvidia.com/object/cuda_sample_data-parallel.html#reduction
/* This reduction interleaves which threads are active by using the modulo
   operator.  This operator is very expensive on GPUs, and the interleaved 
   inactivity means that no whole warps are active, which is also very 
   inefficient */
__global__ void reduction0computemaxc_parallel(struct params *p,   real *wmod, real *wd, int order, int dir)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int tid = threadIdx.x;
  int i,j;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp;
        extern __shared__ real sdata[];
 
  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif

int tnumThreadsPerBlock = 128;
    
int numBlocks = (dimp+tnumThreadsPerBlock-1) / tnumThreadsPerBlock;
  real temp[1024];
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
   //sdata[tid]=0.0;
    if(iindex<1024)
      temp[iindex]=0.0;

     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif

   /*  #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif*/
  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
//	{
 //determin cmax
               //computec3_cdf(wmod+(order*dimp*NVAR),wd,p,ii,dir);
               //p->cmax=0.0;
               
              // if(wd[fencode3_cdf(p,ii,cfast)]>(p->cmax))


         if(iindex<dimp)
                    sdata[tid]=wd[fencode3_cdf(p,ii,cfast)];

              __syncthreads();


    // do reduction in shared mem
    for(unsigned int s=1; s < blockDim.x; s *= 2) {



        // modulo arithmetic is slow!
        if ((tid % (2*s)) == 0) {
            if(sdata[tid+s]>sdata[tid])
                 sdata[tid]=sdata[tid + s];
            
        }
        // strided indexing using sequential addressing is better!
        /*int tindex=2*s*tid;
        if (tindex<blockDim.x) {
            if(sdata[tid+s]>sdata[tid])
                 sdata[tid]=sdata[tid + s];
        }
        __syncthreads();*/
         __syncthreads();
    }

    __syncthreads();
    if(tid==0)
      temp[blockIdx.x]=sdata[0];
__syncthreads();
    if(iindex==0)
       for(int i=0; i<numBlocks; i++)
         if(temp[i]>(p->cmax)) p->cmax=temp[i];
     if (tid == 0 && p->cmax<sdata[0] ) p->cmax = sdata[0];
 


    /* ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
	{
 //determin cmax
               //computec3_cdf(wmod+(order*dimp*NVAR),wd,p,ii,dir);
               //p->cmax=0.0;
               
               if(wd[fencode3_cdf(p,ii,cfast)]>(p->cmax))
                    p->cmax=wd[fencode3_cdf(p,ii,cfast)];
        }


              __syncthreads();*/
 
//        }

//}
//p->cmax=1.0;
 
}




__global__ void newreduction0computemax_parallel(real *cmax, real *temp,int ndimp)
{
  //real *cmax, real *temp, int ndimp

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int tid = threadIdx.x;
  extern __shared__ double partialResult[];

  int i;
   partialResult[tid]=0.0;
   if(iindex<ndimp)
              partialResult[tid]=temp[iindex];
  __syncthreads();


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



__global__ void myreduction0computemaxc_parallel(struct params *p,   real *wmod, real *wd, int order, int dir, real *temp,int ndimp,int s)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int tid = threadIdx.x;
  int i,j;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp;
//        extern __shared__ real sdata[];
 
  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif

//int tnumThreadsPerBlock = numThreadsPerBlock;
//  int tnumThreadsPerBlock = 128; 
//int numBlocks = (dimp+tnumThreadsPerBlock-1) / tnumThreadsPerBlock;
  //real temp[dimp];
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
   //sdata[tid]=0.0;
   // if(iindex<1024)
    //  temp[iindex]=0.0;

     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif
    //int s=1;


   
    //while(((s*=2)<=((ndimp/2)-1)) && ((iindex+s)<ndimp)) {
    if((iindex+s)<ndimp)
            if(temp[iindex+s]>temp[iindex])
                 temp[iindex]=temp[iindex + s];
            
       // }

       //  __syncthreads();
    

   // __syncthreads();

   if(iindex==0)
      p->cmax=temp[0];


 
}



__global__ void myreduction0computemaxcourant_parallel(struct params *p,   real *wmod, real *wd, int order, int dir, real *temp,int ndimp,int s)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int tid = threadIdx.x;
  int i,j;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp;
//        extern __shared__ real sdata[];
 
  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif

//int tnumThreadsPerBlock = 128;
    
//int numBlocks = (dimp+tnumThreadsPerBlock-1) / numThreadsPerBlock;
  //real temp[dimp];
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
   //sdata[tid]=0.0;
   // if(iindex<1024)
    //  temp[iindex]=0.0;

     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif
    //int s=1;


   
    //while(((s*=2)<=((ndimp/2)-1)) && ((iindex+s)<ndimp)) {
    if((iindex+s)<ndimp)
            if(temp[iindex+s]>temp[iindex])
                 temp[iindex]=temp[iindex + s];
            
       // }

       //  __syncthreads();
    

   // __syncthreads();

   if(iindex==0  && (p->maxcourant<temp[0]))
      p->maxcourant=temp[0];


 
}

__global__ void copytotemp_parallel(struct params *p, real *wd, real *wtemp,int ifield)
{
  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int tid = threadIdx.x;
  int i,j;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp;

 
  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif


     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif


     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
	{
 
               
               
                    wtemp[fencode3_cdf(p,ii,0)]=wd[fencode3_cdf(p,ii,ifield)];
        }


              __syncthreads();







}


__global__ void zeropadmaxc_parallel(struct params *p,   real *wmod, real *wd, int order, int dir, real *temp, int ndimp)
{
  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
 
  if(iindex<ndimp)
      temp[iindex]=0.0;

}

__global__ void zeropadmaxcourant_parallel(struct params *p,   real *wmod, real *wd, int order, int dir, real *temp, int ndimp)
{
  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
 
  //if(iindex<ndimp)
  //    temp[iindex]=0.0;

  unsigned int tid = threadIdx.x;
  int i,j;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp;
//        extern __shared__ real sdata[];
 
  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif

    //if(iindex<ndimp)
    //  temp[iindex]=0.0;
  
//int numBlocks = (dimp+tnumThreadsPerBlock-1) / tnumThreadsPerBlock;
  //real temp[dimp];
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
   //sdata[tid]=0.0;
   // if(iindex<1024)
    //  temp[iindex]=0.0;

     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif
    //int s=1;

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
             temp[iindex]=temp[iindex]/(wd[fencode3_cdf(p,ii,delx1+dir)]);
       else
            temp[index]=0.0;







}

__global__ void reduction0computemaxcfast_parallel(struct params *p,   real *wmod, real *wd, int order, int dir)
{


 // int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int tid = threadIdx.x;
  int i,j;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];

    
    unsigned int iindex = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp;
    extern __shared__ real sdata[];



    // perform first level of reduction,
    // reading from global memory, writing to shared memory
   sdata[tid]=0.0;

                   sdata[tid]=wd[blockDim.x+(cfast*dimp)+iindex];

              __syncthreads();


    // do reduction in shared mem
    for(unsigned int s=1; s < blockDim.x; s *= 2) {
        // modulo arithmetic is slow!
        if ((tid % (2*s)) == 0) {
            if(sdata[tid+s]>sdata[tid])
                 sdata[tid]=sdata[tid + s];
        }
        __syncthreads();
    }

    if (tid == 0 && (p->cmax)<sdata[0] ) p->cmax = sdata[0];
 __syncthreads();



 
}




//from http://www.nvidia.com/object/cuda_sample_data-parallel.html#reduction
/* This reduction interleaves which threads are active by using the modulo
   operator.  This operator is very expensive on GPUs, and the interleaved 
   inactivity means that no whole warps are active, which is also very 
   inefficient */
__global__ void reductiona0computemaxc_parallel(struct params *p,   real *wmod, real *wd, int order, int dir)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int tid = threadIdx.x;
  int i,j;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp;
        extern __shared__ real sdata[];
   //__shared__ float sdata[];
 //real sdata[dimp];
  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     



    // perform first level of reduction,
    // reading from global memory, writing to shared memory
   sdata[tid]=0.0;

   if(iindex<dimp)
      sdata[tid]=wd[iindex+(dimp*cfast)];

       /* if(iindex<dimp)
               if(wd[iindex+(dimp*cfast)]>(p->cmax))
                    sdata[tid]=wd[iindex+(dimp*cfast)];*/

              __syncthreads();


    // do reduction in shared mem
    for(unsigned int s=1; s < blockDim.x; s *= 2) {
        // modulo arithmetic is slow!
        if ((tid % (2*s)) == 0) {
            if(sdata[tid+s]>sdata[tid])
                 sdata[tid]=sdata[tid + s];
        }
        __syncthreads();
    }

    if (tid == 0) p->cmax = sdata[0];



     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
	{
 //determin cmax
               //computec3_cdf(wmod+(order*dimp*NVAR),wd,p,ii,dir);
               //p->cmax=0.0;
               
               if(wd[fencode3_cdf(p,ii,cfast)]>(p->cmax))
                    p->cmax=wd[fencode3_cdf(p,ii,cfast)];
        }


              __syncthreads();
  
}


__global__ void computec_parallel(struct params *p,   real *wmod, real *wd, int order, int dir)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
//  real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     






 p->cmax=0.0;


     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
	{
 //determin cmax
               computec3_cdf(wmod+(order*dimp*NVAR),wd,p,ii,dir);
               //p->cmax=0.0;
        }


              __syncthreads();












  
}


__global__ void computedervfields_parallel(struct params *p,   real *wmod, real *wd, int order)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
//  real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     
int dir=0;

 

/*                if(jp==0 && ip==63)
                {
     ii[0]=ip;
     ii[1]=jp;

             // printf("density at 128,128=%12.10f %d\n", wmod[fencode3_cdf(p,ii,rho)], order);    
          printf("density at 128,128=%d %d %12.10f %12.10f  \n",ip,jp, wd[fencode3_cdf(p,ii,delx1)], wd[fencode3_cdf(p,ii,delx2)]);          
}

                if(jp==2 && ip==63)
                {
     ii[0]=ip;
     ii[1]=jp;

             // printf("density at 128,128=%12.10f %d\n", wmod[fencode3_cdf(p,ii,rho)], order);    
          printf("density at 128,128=%d %d %12.10f %12.10f  \n",ip,jp, wd[fencode3_cdf(p,ii,delx1)], wd[fencode3_cdf(p,ii,delx2)]);          
}



                if(jp==255 && ip==63)
                {
     ii[0]=ip;
     ii[1]=jp;

             // printf("density at 128,128=%12.10f %d\n", wmod[fencode3_cdf(p,ii,rho)], order);    
          printf("density at 128,128=%d %d %12.10f %12.10f  \n", ip,jp,wd[fencode3_cdf(p,ii,delx1)], wd[fencode3_cdf(p,ii,delx2)]);          
}



                if(jp==253 && ip==63)
                {
     ii[0]=ip;
     ii[1]=jp;

             // printf("density at 128,128=%12.10f %d\n", wmod[fencode3_cdf(p,ii,rho)], order);    
          printf("density at 128,128=%d %d %12.10f %12.10f  \n", ip,jp,wd[fencode3_cdf(p,ii,delx1)], wd[fencode3_cdf(p,ii,delx2)]);          
}






                if(jp==63 && ip==0)
                {
     ii[0]=ip;
     ii[1]=jp;

             // printf("density at 128,128=%12.10f %d\n", wmod[fencode3_cdf(p,ii,rho)], order);    
          printf("density at 128,128=%d %d %12.10f %12.10f  \n",ip,jp, wd[fencode3_cdf(p,ii,delx1)], wd[fencode3_cdf(p,ii,delx2)]);          
}



                if(jp==63 && ip==2)
                {
     ii[0]=ip;
     ii[1]=jp;

             // printf("density at 128,128=%12.10f %d\n", wmod[fencode3_cdf(p,ii,rho)], order);    
          printf("density at 128,128=%d %d %12.10f %12.10f  \n",ip,jp, wd[fencode3_cdf(p,ii,delx1)], wd[fencode3_cdf(p,ii,delx2)]);          
}






                if(jp==63 && ip==255)
                {
     ii[0]=ip;
     ii[1]=jp;

             // printf("density at 128,128=%12.10f %d\n", wmod[fencode3_cdf(p,ii,rho)], order);    
          printf("density at 128,128=%d %d %12.10f %12.10f  \n", ip,jp,wd[fencode3_cdf(p,ii,delx1)], wd[fencode3_cdf(p,ii,delx2)]);          
}




                if(jp==63 && ip==253)
                {
     ii[0]=ip;
     ii[1]=jp;

             // printf("density at 128,128=%12.10f %d\n", wmod[fencode3_cdf(p,ii,rho)], order);    
          printf("density at 128,128=%d %d %12.10f %12.10f  \n", ip,jp,wd[fencode3_cdf(p,ii,delx1)], wd[fencode3_cdf(p,ii,delx2)]);          
}*/







if(order == 0)
{

     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
	{		
               for(int f=vel1; f<=pkb; f++)
                        wd[fencode3_cdf(p,ii,f)]=0; 

                //Here we set the current field values (order==0)to the
                //values which were updated (order=1)
                //hence the dimp*NVAR term in wmod on RHS of expression 
		#ifdef USE_SAC_3D
		  for(int f=rho; f<=b3; f++)
                  	//wmod[fencode3_cdf(p,ii,f)+dimp*NVAR]=wmod[fencode3_cdf(p,ii,f)];                    
                        wmod[fencode3_cdf(p,ii,f)]=wmod[fencode3_cdf(p,ii,f)+dimp*NVAR]; 




		#else
		  for(int f=rho; f<=b2; f++)
                  	//wmod[fencode3_cdf(p,ii,f)+dimp*NVAR]=wmod[fencode3_cdf(p,ii,f)]; 
                         wmod[fencode3_cdf(p,ii,f)]=wmod[fencode3_cdf(p,ii,f)+dimp*NVAR];


                      
		#endif               



// for(int field=rho;field<=rho ; field++)
//if(  (p->ipe)==0  && ((p)->it)==1 && ( isnan(wmod[fencode3_cdf(p,ii,field)])|| wmod[fencode3_cdf(p,ii,field)]==0 ))
//if(  /*(p->ipe)==0  &&*/ (  wmod[fencode3_cdf(p,ii,field)]==0 ))
//       { 
//    				printf("nant %d %d %d %d %lg %lg \n",ii[0],ii[1],field,dir, wmod[fencode3_cdf(p,ii,rho)],wmod[fencode3_cdf(p,ii,field)+dimp*NVAR] );
//;//wmod[fencode3_cdf(p,ii,rho)]=0.221049;
//;//wmod[fencode3_cdf(p,ii,field)+dimp*NVAR]=0.221049;
//}



        }




}

               __syncthreads();



  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))



     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if( ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
	{		               
             #ifdef ADIABHYDRO
               //computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
               //computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
             #else
               //computevel3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
               //computej3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
               //computepk3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
               //computept3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);

               computebdotv3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);
               //computedivb3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);

             #endif

         }


              __syncthreads();

  
}

/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_cdf(char *label)
{
  // we need to synchronise first to catch errors due to
  // asynchroneous operations that would otherwise
  // potentially go unnoticed

  cudaError_t err;

  err = cudaThreadSynchronize();
  if (err != cudaSuccess)
  {
    char *e = (char*) cudaGetErrorString(err);
    fprintf(stderr, "CUDA Error: %s (at %s)", e, label);
  }

  err = cudaGetLastError();
  if (err != cudaSuccess)
  {
    char *e = (char*) cudaGetErrorString(err);
    fprintf(stderr, "CUDA Error: %s (at %s)", e, label);
  }
}




int cucomputedervfields(struct params **p,  struct params **d_p, real **d_wmod,  real **d_wd, int order)
{
  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   ////cudaSetDevice(selectedDevice);
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif  

 dim3 dimBlock(dimblock, 1);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   // dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;

cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);
     computedervfields_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order);

     cudaThreadSynchronize();
 

    cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);


  //checkErrors("copy data from device");


 


}

int cucomputevels(struct params **p,  struct params **d_p, real **d_wmod,  real **d_wd, int order, int dir)
{
  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   ////cudaSetDevice(selectedDevice);
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 

 //dim3 dimBlock(dimblock, 1);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   // dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;


     computevels_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir);

     cudaThreadSynchronize();
 

   // cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);


  //checkErrors("copy data from device");


 


}

int cucomputemaxc(struct params **p,  struct params **d_p, real **d_wmod,  real **d_wd, int order, int dir, real **wd, real **d_wtemp)
{
  int dimp=(((*p)->n[0]))*(((*p)->n[1]));
  
  real fn,fractn,in;
  int ndimp;

  double *d_cmax;
  double *d_bmax;

////cudaSetDevice(selectedDevice);
   int nit=100;
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 

    fn=log(dimp)/log(2.0);
    fractn=modf(fn,&in);
    
  if(fractn>0)
  {
   fn+=1;
   ndimp=(int)pow(2,fn);
  }
  else
   ndimp=dimp;

  //Number threads per block
  //int NTPB=512;
  int NTPB=tnumThreadsPerBlock;
  //Num blocks is determined by size of zeropadded 2^n size array
  int numBlocks = (ndimp+NTPB-1) / NTPB;

  //Shared memory
  int smemSize = NTPB * sizeof(double);

  //Array to store maximum values for reduction in host memory 
  double *h_cmax = (double*)malloc(numBlocks*sizeof(double));

  cudaMalloc((void**)&d_cmax, numBlocks*sizeof(double)); 
  cudaMalloc((void**)&d_bmax, numBlocks*sizeof(double)); //Array to store maximum values for reduction in GPU global memory

  //set maximum value to zero and update values in GPU memory
  (*p)->cmax=0.0;
  cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);
 
  //determine maximum value of magneto-acoustic fast mode
  zeropadmaxc_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir, *d_wtemp,ndimp);

  //the following cudamemcpy routines are bottlenecks!
  //the copytotemp routine does this on the GPU
  //cudaMemcpy(*wd, *d_wd, NDERV*dimp*sizeof(real), cudaMemcpyDeviceToHost);
  //cudaMemcpy(*d_wtemp, ((*wd)+(cfast*dimp)), dimp*sizeof(real), cudaMemcpyHostToDevice);
  //need GPU memory routine here to transfer data from derived memory (fast wave speeds to temporrary area)
  copytotemp_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wd, *d_wtemp,cfast);
  int i=0;


  //find the maximum in each block
  for(i=0;i<numBlocks;i++)
                h_cmax[i]=0;
  cudaMemcpy(d_bmax, h_cmax, numBlocks*sizeof(double), cudaMemcpyHostToDevice);

  newreduction0computemax_parallel<<<numBlocks,NTPB,smemSize>>>(d_bmax,*d_wtemp,ndimp);
  cudaThreadSynchronize();
  cudaMemcpy(h_cmax, d_bmax, numBlocks*sizeof(double), cudaMemcpyDeviceToHost);

  //compare the maxima for all of the blocks and determine maximum value
  for( i=0;i<numBlocks;i++)          		
                if(h_cmax[i]>((*p)->cmax)) ((*p)->cmax)=h_cmax[i];


 //determine maximum value of sound speed
 //cudaMemcpy(*d_wtemp, ((*wd)+(soundspeed*dimp)), dimp*sizeof(real), cudaMemcpyHostToDevice);
 copytotemp_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wd, *d_wtemp,soundspeed);

 zeropadmaxc_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir, *d_wtemp,ndimp);
 for(i=0;i<numBlocks;i++)
                h_cmax[i]=0;
 cudaMemcpy(d_bmax, h_cmax, numBlocks*sizeof(double), cudaMemcpyHostToDevice);

 newreduction0computemax_parallel<<<numBlocks,NTPB,smemSize>>>(d_bmax,*d_wtemp,ndimp);
 cudaThreadSynchronize();

 cudaMemcpy(h_cmax, d_bmax, numBlocks*sizeof(double), cudaMemcpyDeviceToHost);
  //compare the maxima for all of the blocks and determine maximum value
  for( i=0;i<numBlocks;i++)          		
                if(h_cmax[i]>((*p)->cmax)) ((*p)->cmax)=h_cmax[i];

 //((*p)->cmax)=2.0;
  int oldnumBlocks,newnumBlocks;
  newnumBlocks=numBlocks;

  //printf("loop over blocks %d\n\n\n",newnumBlocks);
  /*while(newnumBlocks>1)
  {

      
       //cudaMemcpy(d_wtemp, h_cmax, newnumBlocks*sizeof(double), cudaMemcpyHostToDevice);
    //cudaMemcpy(h_cmax, d_wtemp, newnumBlocks*sizeof(double), cudaMemcpyDeviceToHost);
  //for (i=0; i<oldnumBlocks; i++)
   // {
    //  fprintf(stdout,"cmax# %d %f\n",i, h_cmax[i]);
   // }


	  for( i=0;i<newnumBlocks;i++)
          {
           //printf("gt10 %d %g\n",i, h_cmax[i]);		
        h_cmax[i]=0;
               
          }

	  cudaMemcpy(d_bmax, h_cmax, newnumBlocks*sizeof(double), cudaMemcpyHostToDevice);



       oldnumBlocks=newnumBlocks;
  	newnumBlocks = (newnumBlocks+NTPB-1) / NTPB;
            // printf("blocsk  %d %d\n",newnumBlocks,oldnumBlocks);

  	newreduction0computemax_parallel<<<newnumBlocks,NTPB,smemSize>>>(d_bmax,*d_wtemp,oldnumBlocks);
       cudaThreadSynchronize();
       cudaMemcpy(h_cmax, d_bmax, oldnumBlocks*sizeof(double), cudaMemcpyDeviceToHost);
     //cudaMemcpy(h_cmax, d_wtemp, oldnumBlocks*sizeof(double), cudaMemcpyDeviceToHost);*/
  /*for (i=0; i<oldnumBlocks; i++)
    {
      fprintf(stdout,"cmax# %d %f\n",i, h_cmax[i]);
    }
       fprintf(stdout,"\n");*/


  //}


//(*p)->cmax=h_cmax[0];


//printf("cmax fast=%g\n",h_cmax[0]);


//reduction0computemaxcfast_parallel<<<numBlocks, numThreadsPerBlock,smemSize>>>(*d_p, *d_wmod,  *d_wd, order, dir);
//myreduction0computemaxcfast_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd,*d_wtemp, order, dir);

 //reductiona0computemaxc_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir);
 //  computemaxc_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir);



    // fastcomputemaxc_parallel<<<numBlocks, numThreadsPerBlock,smemSize>>>(*d_p, *d_wmod,  *d_wd, order, dir);
cudaThreadSynchronize();
//cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);

//printf("cmax slow=%g\n",(*p)->cmax);

//(*p)->cmax=2.0;
cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);
//cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);

//printf("cmax on device %.8f\n",(*p)->cmax);
/*(*p)->cmax=0.0;
cudaMemcpy(*wd, *d_wd, NDERV*dimp*sizeof(real), cudaMemcpyDeviceToHost);
for(int i=0; i<dimp;i++)
{
if(((*wd)[i+(soundspeed*dimp)])>((*p)->cmax))
                    (*p)->cmax=(*wd)[i+(soundspeed*dimp)];
if(((*wd)[i+(cfast*dimp)])>((*p)->cmax))
                    (*p)->cmax=(*wd)[i+(cfast*dimp)];
}*/
/*printf("cmax on cpu %.8f\n",(*p)->cmax);*/
//cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);
 /*for(int i=0; i<nit;i++)
{
 reduction0computemaxc_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir);
   // computemaxc_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir);
     cudaThreadSynchronize();
}*/



//    cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);

  // cudaFree(*d_ttemp);
  //checkErrors("copy data from device");


   free(h_cmax);
  cudaFree(d_bmax);
  cudaFree(d_cmax);


}



int cucomputemaxcourant(struct params **p,  struct params **d_p, real **d_wmod,  real **d_wd, int order, int dir, real **wd, real **d_wtemp)
{
  int dimp=(((*p)->n[0]))*(((*p)->n[1]));
    double *d_cmax;
  double *d_bmax;

  real fn,fractn,in;
  int i,ndimp;
////cudaSetDevice(selectedDevice);
   int nit=100;
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 

    fn=log(dimp)/log(2.0);
    fractn=modf(fn,&in);
    
    if(fractn>0)
    {
       fn+=1;
       ndimp=(int)pow(2,fn);
     }
     else
       ndimp=dimp;
       
       int NTPB=tnumThreadsPerBlock;
  int numBlocks = (ndimp+NTPB-1) / NTPB;

  int smemSize = NTPB * sizeof(double);
 double *h_cmax = (double*)malloc(numBlocks*sizeof(double));

  cudaMalloc((void**)&d_cmax, numBlocks*sizeof(double)); 
  cudaMalloc((void**)&d_bmax, numBlocks*sizeof(double)); 


   ((*p)->maxcourant)=0;
   //(*p)->maxcourant=0.0;
   // int smemSize = numThreadsPerBlock * sizeof(real);
  cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);
 //dim3 dimBlock(dimblock, 1);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   // dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   //int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;




//cudaMemcpy(*d_wtemp, *d_wd, NDERV*dimp*sizeof(real), cudaMemcpyDeviceToHost);
//  zeropadmaxc_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir, *d_wtemp,ndimp);
cudaMemcpy(*wd, *d_wd, NDERV*dimp*sizeof(real), cudaMemcpyDeviceToHost);
cudaMemcpy(*d_wtemp, ((*wd)+(cfast*dimp)), dimp*sizeof(real), cudaMemcpyHostToDevice);
 zeropadmaxcourant_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir, *d_wtemp,ndimp);


   for(i=0;i<numBlocks;i++)
                h_cmax[i]=0;
  cudaMemcpy(d_bmax, h_cmax, numBlocks*sizeof(double), cudaMemcpyHostToDevice);


  newreduction0computemax_parallel<<<numBlocks,NTPB,smemSize>>>(d_bmax,*d_wtemp,ndimp);
  cudaThreadSynchronize();
  cudaMemcpy(h_cmax, d_bmax, numBlocks*sizeof(double), cudaMemcpyDeviceToHost);


  int oldnumBlocks,newnumBlocks;
  newnumBlocks=numBlocks;

  /*while(newnumBlocks>1)
  {
        for(i=0;i<numBlocks;i++)
                h_cmax[i]=0;
        cudaMemcpy(d_bmax, h_cmax, numBlocks*sizeof(double), cudaMemcpyHostToDevice);
        for(i=0;i<numBlocks;i++)
                h_cmax[i]=0;
        cudaMemcpy(d_bmax, h_cmax, numBlocks*sizeof(double), cudaMemcpyHostToDevice);


       //cudaMemcpy(d_wtemp, h_cmax, numBlocks*sizeof(double), cudaMemcpyHostToDevice);
       oldnumBlocks=newnumBlocks;
  	newnumBlocks = (newnumBlocks+NTPB-1) / NTPB;

  	newreduction0computemax_parallel<<<newnumBlocks,NTPB,smemSize>>>(d_bmax,*d_wtemp,oldnumBlocks);
       cudaThreadSynchronize();
       cudaMemcpy(h_cmax, d_bmax, newnumBlocks*sizeof(double), cudaMemcpyDeviceToHost);*/

  /*for (i=0; i<numBlocks; i++)
    {
      fprintf(stdout,"cmax# %d %f\n",i, h_cmax[i]);
    }
       fprintf(stdout,"\n");*/


  //}
  cudaMemcpy(h_cmax, d_bmax, numBlocks*sizeof(double), cudaMemcpyDeviceToHost);

  for( i=0;i<numBlocks;i++)          		
                if(h_cmax[i]>((*p)->maxcourant)) ((*p)->maxcourant)=h_cmax[i];







//(*p)->maxcourant=h_cmax[0];
cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);
















/*int s=1;
myreduction0computemaxcourant_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir, *d_wtemp,ndimp,s);
cudaThreadSynchronize();
while(((s*=2)<=((ndimp/2)-1)) ) 
{
   myreduction0computemaxcourant_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir, *d_wtemp,ndimp,s);
   cudaThreadSynchronize();
}*/
//reduction0computemaxcfast_parallel<<<numBlocks, numThreadsPerBlock,smemSize>>>(*d_p, *d_wmod,  *d_wd, order, dir);
//myreduction0computemaxcfast_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd,*d_wtemp, order, dir);

 //reductiona0computemaxc_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir);
  // computemaxc_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir);
    // fastcomputemaxc_parallel<<<numBlocks, numThreadsPerBlock,smemSize>>>(*d_p, *d_wmod,  *d_wd, order, dir);
cudaThreadSynchronize();

//(*p)->cmax=2.0;
//cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);
cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);

//printf("cmax on device %.8f\n",(*p)->cmax);
//(*p)->cmax=0.0;
//cudaMemcpy(*wd, *d_wd, NDERV*dimp*sizeof(real), cudaMemcpyDeviceToHost);
/*for(int i=0; i<dimp;i++)
{

if(((*wd)[i+(cfast*dimp)])>((*p)->cmax))
                    (*p)->cmax=(*wd)[i+(cfast*dimp)];
}
printf("cmax on cpu %.8f\n",(*p)->cmax);*/
//cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);
 /*for(int i=0; i<nit;i++)
{
 reduction0computemaxc_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir);
   // computemaxc_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir);
     cudaThreadSynchronize();
}*/



//    cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);

  // cudaFree(*d_ttemp);
  //checkErrors("copy data from device");


    free(h_cmax);
  cudaFree(d_bmax);
  cudaFree(d_cmax);



}




int cucomputec(struct params **p,  struct params **d_p, real **d_wmod,  real **d_wd, int order, int dir)
{
  int dimp=(((*p)->n[0]))*(((*p)->n[1]));
////cudaSetDevice(selectedDevice);
   
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 

 //dim3 dimBlock(dimblock, 1);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   // dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;


     computec_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir);

     cudaThreadSynchronize();
 

    cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);


  //checkErrors("copy data from device");


 


}

int cucomputept(struct params **p,  struct params **d_p, real **d_wmod,  real **d_wd, int order, int dir)
{

 int dimp=(((*p)->n[0]))*(((*p)->n[1]));
////cudaSetDevice(selectedDevice);
   
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 

 //dim3 dimBlock(dimblock, 1);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   // dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;

    computeptzero_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir);

     cudaThreadSynchronize();
     computept_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir);

     cudaThreadSynchronize();
 

   // cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);


  //checkErrors("copy data from device");


 


}

int cucomputepk(struct params **p,  struct params **d_p, real **d_wmod,  real **d_wd, int order, int dir)
{

 int dimp=(((*p)->n[0]))*(((*p)->n[1]));
////cudaSetDevice(selectedDevice);
   
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 

 //dim3 dimBlock(dimblock, 1);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   // dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;


     computepk_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir);

     cudaThreadSynchronize();
 

   // cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);


  //checkErrors("copy data from device");


 


}


int cucomputepbg(struct params **p,  struct params **d_p, real **d_wmod,  real **d_wd, int order, int dir)
{

 int dimp=(((*p)->n[0]))*(((*p)->n[1]));
////cudaSetDevice(selectedDevice);
   
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 

 //dim3 dimBlock(dimblock, 1);
    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   // dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;


     computepbg_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir);

     cudaThreadSynchronize();
 

   // cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);


  //checkErrors("copy data from device");


 


}







