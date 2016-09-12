#include "../include/cudapars.h"
#include "../include/paramssteeringtest1.h"

/////////////////////////////////////
// standard imports
/////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include "../include/smaugcukernels.h"

/////////////////////////////////////
// kernel function (CUDA device)
/////////////////////////////////////
#include "../include/gradops_db.cuh"
#include "../include/dervfields_db.cuh"


__device__ __host__
real dbsourcerho (real *dw, real *wd, real *w, struct params *p,int *ii) {

  real src=0;

  
 
  return src;
}

__device__ __host__
real dbsourcemom (real *dw, real *wd, real *w, struct params *p,int *ii,int field, int direction) {

  real src=0;
  switch(direction)
  {
	case 0:
         src= -wd[fencode3_db(p,ii,divb)]*w[fencode3_db(p,ii,b1)];
	break;
	case 1:
         src= -wd[fencode3_db(p,ii,divb)]*w[fencode3_db(p,ii,b2)];
	break;
   #ifdef USE_SAC_3D
	case 2:
         src= -wd[fencode3_db(p,ii,divb)]*w[fencode3_db(p,ii,b3)];
	break;
   #endif
  }

  return(isnan(src)?0:src);


}

__device__ __host__
real dbsourceb (real *dw, real *wd, real *w, struct params *p,int *ii,int field, int direction) {

  real src=0;
  switch(direction)
  {
   #ifdef USE_SAC
	case 0:
         src= -wd[fencode3_db(p,ii,divb)]*w[fencode3_db(p,ii,mom1)]/(w[fencode3_db(p,ii,rho)]+w[fencode3_db(p,ii,rhob)]);
	break;
	case 1:
         src= -wd[fencode3_db(p,ii,divb)]*w[fencode3_db(p,ii,mom2)]/(w[fencode3_db(p,ii,rho)]+w[fencode3_db(p,ii,rhob)]);
	break;
   #endif
   #ifdef USE_SAC_3D
	case 2:
         src= -wd[fencode3_db(p,ii,divb)]*w[fencode3_db(p,ii,mom3)]/(w[fencode3_db(p,ii,rho)]+w[fencode3_db(p,ii,rhob)]);
	break;
   #endif
  }
   return(isnan(src)?0:src);
}

__device__ __host__
real dbsourceenergy (real *dw, real *wd, real *w, struct params *p,int *ii) {

 real src=0;
    src= -wd[fencode3_db(p,ii,divb)]*wd[fencode3_db(p,ii,bdotv)];
 
  return ( src);
}


__device__ __host__
int dbderivsourcerho (real *dw, real *wd, real *w, struct params *p,int *ii) {

  int status=0;
  int field=rho;
        dw[fencode3_db(p,ii,field)]=dw[fencode3_db(p,ii,field)]+dbsourcerho(dw,wd,w,p,ii);
     	//dw[fencode3_db(p,ii,field)]=w[fencode3_db(p,ii,field)]+10;
  return ( status);
}

__device__ __host__
int dbderivsourcemom (real *dw, real *wd, real *w, struct params *p,int *ii,int field, int direction) {

  int status=0;
     	//dw[fencode3_db(p,ii,field)]=w[fencode3_db(p,ii,field)]+20+5*(2*direction+1);
        dw[fencode3_db(p,ii,field)]=dw[fencode3_db(p,ii,field)]+dbsourcemom(dw,wd,w,p,ii,field,direction);
        //dw[fencode3_db(p,ii,field)]=-ddotcurrentmom(dw,wd,w,p,ii,field,direction);

  return ( status);
}

__device__ __host__
int dbderivsourceb (real *dw, real *wd, real *w, struct params *p,int *ii, int field, int direction) {

  int status=0;
        dw[fencode3_db(p,ii,field)]=dw[fencode3_db(p,ii,field)]+dbsourceb(dw,wd,w,p,ii,field,direction);

  return ( status);
}

__device__ __host__
int dbderivsourceenergy (real *dw, real *wd, real *w, struct params *p,int *ii) {

  int status=0;
  int field=energy;
        dw[fencode3_db(p,ii,field)]=dw[fencode3_db(p,ii,field)]+dbsourceenergy(dw,wd,w,p,ii);

  return ( status);
}

//rho, mom1, mom2, mom3, energy, b1, b2, b3
__device__ __host__
void dbderivsource (real *dw, real *wd, real *w, struct params *p,int *ii, int field) {

  //int status=0;
  switch(field)
  {
     case rho:
      dbderivsourcerho(dw,wd,w,p,ii);
     break;
     case mom1:
      dbderivsourcemom(dw,wd,w,p,ii,field,0);
     break;
     case mom2:
      dbderivsourcemom(dw,wd,w,p,ii,field,1);
     break;
   #ifdef USE_SAC_3D
     case mom3:
      dbderivsourcemom(dw,wd,w,p,ii,field,2);
     break;
   #endif
     case energy:
       dbderivsourceenergy(dw,wd,w,p,ii);
     break;
     case b1:
      dbderivsourceb(dw,wd,w,p,ii,field,0);
     break;
     case b2:
      dbderivsourceb(dw,wd,w,p,ii,field,1);
     break;
   #ifdef USE_SAC_3D
     case b3:
      dbderivsourceb(dw,wd,w,p,ii,field,2);
     break;
   #endif
  }
  //return ( status);
}


__global__ void divb_parallel(struct params *p, real *w, real *wmod, 
    real *dwn1, real *wd, int order,int ordero, real dt)
{
  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j,k;

  int ni=p->n[0];
  int nj=p->n[1];

   int ip,jp,ipg,jpg;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int nk,kp,kpg;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   //int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni/((p->npgp[1])*(p->npgp[0])));
   jp=(iindex-(kp*(nj*ni/((p->npgp[1])*(p->npgp[0])))))/(ni/(p->npgp[0]));
   ip=iindex-(kp*nj*ni/((p->npgp[1])*(p->npgp[0])))-(jp*(ni/(p->npgp[0])));
#endif
 #if defined USE_SAC || defined ADIABHYDRO
    jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));
#endif  


int shift=order*NVAR*dimp;

   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   #ifdef USE_SAC_3D
     for(kpg=0;kpg<(p->npgp[2]);kpg++)
   #endif
   {

     iia[0]=ip*(p->npgp[0])+ipg;
     iia[1]=jp*(p->npgp[1])+jpg;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp*(p->npgp[2])+kpg;
           k=iia[2];
     #endif

     #ifdef USE_SAC_3D
      if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
     #else
       if(i<((p->n[0])) && j<((p->n[1])))
     #endif
     for(int f=rho; f<=b2; f++)
                dwn1[fencode3_db(p,iia,f)]=0;
   }
 __syncthreads();

   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   #ifdef USE_SAC_3D
     for(kpg=0;kpg<(p->npgp[2]);kpg++)
   #endif
   {

     iia[0]=ip*(p->npgp[0])+ipg;
     iia[1]=jp*(p->npgp[1])+jpg;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp*(p->npgp[2])+kpg;
           k=iia[2];
     #endif

     #ifdef USE_SAC_3D
      if(i>2 && j>2 && k>2 && i<(ni-2) && j<(nj-2) && k<(nk-2))
     #else
       if(i>2 && j>2 && i<(ni-2) && j<(nj-2))
     #endif
  //if(i>2 && j>2 && i<(ni-2) && j<(nj-2))
	{
           if(p->divbfix)
           {   

               wd[fencode3_db(p,iia,divb)]=grad3d_db(wmod+order*NVAR*dimp,p,iia,b1,0)+grad3d_db(wmod+order*NVAR*dimp,p,iia,b2,1);
               #ifdef USE_SAC
		wd[fencode3_db(p,iia,divb)]+=grad3d_db(wmod+order*NVAR*dimp,p,iia,b1b,0)+grad3d_db(wmod+order*NVAR*dimp,p,iia,b2b,1);
                #endif
               #ifdef USE_SAC_3D
		wd[fencode3_db(p,iia,divb)]+=grad3d_db(wmod+order*NVAR*dimp,p,iia,b3,0)+grad3d_db(wmod+order*NVAR*dimp,p,iia,b3b,1);
                #endif
               for(int f=rho; f<=b2; f++) 
               {              
                  dbderivsource(dwn1,wd,wmod+order*NVAR*dimp,p,iia,f);
 
               }
            }

	}
}
 __syncthreads();


   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   #ifdef USE_SAC_3D
     for(kpg=0;kpg<(p->npgp[2]);kpg++)
   #endif
   {

     iia[0]=ip*(p->npgp[0])+ipg;
     iia[1]=jp*(p->npgp[1])+jpg;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp*(p->npgp[2])+kpg;
           k=iia[2];
     #endif

     #ifdef USE_SAC_3D
      if(i>1 && j >1 && k>1 && i<(ni-2) && j<(nj-2) && k<(nk-2))
     #else
       if(i>1 && j >1 && i<(ni-2) && j<(nj-2))
     #endif
   // if(i>1 && j >1 && i<(ni-2) && j<(nj-2))
                         {
                         if(p->divbfix)
                          { 
                             for(int f=rho; f<=b2; f++) 
                             //                                                  - sign here same as vac maybe a +
                              wmod[fencode3_db(p,iia,f)+(ordero*NVAR*dimp)]=wmod[fencode3_db(p,iia,f)+(ordero*NVAR*dimp)]-dt*dwn1[fencode3_db(p,iia,f)]; 
                          }

                         }
              //  }	
}
  __syncthreads();



  
}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_db(char *label)
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

int cudivb(struct params **p, struct params **d_p, real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd,  int order,int ordero, real dt)
{
    int status=0;
    dim3 dimBlock(dimblock, 1);
    
    dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (((*p)->n[0])*((*p)->n[1])+numThreadsPerBlock-1) / numThreadsPerBlock;


    divb_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w, *d_wmod, *d_dwn1,  *d_wd, order,ordero,dt);
	    //printf("called update\n"); 
    cudaThreadSynchronize();


 return status;


}



