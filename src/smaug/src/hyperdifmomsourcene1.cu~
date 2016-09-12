#include "../include/cudapars.h"
#include "../include/paramssteeringtest1.h"
#include "../include/iobparams.h"
/////////////////////////////////////
// standard imports
/////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include "../include/smaugcukernels.h"

/////////////////////////////////////
// kernel function (CUDA device)
/////////////////////////////////////
#include "../include/gradops_hdmne1.cuh"

__global__ void hyperdifmomsourcene6_parallel(struct params *p,  real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real *wtemp, int field, int dim, int ii, int ii0, real dt)
{

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int ii1;
  real fip,fim1;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  //real dt=p->dt;
   int ip,jp;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   //int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif  

int shift=order*NVAR*dimp;





     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];
     #endif

     #ifdef USE_SAC_3D
      if(i<((p->n[0])) && j<((p->n[1])) && k<((p->n[2])))
     #else
      if(i<((p->n[0])) && j<((p->n[1])))
     #endif

                        //if(i<((p->n[0])) && j<((p->n[1])))
                         {

                             wmod[fencode3_hdmne1(p,iia,mom1+ii0)+(ordero*NVAR*dimp)]=wmod[fencode3_hdmne1(p,iia,mom1+ii0)+(ordero*NVAR*dimp)]+dt*dwn1[fencode3_hdmne1(p,iia,mom1+ii0)];

 
                           //  wmod[fencode3_hdmne1(p,iia,energy)+(ordero*NVAR*dimp)]=wmod[fencode3_hdmne1(p,iia,energy)+(ordero*NVAR*dimp)]+dt*dwn1[fencode3_hdmne1(p,iia,energy)];
    //if(i==127 && j==252)
    //  p->test=dt*dwn1[fencode3_hdmne1(p,iia,mom1+ii0)];

                         }
              //  }	

  //__syncthreads();


  



}


__global__ void hyperdifmomsourcene6a_parallel(struct params *p,  real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real *wtemp, int field, int dim, int ii, int ii0, real dt)
{

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int ii1;
  real fip,fim1;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  //real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
   int ip,jp;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   //int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif  

int shift=order*NVAR*dimp;





     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];
     #endif

     #ifdef USE_SAC_3D
      if(i<((p->n[0])) && j<((p->n[1])) && k<((p->n[2])))
     #else
      if(i<((p->n[0])) && j<((p->n[1])))
     #endif

                        //if(i<((p->n[0])) && j<((p->n[1])))
                         {

                           //  wmod[fencode3_hdmne1(p,iia,mom1+ii0)+(ordero*NVAR*dimp)]=wmod[fencode3_hdmne1(p,iia,mom1+ii0)+(ordero*NVAR*dimp)]+dt*dwn1[fencode3_hdmne1(p,iia,mom1+ii0)];

 
                             wmod[fencode3_hdmne1(p,iia,energy)+(ordero*NVAR*dimp)]=wmod[fencode3_hdmne1(p,iia,energy)+(ordero*NVAR*dimp)]+dt*dwn1[fencode3_hdmne1(p,iia,energy)];

   // if(i==127 && j==252)
   //   p->test=dt*dwn1[fencode3_hdmne1(p,iia,energy)];

   // if(i==127 && j==252)
   //   p->test=wmod[fencode3_hdmne1(p,iia,energy)+(ordero*NVAR*dimp)];

                         }
              //  }	

  //__syncthreads();


  



}






__global__ void hyperdifmomsourcene5_parallel(struct params *p,  real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real *wtemp, int field, int dim, int ii, int ii0, real dt)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int ii1;
  real fip,fim1;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  //real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
   int ip,jp;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   //int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif  

int shift=order*NVAR*dimp;




     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];
     #endif

     #ifdef USE_SAC_3D
      if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
     #else
       if(i<((p->n[0])) && j<((p->n[1])))
     #endif

{
 
    dwn1[fencode3_hdmne1(p,iia,mom1+ii0)]=(grad13n_hdmne1(wtemp,wd,p,iia,tmp7,ii));


  }

 //__syncthreads();


  
  



}



__global__ void hyperdifmomsourcene5a_parallel(struct params *p,  real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real *wtemp, int field, int dim, int ii, int ii0, real dt)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int ii1;
  real fip,fim1;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  //real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
   int ip,jp;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   //int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif  

int shift=order*NVAR*dimp;




 


     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];
     #endif

     #ifdef USE_SAC_3D
      if(i>0 && j >0 && k>0 && i<((p->n[0])-1) && j<((p->n[1])-1) && k<((p->n[2])-1))
     #else
      if(i>0 && j >0 && i<((p->n[0])-1) && j<((p->n[1])-1))
     #endif

     dwn1[fencode3_hdmne1(p,iia,energy)]=(grad13n_hdmne1(wtemp,wd,p,iia,tmp8,ii));
     //dwn1[fencode3_hdmne1(p,iia,energy)]=-2.9e-4;
    //if(i==127 && j==252)
    //  p->test=(p->dt)*dwn1[fencode3_hdmne1(p,iia,energy)];



 //__syncthreads();


  



}

__global__ void hyperdifmomsourcene4_parallel(struct params *p,  real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real *wtemp, int field, int dim, int ii, int ii0, real dt)
{
  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int ii1;
  real fip,fim1;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  //real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
   int ip,jp;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   //int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif  

int shift=order*NVAR*dimp;




     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];
     #endif

     #ifdef USE_SAC_3D
      if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
     #else
       if(i<((p->n[0])) && j<((p->n[1])))
     #endif
  //if(i<((p->n[0])) && j<((p->n[1])))
	{		               
     wtemp[fencode3_hdmne1(p,iia,tmp7)]=wtemp[fencode3_hdmne1(p,iia,tmp1)]*wtemp[fencode3_hdmne1(p,iia,tmp6)];

     wtemp[fencode3_hdmne1(p,iia,tmp8)]=wtemp[fencode3_hdmne1(p,iia,tmp6)]*wmod[(shift)+fencode3_hdmne1(p,iia,mom1+ii0)];

  //  if(i==127 && j==252)
  //    p->test=wtemp[fencode3_hdmne1(p,iia,tmp7)];



   }

 //__syncthreads();




  



}


__global__ void hyperdifmomsourcene3_parallel(struct params *p,  real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real *wtemp, int field, int dim, int ii, int ii0, real dt)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int ii1;
  real fip,fim1;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  //real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];

   int ip,jp;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   //int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif  

int shift=order*NVAR*dimp;




     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];
     #endif

     #ifdef USE_SAC_3D
      if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
     #else
       if(i<((p->n[0])) && j<((p->n[1])))
     #endif
  //if(i<((p->n[0])) && j<((p->n[1])))
  {

     //wtemp[fencode3_hdmne1(p,iia,tmp6)]=wtemp[fencode3_hdmne1(p,iia,tmp5)]*((wd[fencode3_hdmne1(p,iia,hdnur)]+wd[fencode3_hdmne1(p,iia,hdnul)]+2.0*wd[fencode3_hdmne1(p,iia,nushk1+dim)]))/4.0;
     //wtemp[fencode3_hdmne1(p,iia,tmp6)]=1.0;
     // wtemp[fencode3_hdmne1(p,iia,tmp6)]=-1.0e-41*((wd[fencode3_hdmne1(p,iia,hdnur)]+wd[fencode3_hdmne1(p,iia,hdnul)]+2.0*wd[fencode3_hdmne1(p,iia,nushk1+dim)]))/4.0;
    wtemp[fencode3_hdmne1(p,iia,tmp6)]=wtemp[fencode3_hdmne1(p,iia,tmp5)]*((wd[fencode3_hdmne1(p,iia,hdnur)]+wd[fencode3_hdmne1(p,iia,hdnul)]+2.0*wd[fencode3_hdmne1(p,iia,nushk1+dim)]))/4.0;


 //wtemp[fencode3_hdmne1(p,iia,tmp6)]=wtemp[fencode3_hdmne1(p,iia,tmp5)]*((1.4e-4))/4.0;
  //  if(i==127 && j==252)
  //    p->test=((wd[fencode3_hdmne1(p,iia,hdnur)]+wd[fencode3_hdmne1(p,iia,hdnul)]+2.0*wd[fencode3_hdmne1(p,iia,nushk1+dim)]))/4.0;
;//   if(i==127 && j==252)
;//      p->test=((wtemp[fencode3_hdmne1(p,iia,tmp5)]));


   }

//__syncthreads();

/*if(iindex==0)
{

    for(iia[0]=0;iia[0]<((p->n[0]));iia[0]++)
      for(iia[1]=0;iia[1]<((p->n[1]));iia[1]++)
        //(p->test)+=wtemp[fencode3_hdmne1(p,iia,tmp5)];
        if((wtemp[fencode3_hdmne1(p,iia,tmp5)])>(p->test))
               p->test=wtemp[fencode3_hdmne1(p,iia,tmp5)];

//   p->test/=((p->n[0])*(p->n[1]));


}*/




  



}

__global__ void hyperdifmomsourcene2_parallel(struct params *p,  real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real *wtemp, int field, int dim, int ii, int ii0, real dt)
{

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int ii1;
  real fip,fim1;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  //real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];

   int ip,jp;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   //int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif  

int shift=order*NVAR*dimp;




     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];
     #endif

     #ifdef USE_SAC_3D
      if(i>0 && j >0 && k>0 && i<((p->n[0])-1) && j<((p->n[1])-1)  && k<((p->n[2])-1))
     #else
       if(i>0 && j >0 && i<((p->n[0])-1) && j<((p->n[1])-1))
     #endif
  //if(i>0 && j >0 && i<((p->n[0])-1) && j<((p->n[1])-1))
      wtemp[fencode3_hdmne1(p,iia,tmp5)]=(grad13n_hdmne1(wtemp,wd,p,iia,tmp4,dim));
        // if(i==127 && j==252)
        //       p->test=wtemp[fencode3_hdmne1(p,iia,tmp5)];



//__syncthreads();

}


__global__ void hyperdifmomsourcene1_parallel(struct params *p,  real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real *wtemp, int field, int dim, int ii, int ii0, real dt)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int ii1;
  real fip,fim1;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  //real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];



  

   int ip,jp;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   //int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif  

int shift=order*NVAR*dimp;



  //init rhol and rhor

     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];
     #endif

     #ifdef USE_SAC_3D
      if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
     #else
       if(i<((p->n[0])) && j<((p->n[1])))
     #endif
  //if(i<((p->n[0])) && j<((p->n[1])))
  {
    for(int f=tmp1; f<=tmp8; f++)	
        wtemp[fencode3_hdmne1(p,iia,f)]=0.0;

     dwn1[fencode3_hdmne1(p,iia,energy)]=0.0;
     dwn1[fencode3_hdmne1(p,iia,mom1+ii0)]=0.0;

   }

 //__syncthreads();



     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];
     #endif

     #ifdef USE_SAC_3D
      if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
     #else
       if(i<((p->n[0])) && j<((p->n[1])))
     #endif
  //if(i<((p->n[0])) && j<((p->n[1])))
  {

     #ifdef ADIABHYDRO
;
    #else
     wtemp[fencode3_hdmne1(p,iia,tmp1)]=wmod[(shift)+fencode3_hdmne1(p,iia,rho)]+wmod[(shift)+fencode3_hdmne1(p,iia,rhob)];

     wtemp[fencode3_hdmne1(p,iia,tmp4)]=wmod[(shift)+fencode3_hdmne1(p,iia,mom1+field)]/(wmod[(shift)+fencode3_hdmne1(p,iia,rho)]+wmod[(shift)+fencode3_hdmne1(p,iia,rhob)]);
    #endif



   }

//__syncthreads();


}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_hdmne1ne(char *label)
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





int cuhyperdifmomsourcene1(struct params **p, struct params **d_p, real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero, real **d_wtemp, int field, int dim, int ii, int ii0, real dt)
{
  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   //cudaSetDevice(selectedDevice);
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 
   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;


     hyperdifmomsourcene1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim,ii,ii0,dt);
     cudaThreadSynchronize();
     hyperdifmomsourcene2_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim,ii,ii0,dt);
     cudaThreadSynchronize();
     hyperdifmomsourcene3_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim,ii,ii0,dt);
     cudaThreadSynchronize();





     hyperdifmomsourcene4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim,ii,ii0,dt);
     cudaThreadSynchronize();
     hyperdifmomsourcene5_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim,ii,ii0,dt);
     cudaThreadSynchronize();
     hyperdifmomsourcene6_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim,ii,ii0,dt);
     cudaThreadSynchronize();

     hyperdifmomsourcene5a_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim,ii,ii0,dt);
     cudaThreadSynchronize();
     hyperdifmomsourcene6a_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim,ii,ii0,dt);
     cudaThreadSynchronize();
cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);
//printf("tmp2 %d %d %10.20g\n",ii,dim,(*p)->test);
//cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);
////printf("test %g\n",(*p)->test);

}







