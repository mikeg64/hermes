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
#include "../include/gradops_hdm1.cuh"


__global__ void hyperdifmomsource3_parallel(struct params *p,  real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real *wtemp, int field, int dim, int ii, int ii0, real dt)
{
  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int ii1;
  real fip,fim1,tmpc;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  //real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
  real rdx;
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
  rdx=(((wd[encode3_hdm1(p,i,j,k,delx1)])*(dim==0))+(wd[encode3_hdm1(p,i,j,k,delx2)])*(dim==1)+(wd[encode3_hdm1(p,i,j,k,delx3)])*(dim==2));
#else
  rdx=(((wd[encode3_hdm1(p,i,j,k,delx1)])*(dim==0))+  (wd[encode3_hdm1(p,i,j,k,delx2)])*(dim==1)  );
#endif


     #ifdef USE_SAC_3D
      if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
     #else
       if(i<((p->n[0])) && j<((p->n[1])))
     #endif

  //if(i<((p->n[0])) && j<((p->n[1])))
	{		               

;//dwn1[fencode3_hdm1(p,iia,energy)]=(wtemp[fencode3_hdm1(p,iia,tmp6)]*(wd[fencode3_hdm1(p,iia,hdnur)]+wd[fencode3_hdm1(p,iia,nushk1+dim)])*wtemp[fencode3_hdm1(p,iia,tmp8)]-wtemp[fencode3_hdm1(p,iia,tmp5)]*(wd[fencode3_hdm1(p,iia,hdnul)]+wd[fencode3_hdm1(p,iia,nushk1+dim)])*wtemp[fencode3_hdm1(p,iia,tmp7)])/(rdx)/2;

dwn1[fencode3_hdm1(p,iia,mom1+ii0)]=(wtemp[fencode3_hdm1(p,iia,tmp3)]*(wd[fencode3_hdm1(p,iia,hdnur)]+wd[fencode3_hdm1(p,iia,nushk1+dim)])*wtemp[fencode3_hdm1(p,iia,tmp8)]-wtemp[fencode3_hdm1(p,iia,tmp2)]*(wd[fencode3_hdm1(p,iia,hdnul)]+wd[fencode3_hdm1(p,iia,nushk1+dim)])*wtemp[fencode3_hdm1(p,iia,tmp7)])/(rdx)/2;

                              wmod[fencode3_hdm1(p,iia,mom1+ii0)+(ordero*NVAR*dimp)]=wmod[fencode3_hdm1(p,iia,mom1+ii0)+(ordero*NVAR*dimp)]+dt*dwn1[fencode3_hdm1(p,iia,mom1+ii0)];
                             //wmod[fencode3_hdm1(p,iia,mom1+ii0)+(ordero*NVAR*dimp)]=wmod[fencode3_hdm1(p,iia,mom1+ii0)+(ordero*NVAR*dimp)]+dt*((wtemp[fencode3_hdm1(p,iia,tmp3)]*(wd[fencode3_hdm1(p,iia,hdnur)]+wd[fencode3_hdm1(p,iia,nushk1+dim)])*wtemp[fencode3_hdm1(p,iia,tmp8)]-wtemp[fencode3_hdm1(p,iia,tmp2)]*(wd[fencode3_hdm1(p,iia,hdnul)]+wd[fencode3_hdm1(p,iia,nushk1+dim)])*wtemp[fencode3_hdm1(p,iia,tmp7)])/(rdx)/2);

   //del=wmod[fencode3_hdm1(p,iia,energy)+(ordero*NVAR*dimp)]+dt*dwn1[fencode3_hdm1(p,iia,energy)]; 
   //if(del<0.011 && del>0.009)
           //  wmod[fencode3_hdm1(p,iia,energy)+(ordero*NVAR*dimp)]=del; 
                            ;//wmod[fencode3_hdm1(p,iia,energy)+(ordero*NVAR*dimp)]=wmod[fencode3_hdm1(p,iia,energy)+(ordero*NVAR*dimp)]+dt*dwn1[fencode3_hdm1(p,iia,energy)]; 
                               
                             // wmod[fencode3_hdm1(p,iia,mom1+ii0)+(ordero*NVAR*dimp)]=wmod[fencode3_hdm1(p,iia,mom1+ii0)+(ordero*NVAR*dimp)]+dt*((wtemp[fencode3_hdm1(p,iia,tmp3)]*wd[fencode3_hdm1(p,iia,hdnur)]*wtemp[fencode3_hdm1(p,iia,tmp8)]-wtemp[fencode3_hdm1(p,iia,tmp2)]*wd[fencode3_hdm1(p,iia,hdnul)]*wtemp[fencode3_hdm1(p,iia,tmp7)])/(rdx)/2); 

                            // wmod[fencode3_hdm1(p,iia,energy)+(ordero*NVAR*dimp)]=wmod[fencode3_hdm1(p,iia,energy)+(ordero*NVAR*dimp)]+dt*((wtemp[fencode3_hdm1(p,iia,tmp6)]*wd[fencode3_hdm1(p,iia,hdnur)]*wtemp[fencode3_hdm1(p,iia,tmp8)]-wtemp[fencode3_hdm1(p,iia,tmp5)]*wd[fencode3_hdm1(p,iia,hdnul)]*wtemp[fencode3_hdm1(p,iia,tmp7)])/(rdx)/2); 


   }

 //__syncthreads();






   



  
}



__global__ void hyperdifmomsource3a_parallel(struct params *p,  real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real *wtemp, int field, int dim, int ii, int ii0, real dt)
{
  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int ii1;
  real fip,fim1,tmpc;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  //real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
  real rdx;
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
  rdx=(((wd[encode3_hdm1(p,i,j,k,delx1)])*(dim==0))+(wd[encode3_hdm1(p,i,j,k,delx2)])*(dim==1)+(wd[encode3_hdm1(p,i,j,k,delx3)])*(dim==2));
#else
  rdx=(((wd[encode3_hdm1(p,i,j,k,delx1)])*(dim==0))+  (wd[encode3_hdm1(p,i,j,k,delx2)])*(dim==1)  );
#endif

     #ifdef USE_SAC_3D
      if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
     #else
       if(i<((p->n[0])) && j<((p->n[1])))
     #endif

  //if(i<((p->n[0])) && j<((p->n[1])))
	{		               

dwn1[fencode3_hdm1(p,iia,energy)]=(wtemp[fencode3_hdm1(p,iia,tmp6)]*(wd[fencode3_hdm1(p,iia,hdnur)]+wd[fencode3_hdm1(p,iia,nushk1+dim)])*wtemp[fencode3_hdm1(p,iia,tmp8)]-wtemp[fencode3_hdm1(p,iia,tmp5)]*(wd[fencode3_hdm1(p,iia,hdnul)]+wd[fencode3_hdm1(p,iia,nushk1+dim)])*wtemp[fencode3_hdm1(p,iia,tmp7)])/(rdx)/2;

;//dwn1[fencode3_hdm1(p,iia,mom1+ii0)]=(wtemp[fencode3_hdm1(p,iia,tmp3)]*(wd[fencode3_hdm1(p,iia,hdnur)]+wd[fencode3_hdm1(p,iia,nushk1+dim)])*wtemp[fencode3_hdm1(p,iia,tmp8)]-wtemp[fencode3_hdm1(p,iia,tmp2)]*(wd[fencode3_hdm1(p,iia,hdnul)]+wd[fencode3_hdm1(p,iia,nushk1+dim)])*wtemp[fencode3_hdm1(p,iia,tmp7)])/(rdx)/2;

                             ;// wmod[fencode3_hdm1(p,iia,mom1+ii0)+(ordero*NVAR*dimp)]=wmod[fencode3_hdm1(p,iia,mom1+ii0)+(ordero*NVAR*dimp)]+dt*dwn1[fencode3_hdm1(p,iia,mom1+ii0)];

   //del=wmod[fencode3_hdm1(p,iia,energy)+(ordero*NVAR*dimp)]+dt*dwn1[fencode3_hdm1(p,iia,energy)]; 
   //if(del<0.011 && del>0.009)
           //  wmod[fencode3_hdm1(p,iia,energy)+(ordero*NVAR*dimp)]=del; 
                            wmod[fencode3_hdm1(p,iia,energy)+(ordero*NVAR*dimp)]=wmod[fencode3_hdm1(p,iia,energy)+(ordero*NVAR*dimp)]+dt*dwn1[fencode3_hdm1(p,iia,energy)]; 
                              // wmod[fencode3_hdm1(p,iia,energy)+(ordero*NVAR*dimp)]=wmod[fencode3_hdm1(p,iia,energy)+(ordero*NVAR*dimp)]+dt*((wtemp[fencode3_hdm1(p,iia,tmp6)]*(wd[fencode3_hdm1(p,iia,hdnur)]+wd[fencode3_hdm1(p,iia,nushk1+dim)])*wtemp[fencode3_hdm1(p,iia,tmp8)]-wtemp[fencode3_hdm1(p,iia,tmp5)]*(wd[fencode3_hdm1(p,iia,hdnul)]+wd[fencode3_hdm1(p,iia,nushk1+dim)])*wtemp[fencode3_hdm1(p,iia,tmp7)])/(rdx)/2); 
                             // wmod[fencode3_hdm1(p,iia,mom1+ii0)+(ordero*NVAR*dimp)]=wmod[fencode3_hdm1(p,iia,mom1+ii0)+(ordero*NVAR*dimp)]+dt*((wtemp[fencode3_hdm1(p,iia,tmp3)]*wd[fencode3_hdm1(p,iia,hdnur)]*wtemp[fencode3_hdm1(p,iia,tmp8)]-wtemp[fencode3_hdm1(p,iia,tmp2)]*wd[fencode3_hdm1(p,iia,hdnul)]*wtemp[fencode3_hdm1(p,iia,tmp7)])/(rdx)/2); 

                            // wmod[fencode3_hdm1(p,iia,energy)+(ordero*NVAR*dimp)]=wmod[fencode3_hdm1(p,iia,energy)+(ordero*NVAR*dimp)]+dt*((wtemp[fencode3_hdm1(p,iia,tmp6)]*wd[fencode3_hdm1(p,iia,hdnur)]*wtemp[fencode3_hdm1(p,iia,tmp8)]-wtemp[fencode3_hdm1(p,iia,tmp5)]*wd[fencode3_hdm1(p,iia,hdnul)]*wtemp[fencode3_hdm1(p,iia,tmp7)])/(rdx)/2); 


   }

 //__syncthreads();






   



  
}


__global__ void hyperdifmomsource2_parallel(struct params *p,  real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real *wtemp, int field, int dim, int ii, int ii0, real dt)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int ii1;
  real fip,fim1,tmpc;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  //real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
  real rdx;

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

/*#ifdef USE_SAC_3D
  rdx=(((p->dx[0])*(dim==0))+(p->dx[1])*(dim==1)+(p->dx[2])*(dim==2));
#else
  rdx=(((p->dx[0])*(dim==0))+  (p->dx[1])*(dim==1)  );
#endif*/

 

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
  {

     wtemp[fencode3_hdm1(p,iia,tmp8)]=grad1r3n_hdm1(wtemp,wd,p,iia,tmp4,dim);
     wtemp[fencode3_hdm1(p,iia,tmp7)]=grad1l3n_hdm1(wtemp,wd,p,iia,tmp4,dim);

   }


//__syncthreads();  //can remove?



  
}




__global__ void hyperdifmomsource1_parallel(struct params *p,  real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real *wtemp, int field, int dim, int ii, int ii0, real dt)
{
  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int ii1;
  real fip,fim1,tmpc;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  //real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
  real rdx;
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

/*#ifdef USE_SAC_3D
  rdx=(((p->dx[0])*(dim==0))+(p->dx[1])*(dim==1)+(p->dx[2])*(dim==2));
#else
  rdx=(((p->dx[0])*(dim==0))+  (p->dx[1])*(dim==1)  );
#endif*/



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
  //init rhol and rhor
  //if(i<((p->n[0])) && j<((p->n[1])))
  {
    for(int f=tmp1; f<=tmp8; f++)	
        wtemp[fencode3_hdm1(p,iia,f)]=0.0;


dwn1[fencode3_hdm1(p,iia,energy)]=0.0;
dwn1[fencode3_hdm1(p,iia,mom1+ii0)]=0.0;
   }



 //__syncthreads();

//tmp2  rhor
//tmp3  rhol
//tmp1  mom+field/rho

//tmp4  rhoc


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
//if( i<((p->n[0])) && j<((p->n[1])))
  {
#ifdef ADIABHYDRO
;
#else
    wtemp[fencode3_hdm1(p,iia,tmp4)]=wmod[(shift)+fencode3_hdm1(p,iia,mom1+field)]/(wmod[(shift)+fencode3_hdm1(p,iia,rho)]+wmod[(shift)+fencode3_hdm1(p,iia,rhob)]);
#endif

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
       if(i>0 && j >0 && k>0 && i<((p->n[0])-1) && j<((p->n[1])-1)  && k<((p->n[2])-1))
     #else
       if(i>0 && j >0 && i<((p->n[0])-1) && j<((p->n[1])-1))
     #endif
//if(i>0 && j >0 && i<((p->n[0])-1) && j<((p->n[1])-1))
  {

     #ifdef USE_SAC_3D
       wtemp[fencode3_hdm1(p,iia,tmp2)]=(wmod[(shift)+fencode3_hdm1(p,iia,rho)]+wmod[(shift)+fencode3_hdm1(p,iia,rhob)]+wmod[(shift)+encode3_hdm1(p,i-(dim==0),j-(dim==1),k-(dim==2),rho)]+wmod[(shift)+encode3_hdm1(p,i-(dim==0),j-(dim==1),k-(dim==2),rhob)])/2;
       wtemp[fencode3_hdm1(p,iia,tmp3)]=(wmod[(shift)+fencode3_hdm1(p,iia,rho)]+wmod[(shift)+fencode3_hdm1(p,iia,rhob)]+wmod[(shift)+encode3_hdm1(p,i+(dim==0),j+(dim==1),k+(dim==2),rho)]+wmod[(shift)+encode3_hdm1(p,i+(dim==0),j+(dim==1),k+(dim==2),rhob)])/2;
     #endif

     #ifdef USE_SAC
       wtemp[fencode3_hdm1(p,iia,tmp2)]=(wmod[(shift)+fencode3_hdm1(p,iia,rho)]+wmod[(shift)+fencode3_hdm1(p,iia,rhob)]+wmod[(shift)+fencode_hdm1(p,i-(dim==0),j-(dim==1),rho)]+wmod[(shift)+fencode_hdm1(p,i-(dim==0),j-(dim==1),rhob)])/2;
       wtemp[fencode3_hdm1(p,iia,tmp3)]=(wmod[(shift)+fencode3_hdm1(p,iia,rho)]+wmod[(shift)+fencode3_hdm1(p,iia,rhob)]+wmod[(shift)+fencode_hdm1(p,i+(dim==0),j+(dim==1),rho)]+wmod[(shift)+fencode_hdm1(p,i+(dim==0),j+(dim==1),rhob)])/2;
     #endif


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
       if(i>0 && j >0 && k>0 && i<((p->n[0])-1) && j<((p->n[1])-1)  && k<((p->n[2])-1))
     #else
       if(i>0 && j >0 && i<((p->n[0])-1) && j<((p->n[1])-1))
     #endif
//if(i>0 && j >0 && i<((p->n[0])-1) && j<((p->n[1])-1))
  {
     #ifdef USE_SAC_3D
     wtemp[fencode3_hdm1(p,iia,tmp5)]=(wmod[(shift)+fencode3_hdm1(p,iia,mom1+ii0)]+wmod[(shift)+encode3_hdm1(p,i-(dim==0),j-(dim==1),k-(dim==2),mom1+ii0)])/2;
     wtemp[fencode3_hdm1(p,iia,tmp6)]=(wmod[(shift)+fencode3_hdm1(p,iia,mom1+ii0)]+wmod[(shift)+encode3_hdm1(p,i+(dim==0),j+(dim==1),k+(dim==2),mom1+ii0)])/2;
     #else
     wtemp[fencode3_hdm1(p,iia,tmp5)]=(wmod[(shift)+fencode3_hdm1(p,iia,mom1+ii0)]+wmod[(shift)+fencode_hdm1(p,i-(dim==0),j-(dim==1),mom1+ii0)])/2;
     wtemp[fencode3_hdm1(p,iia,tmp6)]=(wmod[(shift)+fencode3_hdm1(p,iia,mom1+ii0)]+wmod[(shift)+fencode_hdm1(p,i+(dim==0),j+(dim==1),mom1+ii0)])/2;
     #endif
   }


//__syncthreads();



  
}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_hdm1(char *label)
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





int cuhyperdifmomsource1(struct params **p,  struct params **d_p,   real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero, real **d_wtemp, int field, int dim, int ii, int ii0, real dt)
{
  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 
   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;

     hyperdifmomsource1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim,ii,ii0,dt);
     cudaThreadSynchronize();

     hyperdifmomsource2_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim,ii,ii0,dt);
     cudaThreadSynchronize();

     hyperdifmomsource3_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim,ii,ii0,dt);
     cudaThreadSynchronize();
     hyperdifmomsource3a_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim,ii,ii0,dt);
     cudaThreadSynchronize();

}







