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
#include "../include/gradops_hdb1.cuh"






__global__ void hyperdifbsource4_parallel(struct params *p,  real *wmod, 
    real *dwn1, real *wd, int order,int ordero, real *wtemp, int field, int dim, int jj, int ii0,int mm,real sb,real dt)
{

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int m,ii1;
  real fip,fim1,tmpc;
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
  real rdx=(((wd[encode3_hdb1(p,i,j,k,delx1)])*(dim==0))+(wd[encode3_hdb1(p,i,j,k,delx2)])*(dim==1)+(wd[encode3_hdb1(p,i,j,k,delx3)])*(dim==2));
#else
  real rdx=(((wd[encode3_hdb1(p,i,j,k,delx1)])*(dim==0))+  (wd[encode3_hdb1(p,i,j,k,delx2)])*(dim==1)  );
#endif


     #ifdef USE_SAC_3D
     if(i<((p->n[0])) && j<((p->n[1])) && k<((p->n[2])))
     #else
    if(i<((p->n[0])) && j<((p->n[1])))
     #endif

                         {
                              //                                                                                  - sign here same as vac maybe a +
                              wmod[fencode3_hdb1(p,iia,b1+ii0)+(ordero*NVAR*dimp)]=wmod[fencode3_hdb1(p,iia,b1+ii0)+(ordero*NVAR*dimp)]+dt*dwn1[fencode3_hdb1(p,iia,b1+ii0)]; 
                             wmod[fencode3_hdb1(p,iia,energy)+(ordero*NVAR*dimp)]=wmod[fencode3_hdb1(p,iia,energy)+(ordero*NVAR*dimp)]+dt*dwn1[fencode3_hdb1(p,iia,energy)]; 

                         }
              //  }	

  //__syncthreads();  
}




__global__ void hyperdifbsource3_parallel(struct params *p, real *wmod, 
    real *dwn1, real *wd, int order,int ordero, real *wtemp, int field, int dim, int jj, int ii0,int mm,real sb,real dt)
{
 
  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int m,ii1;
  real fip,fim1,tmpc;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];

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
  real rdx=(((wd[encode3_hdb1(p,i,j,k,delx1)])*(dim==0))+(wd[encode3_hdb1(p,i,j,k,delx2)])*(dim==1)+(wd[encode3_hdb1(p,i,j,k,delx3)])*(dim==2));
#else
  real rdx=(((wd[encode3_hdb1(p,i,j,k,delx1)])*(dim==0))+  (wd[encode3_hdb1(p,i,j,k,delx2)])*(dim==1)  );
#endif



     #ifdef USE_SAC_3D
     if(i<((p->n[0])) && j<((p->n[1])) && k<((p->n[2])))
     #else
    if(i<((p->n[0])) && j<((p->n[1])))
     #endif
	{		               



dwn1[fencode3_hdb1(p,iia,b1+ii0)]=sb*(wtemp[fencode3_hdb1(p,iia,tmp5)]*wd[fencode3_hdb1(p,iia,hdnur)]-wtemp[fencode3_hdb1(p,iia,tmp4)]*wd[fencode3_hdb1(p,iia,hdnul)])/rdx;

dwn1[fencode3_hdb1(p,iia,energy)]=sb*(wtemp[fencode3_hdb1(p,iia,tmp3)]*wtemp[fencode3_hdb1(p,iia,tmp5)]*wd[fencode3_hdb1(p,iia,hdnur)]-wtemp[fencode3_hdb1(p,iia,tmp2)]*wtemp[fencode3_hdb1(p,iia,tmp4)]*wd[fencode3_hdb1(p,iia,hdnul)])/rdx;

//if(i==0 && j==139)
//           printf("b1 e %d %10.20g %10.20g\n",ii0,dwn1[fencode3_hdb1(p,iia,b1+ii0)],dwn1[fencode3_hdb1(p,iia,energy)]);



//    if(i==127 && j==2)
//           printf("tmpL R %10.20g %10.20g\n",wtemp[fencode3_hdb1(p,iia,tmp4)]*wd[fencode3_hdb1(p,iia,hdnul)],wtemp[fencode3_hdb1(p,iia,tmp5)]*wd[fencode3_hdb1(p,iia,hdnur)]);
//    if(i==127 && j==2)
//           printf("nuL R %10.20g %10.20g\n",wd[fencode3_hdb1(p,iia,hdnul)],wd[fencode3_hdb1(p,iia,hdnur)]);

   }

 //__syncthreads();


 
}






__global__ void hyperdifbsource2_parallel(struct params *p,  real *wmod, 
    real *dwn1, real *wd, int order,int ordero, real *wtemp, int field, int dim, int jj, int ii0,int mm,real sb)
{

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int m,ii1;
  real fip,fim1,tmpc;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];

  real dt=p->dt;
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
     //if(i<((p->n[0])) && j<((p->n[1])) && k<((p->n[2])))
     if(i>0 && j >0 && k>0 && i<((p->n[0])-1) && j<((p->n[1])-1)  && k < ((p->n[2])-1))
     #else
    //if(i<((p->n[0])) && j<((p->n[1])))
    if(i>0 && j >0 && i<((p->n[0])-1) && j<((p->n[1])-1))
     #endif
 // if(i>0 && j >0 && i<((p->n[0])-1) && j<((p->n[1])-1))
	{		               

//wtemp[fencode3_hdb1(p,iia,tmp4)]=grad1l3n_hdb1(wtemp,wd,p,iia,tmp1,dim);
//wtemp[fencode3_hdb1(p,iia,tmp5)]=grad1r3n_hdb1(wtemp,wd,p,iia,tmp1,dim);


wtemp[fencode3_hdb1(p,iia,tmp4)]=grad1l3n_hdb1(wtemp,wd,p,iia,tmp1,dim);
wtemp[fencode3_hdb1(p,iia,tmp5)]=grad1r3n_hdb1(wtemp,wd,p,iia,tmp1,dim);

   }

 //__syncthreads();   

  //   if(i==127 && j==2)
  //         printf("L R %d %10.20g %10.20g\n",dim,wtemp[fencode3_hdb1(p,iia,tmp4)],wtemp[fencode3_hdb1(p,iia,tmp5)]);


   



}



__global__ void hyperdifbsource1_parallel(struct params *p,  real *wmod, 
    real *dwn1, real *wd, int order,int ordero, real *wtemp, int field, int dim, int jj, int ii0,int mm,real sb)
{
  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int m,ii1;
  real fip,fim1,tmpc;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];

  real dt=p->dt;
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

  //init rhol and rhor

  {
    for(int f=tmp1; f<=tmp8; f++)	
        wtemp[fencode3_hdb1(p,iia,f)]=0.0;

   dwn1[fencode3_hdb1(p,iia,energy)]=0.0;
   dwn1[fencode3_hdb1(p,iia,b1+ii0)]=0.0;
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
     if(i>0 && j >0 && k>0 && i<((p->n[0])-1) && j<((p->n[1])-1) && k<((p->n[2])-1))
     #else
    if(i>0 && j >0 && i<((p->n[0])-1) && j<((p->n[1])-1))
     #endif

  //if(i>0 && j >0 && i<((p->n[0])-1) && j<((p->n[1])-1))
  {

     #ifdef USE_SAC_3D
       wtemp[fencode3_hdb1(p,iia,tmp2)]=(wmod[shift+fencode3_hdb1(p,iia,b1+jj)]+wmod[shift+encode3_hdb1(p,i-(dim==0),j-(dim==1),k-(dim==2),b1+jj)])/2;
       wtemp[fencode3_hdb1(p,iia,tmp3)]=(wmod[shift+fencode3_hdb1(p,iia,b1+jj)]+wmod[shift+encode3_hdb1(p,i+(dim==0),j+(dim==1),k+(dim==2),b1+jj)])/2;
     #else
       wtemp[fencode3_hdb1(p,iia,tmp2)]=(wmod[shift+fencode3_hdb1(p,iia,b1+jj)]+wmod[shift+encode3_hdb1(p,i-(dim==0),j-(dim==1),k,b1+jj)])/2;
       wtemp[fencode3_hdb1(p,iia,tmp3)]=(wmod[shift+fencode3_hdb1(p,iia,b1+jj)]+wmod[shift+encode3_hdb1(p,i+(dim==0),j+(dim==1),k,b1+jj)])/2;
     #endif
     wtemp[fencode3_hdb1(p,iia,tmp1)]=wmod[shift+fencode3_hdb1(p,iia,b1+field)];

     //if(i==127 && j==2)
     //      printf("tmp1 %d %10.20g\n",field,wmod[shift+fencode3_hdb1(p,iia,b1+field)]);

   }


//__syncthreads();






   

}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_hdb1(char *label)
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





int cuhyperdifbsource1(struct params **p,  struct params **d_p,   real **d_wmod, real **d_dwn1, real **d_wd, int order,int ordero, real **d_wtemp, int field, int dim, int jj, int ii0,int mm,real sb, real dt)
{
  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 
   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;


     hyperdifbsource1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim,jj,ii0,mm,sb); 
     cudaThreadSynchronize();
     hyperdifbsource2_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim,jj,ii0,mm,sb); 
     cudaThreadSynchronize();
     hyperdifbsource3_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim,jj,ii0,mm,sb,dt); 
     cudaThreadSynchronize();
     hyperdifbsource4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim,jj,ii0,mm,sb,dt); 
     cudaThreadSynchronize();

}







