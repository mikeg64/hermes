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
#include "../include/gradops_hde1.cuh"



__global__ void hyperdifesource4_parallel(struct params *p,  real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real *wtemp, int field, int dim, real dt)
{

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int ii1,ii0;
  real fip,fim1;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  //real dt=p->dt;
  real rdx;
  real dy=p->dx[1];
  real dx=p->dx[0];

   int ip,jp;
  int ii[NDIM];
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


real del;


     ii[0]=ip;
     ii[1]=jp;
     i=ii[0];
     j=ii[1];
     k=0;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
           k=ii[2];
     #endif

#ifdef USE_SAC_3D
  rdx=(((wd[encode3_hde1(p,i,j,k,delx1)])*(dim==0))+(wd[encode3_hde1(p,i,j,k,delx2)])*(dim==1)+(wd[encode3_hde1(p,i,j,k,delx3)])*(dim==2));
#else
  rdx=(((wd[encode3_hde1(p,i,j,k,delx1)])*(dim==0))+  (wd[encode3_hde1(p,i,j,k,delx2)])*(dim==1)  );
#endif
     #ifdef USE_SAC_3D
       if(i<((p->n[0])) && j<((p->n[1])) && k<((p->n[2])))
     #else
       if(i<((p->n[0])) && j<((p->n[1])))
     #endif
  {


//dwn1[fencode3_hde1(p,ii,field)]=( wtemp[fencode3_hde1(p,ii,hdnur)] *wtemp[fencode3_hde1(p,ii,tmp3)] - wtemp[fencode3_hde1(p,ii,hdnul)] *wtemp[fencode3_hde1(p,ii,tmp2)])/rdx;

   // wmod[fencode3_hde1(p,ii,field)+(ordero*NVAR*dimp)]=wmod[fencode3_hde1(p,ii,field)+(ordero*NVAR*dimp)]+dt*dwn1[fencode3_hde1(p,ii,field)]; 
   //del=wmod[fencode3_hde1(p,ii,field)+(ordero*NVAR*dimp)]+dt*dwn1[fencode3_hde1(p,ii,field)]; 
  // if(del<0.011 && del>0.009)
   //          wmod[fencode3_hde1(p,ii,field)+(ordero*NVAR*dimp)]=del;

   wmod[fencode3_hde1(p,ii,field)+(ordero*NVAR*dimp)]=wmod[fencode3_hde1(p,ii,field)+(ordero*NVAR*dimp)]+dt*( (wd[fencode3_hde1(p,ii,hdnur)]+wd[fencode3_hde1(p,ii,nushk1+dim)]) *wtemp[fencode3_hde1(p,ii,tmp3)] - (wd[fencode3_hde1(p,ii,hdnul)]+wd[fencode3_hde1(p,ii,nushk1+dim)]) *wtemp[fencode3_hde1(p,ii,tmp2)])/rdx;

  }

//__syncthreads();


/*if(iindex==0)
{
  p->hdmean=0.0;
  p->hdmax=0;

    for(ii[0]=0;ii[0]<((p->n[0]));ii[0]++)
      for(ii[1]=0;ii[1]<((p->n[1]));ii[1]++)
     #ifdef USE_SAC_3D
        for(ii[2]=0;ii[2]<((p->n[2]));ii[2]++)
     #endif
	{ 

             if((wtemp[encode3_hde1(p,ii[0],ii[1],0,tmp2)])>(p->hdmax))
                    p->hdmax=(wtemp[encode3_hde1(p,ii[0],ii[1],0,tmp2)]);
              p->hdmean=(p->hdmean)+wtemp[encode3_hde1(p,ii[0],ii[1],0,tmp2)];
	}
       p->hdmean=(p->hdmean)/(((p->n[0]))*((p->n[1])));

}
 //__syncthreads();*/


 
}


__global__ void hyperdifesource3_parallel(struct params *p,  real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real *wtemp, int field, int dim)
{

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int ii1,ii0;
  real fip,fim1;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real rdx;
  real dy=p->dx[1];
  real dx=p->dx[0];

   int ip,jp;
  int ii[NDIM];
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







     ii[0]=ip;
     ii[1]=jp;
     i=ii[0];
     j=ii[1];
     k=0;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
           k=ii[2];
     #endif


	#ifdef USE_SAC_3D
	  rdx=(((wd[encode3_hde1(p,i,j,k,delx1)])*(dim==0))+(wd[encode3_hde1(p,i,j,k,delx2)])*(dim==1)+(wd[encode3_hde1(p,i,j,k,delx3)])*(dim==2));
	#else
	  rdx=(((wd[encode3_hde1(p,i,j,k,delx1)])*(dim==0))+  (wd[encode3_hde1(p,i,j,k,delx2)])*(dim==1)  );
	#endif




     #ifdef USE_SAC_3D
       if(i<((p->n[0])) && j<((p->n[1])) && k<((p->n[2])))
     #else
       if(i<((p->n[0])) && j<((p->n[1])))
     #endif  
  {


 
dwn1[fencode3_hde1(p,ii,field)]=( (wd[fencode3_hde1(p,ii,hdnur)]+wd[fencode3_hde1(p,ii,nushk1+dim)]) *wtemp[fencode3_hde1(p,ii,tmp3)] - (wd[fencode3_hde1(p,ii,hdnul)]+wd[fencode3_hde1(p,ii,nushk1+dim)]) *wtemp[fencode3_hde1(p,ii,tmp2)])/rdx;
   


  }

//__syncthreads();



   
/*   for(ipg=0;ipg<(p->npgp[0]);ipg++)
   for(jpg=0;jpg<(p->npgp[1]);jpg++)
   #ifdef USE_SAC_3D
     for(kpg=0;kpg<(p->npgp[2]);kpg++)
   #endif
   {

     ii[0]=ip*(p->npgp[0])+ipg;
     ii[1]=jp*(p->npgp[1])+jpg;
     i=ii[0];
     j=ii[1];
     k=0;
     #ifdef USE_SAC_3D
	   ii[2]=kp*(p->npgp[2])+kpg;
           k=ii[2];
     #endif

     #ifdef USE_SAC_3D
       if(i<((p->n[0])) && j<((p->n[1])) && k<((p->n[2])))
     #else
       if(i<((p->n[0])) && j<((p->n[1])))
     #endif  
                         //if(i<((p->n[0])) && j<((p->n[1])))
                         {
                              //                                                                                  - sign here same as vac maybe a +
                              wmod[fencode3_hde1(p,ii,field)+(ordero*NVAR*dimp)]=wmod[fencode3_hde1(p,ii,field)+(ordero*NVAR*dimp)]+dt*dwn1[fencode3_hde1(p,ii,field)]; 

                         }
              //  }	
}
  //__syncthreads();*/



 
}

__global__ void hyperdifesource2_parallel(struct params *p,  real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real *wtemp, int field, int dim)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  //int i = blockIdx.x * blockDim.x + threadIdx.x;
  //int j = blockIdx.y * blockDim.y + threadIdx.y;

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int ii1,ii0;
  real fip,fim1;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real rdx;
  real dy=p->dx[1];
  real dx=p->dx[0];
  //real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

   int ip,jp;
  int ii[NDIM];
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


     ii[0]=ip;
     ii[1]=jp;
     i=ii[0];
     j=ii[1];
     k=0;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
           k=ii[2];
     #endif

     #ifdef USE_SAC_3D
       if(i>0 && j >0 && k>0 && i<((p->n[0])-1) && j<((p->n[1])-1)   && k<((p->n[2])-1))
     #else
       if(i>0 && j >0 && i<((p->n[0])-1) && j<((p->n[1])-1))
     #endif  
  //if(i>0 && j >0 && i<((p->n[0])-1) && j<((p->n[1])-1))
  {
	wtemp[fencode3_hde1(p,ii,tmp2)]= grad1l3n_hde1(wtemp,wd,p,ii,tmp1,dim) ;
	wtemp[fencode3_hde1(p,ii,tmp3)]= grad1r3n_hde1(wtemp,wd,p,ii,tmp1,dim) ;
	//wtemp[fencode3_hde1(p,ii,tmp2)]= -0.0007 ;
	//wtemp[fencode3_hde1(p,ii,tmp3)]= -0.00005 ;
	//wtemp[fencode3_hde1(p,ii,tmp2)]= (  ( wtemp[encode3_hde1(p,i,j,k,rho)]-wtemp[encode3_hde1(p,i-(dim==0),j-(dim==1),k,rho)]) /((p->dx[0]))    ) ;
	//wtemp[fencode3_hde1(p,ii,tmp3)]= (  ( wtemp[encode3_hde1(p,i+(dim==0),j+(dim==1),k,rho)]-wtemp[encode3_hde1(p,i,j,k,rho)]) /((p->dx[0]))    ) ;

  }

//__syncthreads();




 
}



__global__ void hyperdifesource1a_parallel(struct params *p,  real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real *wtemp, int field, int dim)
{

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int ii1,ii0;
  real fip,fim1;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real rdx;
  real dy=p->dx[1];
  real dx=p->dx[0];
   int ip,jp;
  int ii[NDIM];
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


   

     ii[0]=ip;
     ii[1]=jp;
     i=ii[0];
     j=ii[1];
     k=0;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
           k=ii[2];
     #endif

     #ifdef USE_SAC_3D
       if(i<((p->n[0])) && j<((p->n[1])) && k<((p->n[2])))
     #else
       if(i<((p->n[0])) && j<((p->n[1])))
     #endif 
//if(i<((p->n[0])) && j<((p->n[1])))
  {

#ifdef USE_SAC
     wtemp[fencode3_hde1(p,ii,tmp1)]=/*100**/((wmod[shift+fencode3_hde1(p,ii,energy)]-0.5*(


(wmod[shift+fencode3_hde1(p,ii,b1)]*wmod[shift+fencode3_hde1(p,ii,b1)]+wmod[shift+fencode3_hde1(p,ii,b2)]*wmod[shift+fencode3_hde1(p,ii,b2)])

+((wmod[shift+fencode3_hde1(p,ii,mom1)]*wmod[shift+fencode3_hde1(p,ii,mom1)]+wmod[shift+fencode3_hde1(p,ii,mom2)]*wmod[shift+fencode3_hde1(p,ii,mom2)])/(wmod[shift+fencode3_hde1(p,ii,rho)]+wmod[shift+fencode3_hde1(p,ii,rhob)])))));
#endif
#ifdef USE_SAC_3D
     wtemp[fencode3_hde1(p,ii,tmp1)]=wmod[shift+fencode3_hde1(p,ii,energy)]-0.5*((wmod[shift+fencode3_hde1(p,ii,b1)]*wmod[shift+fencode3_hde1(p,ii,b1)]+wmod[shift+fencode3_hde1(p,ii,b2)]*wmod[shift+fencode3_hde1(p,ii,b2)]+wmod[shift+fencode3_hde1(p,ii,b3)]*wmod[shift+fencode3_hde1(p,ii,b3)])

+((wmod[shift+fencode3_hde1(p,ii,mom1)]*wmod[shift+fencode3_hde1(p,ii,mom1)]+wmod[shift+fencode3_hde1(p,ii,mom2)]*wmod[shift+fencode3_hde1(p,ii,mom2)]+wmod[shift+fencode3_hde1(p,ii,mom3)]*wmod[shift+fencode3_hde1(p,ii,mom3)])/(wmod[shift+fencode3_hde1(p,ii,rho)]+wmod[shift+fencode3_hde1(p,ii,rhob)]))
);

#endif
 


  }

//__syncthreads();




 
}





__global__ void hyperdifesource1_parallel(struct params *p,  real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real *wtemp, int field, int dim)
{

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int ii1,ii0;
  real fip,fim1;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real rdx;
  real dy=p->dx[1];
  real dx=p->dx[0];
   int ip,jp;
  int ii[NDIM];
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


     ii[0]=ip;
     ii[1]=jp;
     i=ii[0];
     j=ii[1];
     k=0;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
           k=ii[2];
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
        wtemp[fencode3_hde1(p,ii,f)]=0.0;
    dwn1[fencode3_hde1(p,ii,field)]=0.0;
   }

 //__syncthreads();

 




 
}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_hde1(char *label)
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





int cuhyperdifesource1(struct params **p,  struct params **d_p,   real **d_wmod, real **d_dwn1, real **d_wd, int order,int ordero, real **d_wtemp, int field, int dim,real dt)
{
  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   //cudaSetDevice(selectedDevice);
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 
   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;


     hyperdifesource1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim);
      cudaThreadSynchronize();

     hyperdifesource1a_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim);
      cudaThreadSynchronize();



     hyperdifesource2_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim);
      cudaThreadSynchronize();

     //hyperdifesource3_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim);
    //  cudaThreadSynchronize();

     hyperdifesource4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod, *d_dwn1,  *d_wd, order,ordero,*d_wtemp, field, dim,dt);
      cudaThreadSynchronize();

    /*cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);
    printf("dim hdmean hdmax %d %8.8g %8.8g \n",dim, (*p)->hdmean, (*p)->hdmax);*/
}







