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
#include "../include/gradops_b.cuh"


__global__ void boundary_parallel(struct params *p, struct bparams *bp, struct state *s,  real *wmod, int order, int dir, int field)
{

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;

  int ni=p->n[0];
  int nj=p->n[1];
  
  int ip,jp,ipg,jpg;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));

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
     k=0;

      //for( f=rho; f<=b2; f++)
      //{  
      f=field;
	       if(i<((p->n[0])) && j<((p->n[1]))) 
		{
		       bc3_periodic1_dir_b(wmod+order*NVAR*dimp,p,iia,f,dir);  //for OZT
		       //  bc3_cont_cd4_b(wmod+order*NVAR*dimp,p,iia,f);  //for BW
		       //  bc3_fixed_b(wmod+order*NVAR*dimp,p,iia,f,0.0);
		}
       //}

 __syncthreads();





  
}

__global__ void boundary0_parallel(struct params *p, struct bparams *bp, struct state *s,  real *wmod, int order, int dir, int field)
{

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;

  int ni=p->n[0];
  int nj=p->n[1];
  
  int ip,jp,ipg,jpg,kp;
  //int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif   
  
kp=0;

int shift=order*NVAR*dimp;



    // iia[0]=ip;
    // iia[1]=jp;
     k=0;

      //for( f=rho; f<=b2; f++)
      //{  
      //f=field;
	      // if(i<((p->n[0])) && j<((p->n[1]))) 
		//{
		      // bc3_periodic1_dir1_b(wmod+order*NVAR*dimp,p,iia,f,dir);  //for OZT
		       //  bc3_cont_cd4_b(wmod+order*NVAR*dimp,p,iia,f);  //for BW
		       //  bc3_fixed_b(wmod+order*NVAR*dimp,p,iia,f,0.0);
		//}
       //}





kp=0;
f=field;

                if((ip==0 || ip==1) )                
                    wmod[encode3_b(p,ip,jp,kp,f)+shift]=wmod[encode3_b(p,(p->n[0])-4+ip,jp,kp,f)+shift];

                else if(((ip==((p->n[0])-1)) || (ip==((p->n[0])-2))) )                
                    wmod[encode3_b(p,ip,jp,kp,f)+shift]=wmod[encode3_b(p,4-(p->n[0])+ip,jp,kp,f)+shift];




 __syncthreads();




  
}


__global__ void boundary1_parallel(struct params *p, struct bparams *bp, struct state *s,  real *wmod, int order, int dir, int field)
{

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int index,k;
  int f;

  int ni=p->n[0];
  int nj=p->n[1];
  
  int ip,jp,kp,ipg,jpg;
  int iia[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif   
  


int shift=order*NVAR*dimp;



   //  iia[0]=ip;
  //   iia[1]=jp;
   //  k=0;

      //for( f=rho; f<=b2; f++)
      //{  
      
	      // if(i<((p->n[0])) && j<((p->n[1]))) 
		//{
		       //bc3_periodic1_dir2_b(wmod+order*NVAR*dimp,p,iia,f,dir);  //for OZT
		       //  bc3_cont_cd4_b(wmod+order*NVAR*dimp,p,iia,f);  //for BW
		       //  bc3_fixed_b(wmod+order*NVAR*dimp,p,iia,f,0.0);
		//}
       //}



f=field;
kp=0;


                if((jp==0 || jp==1)  )                
                  wmod[shift+encode3_b(p,ip,jp,kp,f)]=wmod[shift+encode3_b(p,ip,(p->n[1])-4+jp,kp,f)];

                else if(((jp==((p->n[1])-1)) || (jp==((p->n[1])-2))) )                 
                  wmod[shift+encode3_b(p,ip,jp,kp,f)]=wmod[shift+encode3_b(p,ip,4-(p->n[1])+jp,kp,f)];




 __syncthreads();





  
}





int cuboundary(struct params **p, struct bparams **bp,struct params **d_p, struct bparams **d_bp, struct state **d_s,  real **d_wmod,  int order,int idir,int field)
{
//printf("bound\n");

 dim3 dimBlock(dimblock, 1);

int numBlocks = ((dimproduct_b(*p)+numThreadsPerBlock-1)) / numThreadsPerBlock;
if((*p)->boundtype[field][0][0]==0)
{
//printf("bound 0\n");
      // boundary_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);
	boundary0_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,0,field);
}
    cudaThreadSynchronize();
if((*p)->boundtype[field][1][0]==0)
{
//printf("bound 1\n");
    //boundary_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
        boundary1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_bp,*d_s, *d_wmod, order,1,field);
}
    cudaThreadSynchronize();
}
