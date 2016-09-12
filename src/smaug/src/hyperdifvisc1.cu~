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
#include "../include/gradops_hdv1.cuh"

__device__ __host__
void bc_hyperdif(real *wt, struct params *p,int *ii, int f,int dir) {

   int i=ii[0];
   int j=ii[1];
   int k=0;
 #ifdef USE_SAC_3D
	k=ii[2];
 #endif

int is=1;
 #ifdef USE_SAC
   if(  (dir == 0) && (i==(p->n[0])-1)   && j>=0   && j<(p->n[1])           )
   {
      //for(int is=0;is<3-2*(j<((p->n[1])-1));is++)   
         wt[encode3p2_hdv1(p,i+2,j+is,k,f)]=wt[encode3p2_hdv1(p,(p->n[0])-5,j+is,k,f)];
         
   }
   else if((dir == 1) && (j==(p->n[1])-1)    && i>0   && i<((p->n[0]))  )
    //for(int is=0;is<3-2*(j<((p->n[1])-1));is++)
       wt[encode3p2_hdv1(p,i+is,j+2,k,f)]=wt[encode3p2_hdv1(p,i+is,(p->n[1])-5,k,f)];
  else if((dir == 0) && (i==0)    && j>0   && j<((p->n[1]))   )
    //for(int is=0;is<3-2*(j<((p->n[1])-1));is++)
       wt[encode3p2_hdv1(p,0,j+is,k,f)]=wt[encode3p2_hdv1(p,6,j+is,k,f)];
   else if((dir == 1) && (j==0)    && i>0   && i<((p->n[0]))   )
    //for(int is=0;is<3-2*(j<((p->n[1])-1));is++)
       wt[encode3p2_hdv1(p,i+is,0,k,f)]=wt[encode3p2_hdv1(p,i+is,6,k,f)];

#endif
 #ifdef USE_SAC_3D
   if(  (dir == 0) && (i==(p->n[0])-1)   && j>0   && j<(p->n[1])      && k>0   && k<(p->n[2])     )
         wt[encode3p2_hdv1(p,i+2,j,k,f)]=wt[encode3p2_hdv1(p,(p->n[0])-5,j,k,f)];
   else if((dir == 1) && (j==(p->n[1])-1)    && i>0   && i<((p->n[0])) && k>0   && k<((p->n[2]))  )
       wt[encode3p2_hdv1(p,i,j+2,k,f)]=wt[encode3p2_hdv1(p,i,(p->n[1])-5,k,f)];
   else if((dir == 2) && (k==(p->n[2])-1)    && i>0   && i<((p->n[0])) && j>0   && j<((p->n[1]))  )
       wt[encode3p2_hdv1(p,i,j,k+2,f)]=wt[encode3p2_hdv1(p,i,j,(p->n[2])-5,f)];
  else if((dir == 0) && (i==0)    && j>0   && j<((p->n[1])) && k>0   && k<((p->n[2]))  )
       wt[encode3p2_hdv1(p,0,j,k,f)]=wt[encode3p2_hdv1(p,6,j,k,f)];
   else if((dir == 1) && (j==0)    && i>0   && i<((p->n[0]))  && k>0   && k<((p->n[2]))  )
       wt[encode3p2_hdv1(p,i,0,k,f)]=wt[encode3p2_hdv1(p,i,6,k,f)];
   else if((dir == 2) && (k==0)    && i>0   && i<((p->n[0])) && j>0   && j<((p->n[1]))   )
       wt[encode3p2_hdv1(p,i,j,0,f)]=wt[encode3p2_hdv1(p,i,j,6,f)];
#endif




 
}


/*__device__ __host__
void bc_periodic1_temp2(real *wt, struct params *p,int i, int j, int f) {

                if(i==1 )                
                    wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,6,j,f)];
                else if((i==((p->n[0]))) )                
                    wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,i-4,j,f)];
                else if(j==1  )                
                  wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,i,6,f)];
                else if((j==((p->n[1]))) )                
                  wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,i,j-4,f)];
}*/

/*__device__ __host__
void bc_periodic2_temp2(real *wt, struct params *p,int i, int j, int f) {


               if(i<1 && j<1)
                {
                  if(i==j)
                    //wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,(p->n[0])-3+i,j,f)];
                    wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,i,6,f)];
                  else                  
                    //wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,i,(p->n[1])-3+j,f)];
                    wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,6,j,f)];                                    
                }
                else if(i<1 && j>((p->n[1])-1))
                {
                  if(i==(j-(p->n[1])-1))                  
                    //wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,(p->n[0])-3+i,4-(p->n[1])+j,f)];
                    wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,6,j,f)];                                     
                  else                  
                    wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,i,j-6,f)];                                     
                }
                else if(i>((p->n[0])-1) && j<1)
                {
                  if((i-(p->n[0])+1)==j)                  
                    wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,i-5,j,f)];                                    
                  else                  
                   wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,i,4,f)];                                    
                }
                else if(i>((p->n[0])-1) && j>((p->n[1])-1))
                {
                  if(i==j)                  
                    wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,i,j-5,f)];                                    
                  else                  
                    wt[fencode_hdv1(p,i,j,f)]=wt[fencode_hdv1(p,i-5,j,f)];                                    
                }                       
                 
                




}*/


__global__ void hyperdifvisc5_parallel(struct params *p,real *wmod, 
     real *wd, int order, real *wtemp, real *wtemp1, real *wtemp2, int field, int dim,int hand)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int is,js;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
  //real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

  real maxt=0,max3=0, max1=0;
  
   int ip,jp,ipg,jpg;
  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp,kpg;
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

   
if(iindex==0)
{
  p->hdmean=0.0;
  p->hdmax=0;
 //  for(ipg=0;ipg<(p->npgp[0]);ipg++)
 //  for(jpg=0;jpg<(p->npgp[1]);jpg++)
  // {

  //   i=ip*(p->npgp[0])+ipg;
 //    j=jp*(p->npgp[1])+jpg;
   //if( i<((p->n[0])) && j<((p->n[1])))
  //if(i>1 && j >1 && i<((p->n[0])-2) && j<((p->n[1])-2))
    //p->cmax=0.0;
    for(ii[0]=1;ii[0]<((p->n[0])-1);ii[0]++)
      for(ii[1]=1;ii[1]<((p->n[1])-1);ii[1]++)
     #ifdef USE_SAC_3D
        for(ii[2]=1;ii[2]<((p->n[2])-1);ii[2]++)
     #endif
	{ 
              // computecmax3_cdf(wmod+(order*dimp*NVAR),wd,p,ii);

             
                    // atomicExch(&(p->cmax),(wd[fencode3_MODID(p,ii,soundspeed)]));
               #ifdef USE_SAC_3D
                if(wd[encode3_hdv1(p,ii[0],ii[1],ii[2],hdnur+hand)]>(p->maxviscoef))
                    p->maxviscoef=(wd[encode3_hdv1(p,ii[0],ii[1],ii[2],hdnur+hand)]);
               #else
                 if(wd[encode3_hdv1(p,ii[0],ii[1],0,hdnur+hand)]>(p->maxviscoef))
                    p->maxviscoef=(wd[encode3_hdv1(p,ii[0],ii[1],0,hdnur+hand)]);
               #endif

              if(wd[encode3_hdv1(p,ii[0],ii[1],0,hdnur+hand)]>(p->hdmax))
                    p->hdmax=(wd[encode3_hdv1(p,ii[0],ii[1],0,hdnur+hand)]);

              p->hdmean=(p->hdmean)+wd[encode3_hdv1(p,ii[0],ii[1],0,hdnur+hand)];

	}
p->hdmean=(p->hdmean)/(((p->n[0])-2)*((p->n[1]))-2);
 //  }
}
 __syncthreads();



 
}







__global__ void hyperdifvisc4_parallel(struct params *p,real *wmod, 
     real *wd, int order, real *wtemp, real *wtemp1, real *wtemp2, int field, int dim,int hand)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int is,js;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
  //real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

  real maxt=0,max3=0, max1=0;
  
   int ip,jp,ipg,jpg;
  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp,kpg;
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


   //tmp1  tmp_nuI
   //tmp2  d3r
    //tmp3 d1r
//tmp4    md3r
//tmp5    md1r
//tmp6    d3l
//tmp7    d1l
//tmp8    md3l
//tmp9    md1l







 //  p->maxviscoef=0;
//  p->cmax=1.0;

    //finally update nur and nul
//tmp4    md3r
//tmp5    md1r
   
   for(ipg=0;ipg<(p->npgp[0]);ipg++)
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
       if(i>1 && i<((p->n[0])-2) && j>1 && j<((p->n[1])-2) && k>1 && k<((p->n[2])-2))
     #else
       if(i>1 && i<((p->n[0])-2) && j>1 && j<((p->n[1])-2))
     #endif
   //if(i>1 && i<((p->n[0])-2) && j>1 && j<((p->n[1])-2))
   {
     //wd[encode3_hdv1(p,i,j,hdnur+hand)]=wtemp2[encode3_hdv1(p,i+1,j+1,tmpnui)];
     if(wtemp[encode3_hdv1(p,i,j,k,tmp5)]>0)
{
//p->cmax=1.0;
     #ifdef USE_SAC_3D
	wd[encode3_hdv1(p,i,j,k,hdnur+hand)]=((dim==0)*(p->dx[0])+(dim==1)*(p->dx[1])+(dim==2)*(p->dx[2]))*(p->cmax)*(p->chyp[field])*wtemp[encode3_hdv1(p,i,j,k,tmp4)]/wtemp[encode3_hdv1(p,i,j,k,tmp5)];
     #else
	wd[encode3_hdv1(p,i,j,k,hdnur+hand)]=((dim==0)*(p->dx[0])+(dim==1)*(p->dx[1]))*(p->cmax)*(p->chyp[field])*wtemp[encode3_hdv1(p,i,j,k,tmp4)]/wtemp[encode3_hdv1(p,i,j,k,tmp5)];
     #endif
       // wd[encode3_hdv1(p,i,j,k,hdnur+hand)]=1.0e-2; 
          //wd[encode3_hdv1(p,i,j,hdnur+hand)]=wtemp[encode3_hdv1(p,i,j,tmp4)];
	//wd[encode3_hdv1(p,i,j,hdnul+hand)]=0.01;
}
     else
        wd[encode3_hdv1(p,i,j,k,hdnur+hand)]=0;


//        wd[encode3_hdv1(p,i,j,k,hdnur+hand)]=1.0e-2; 


   }
}
 __syncthreads();



 
}





__global__ void hyperdifvisc3_parallel(struct params *p,real *wmod, 
     real *wd, int order, real *wtemp, real *wtemp1, real *wtemp2, int field, int dim,int hand)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int is,js,ks;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
  //real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

  real maxt1=0,max3=0, maxt2=0;
  
   int ip,jp,ipg,jpg;
  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp,kpg;
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


   //tmp1  tmp_nuI
   //tmp2  d3r
    //tmp3 d1r
//tmp4    md3r
//tmp5    md1r
//tmp6    d3l
//tmp7    d1l
//tmp8    md3l
//tmp9    md1l





  //compute md3r and md1r
//tmp4    md3r
//tmp5    md1r
  //js=0;
 // is=0;
   for(ipg=0;ipg<(p->npgp[0]);ipg++)
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
       //if(ii[0]>1 && ii[1]>1 && ii[2]>1 && ii[0]<p->n[0] && ii[1]<p->n[1]  && ii[2]<p->n[2])
       if(i>1 && j>1 && k>1 && i<((p->n[0])-2) && j<((p->n[1])-2)   && k<((p->n[2]))-2)
     #else
       //if(ii[0]>1 && ii[1]>1 && ii[0]<p->n[0] && ii[1]<p->n[1])
       if(i>1 && j>1 && i<((p->n[0])-2) && j<((p->n[1])-2))
     #endif

 // if( i>1 && j>1 && i<((p->n[0])-2) && j<((p->n[1])-2))            
   {
         maxt1=0;

     #ifdef USE_SAC_3D
         for(is=-(dim==0); is<=(dim==0); is++)
                for(js=-(dim==1); js<=(dim==1); js++)
                   for(ks=-(dim==2); ks<=(dim==2); ks++)
                {
                   if(wtemp1[encode3p1_hdv1(p,i+1+is,j+1+js,k+1+ks,d3)]>maxt1)
                         maxt1=wtemp1[encode3p1_hdv1(p,i+1+is,j+1+js,k+1+ks,d3)];

                }
	#else
         for(is=-(dim==0); is<=(dim==0); is++)
                for(js=-(dim==1); js<=(dim==1); js++)
                {
                   if(wtemp1[encode3p1_hdv1(p,i+1+is,j+1+js,k,d3)]>maxt1)
                         maxt1=wtemp1[encode3p1_hdv1(p,i+1+is,j+1+js,k,d3)];

                }
	#endif
          wtemp[encode3_hdv1(p,i,j,k,tmp4)]=maxt1;

         maxt2=0;

     #ifdef USE_SAC_3D
         for(is=-2*(dim==0); is<=2*(dim==0); is++)
                for(js=-2*(dim==1); js<=2*(dim==1); js++)
                   for(ks=-(dim==2); ks<=(dim==2); ks++)
                {
                   if(wtemp1[encode3p1_hdv1(p,i+1+is,j+1+js,k+1+ks,d1)]>maxt2)
                        maxt2=wtemp1[encode3p1_hdv1(p,i+1+is,j+1+js,k+1+ks,d1)];

                }
	#else
         for(is=-2*(dim==0); is<=2*(dim==0); is++)
                for(js=-2*(dim==1); js<=2*(dim==1); js++)
                {
                   if(wtemp1[encode3p1_hdv1(p,i+1+is,j+1+js,k,d1)]>maxt2)
                        maxt2=wtemp1[encode3p1_hdv1(p,i+1+is,j+1+js,k,d1)];

                }
	#endif
          wtemp[encode3_hdv1(p,i,j,k,tmp5)]=maxt2;
   }
}
   __syncthreads();







 
}




__global__ void hyperdifvisc2_parallel(struct params *p,real *wmod, 
     real *wd, int order, real *wtemp, real *wtemp1, real *wtemp2, int field, int dim,int hand)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  //int i = blockIdx.x * blockDim.x + threadIdx.x;
  //int j = blockIdx.y * blockDim.y + threadIdx.y;

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j;
  int is,js;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
  //real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

  real maxt=0,max3=0, max1=0;
  
   int ip,jp,ipg,jpg;
  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp,kpg;
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






   //tmp1  tmp_nuI
 
//compute d3r and d1r
   //tmp2  d3r
    //tmp3 d1r

   for(ipg=0;ipg<(p->npgp[0]);ipg++)
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
       if(ii[0]>1 && ii[1]>1 && ii[2]>1 && ii[0]<p->n[0] && ii[1]<p->n[1]  && ii[2]<p->n[2])
     #else
       if(ii[0]>1 && ii[1]>1 && ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
 
   //if(i>1 && j>1 && i<((p->n[0])) && j<((p->n[1])))       
   { 
     if(hand==0)
     {
	#ifdef USE_SAC_3D
		   wtemp1[encode3p1_hdv1(p,i,j,k,d3)]=fabs(3.0*(wtemp2[encode3p2_hdv1(p,i+(dim==0),j+(dim==1),k+(dim==2),tmpnui)] - wtemp2[encode3p2_hdv1(p,i,j,k,tmpnui)] ) - (wtemp2[encode3p2_hdv1(p,i+2*(dim==0),j+2*(dim==1),k+2*(dim==2),tmpnui)] - wtemp2[encode3p2_hdv1(p,i-(dim==0),j-(dim==1),k-(dim==2),tmpnui)]    ));
	#else
		   wtemp1[encode3p1_hdv1(p,i,j,k,d3)]=fabs(3.0*(wtemp2[encode3p2_hdv1(p,i+(dim==0),j+(dim==1),k,tmpnui)] - wtemp2[encode3p2_hdv1(p,i,j,k,tmpnui)] ) - (wtemp2[encode3p2_hdv1(p,i+2*(dim==0),j+2*(dim==1),k,tmpnui)] - wtemp2[encode3p2_hdv1(p,i-(dim==0),j-(dim==1),k,tmpnui)]    ));
	#endif
     }
     else
     {
	#ifdef USE_SAC_3D
		   wtemp1[encode3p1_hdv1(p,i,j,k,d3)]=fabs(3.0*(wtemp2[encode3p2_hdv1(p,i,j,k,tmpnui)] - wtemp2[encode3p2_hdv1(p,i-(dim==0),j-(dim==1),k-(dim==2),tmpnui)]) - (wtemp2[encode3p2_hdv1(p,i+(dim==0),j+(dim==1),k+(dim==2),tmpnui)] - wtemp2[encode3p2_hdv1(p,i-2*(dim==0),j-2*(dim==1),k-2*(dim==2),tmpnui)]    ));
		   //wtemp1[encode3_hdv1(p,i,j,k,d3)]=fabs(3.0*(wtemp2[encode3_hdv1(p,i,j,k,tmpnui)] - wtemp2[encode3_hdv1(p,i-(dim==0),j-(dim==1),k-(dim==2),tmpnui)]) - (wtemp2[encode3_hdv1(p,i+(dim==0),j+(dim==1),k+(dim==2),tmpnui)] - wtemp2[encode3_hdv1(p,i-2*(dim==0),j-2*(dim==1),k-2*(dim==2),tmpnui)]    ));
	#else
		   wtemp1[encode3p1_hdv1(p,i,j,k,d3)]=fabs(3.0*(wtemp2[encode3p2_hdv1(p,i,j,k,tmpnui)] - wtemp2[encode3p2_hdv1(p,i-(dim==0),j-(dim==1),k,tmpnui)] ) - (wtemp2[encode3p2_hdv1(p,i+(dim==0),j+(dim==1),k,tmpnui)] - wtemp2[encode3p2_hdv1(p,i-2*(dim==0),j-2*(dim==1),k,tmpnui)]    ));
	#endif
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
       //if(ii[0]>1 && ii[1]>1 && ii[2]>1 && ii[0]<p->n[0] && ii[1]<p->n[1]  && ii[2]<p->n[2])
       if(i>0 && j>0 && k>0 && i<=((p->n[0])) && j<=((p->n[1]))   && k<=((p->n[2])))
     #else
       //if(ii[0]>1 && ii[1]>1 && ii[0]<p->n[0] && ii[1]<p->n[1])
       if(i>0 && j>0 && i<=((p->n[0])) && j<=((p->n[1])))
     #endif

   //if(i>0 && j>0 && i<=((p->n[0])) && j<=((p->n[1])))            
   { 
     if(hand==0)
     {
     #ifdef USE_SAC_3D
           wtemp1[encode3p1_hdv1(p,i,j,k,d1)]=fabs((wtemp2[encode3p2_hdv1(p,i+(dim==0),j+(dim==1),k+(dim==2),tmpnui)] - wtemp2[encode3p2_hdv1(p,i,j,k,tmpnui)] ));
     #else
           wtemp1[encode3p1_hdv1(p,i,j,k,d1)]=fabs((wtemp2[encode3p2_hdv1(p,i+(dim==0),j+(dim==1),k,tmpnui)] - wtemp2[encode3p2_hdv1(p,i,j,k,tmpnui)] ));
     #endif
     }
     else
     {
     #ifdef USE_SAC_3D
           wtemp1[encode3_hdv1(p,i,j,k,d1)]=fabs((wtemp2[encode3p2_hdv1(p,i,j,k,tmpnui)] - wtemp2[encode3p2_hdv1(p,i-(dim==0),j-(dim==1),k-(dim==2),tmpnui)] ));
     #else
           wtemp1[encode3p1_hdv1(p,i,j,k,d1)]=fabs(wtemp2[encode3p2_hdv1(p,i,j,k,tmpnui)]-(wtemp2[encode3p2_hdv1(p,i-(dim==0),j-(dim==1),k,tmpnui)]   ));
     #endif
     }
   }
}
   __syncthreads();



}



__global__ void hyperdifvisc1a_parallel(struct params *p,real *wmod, 
     real *wd, int order, real *wtemp, real *wtemp1, real *wtemp2, int field, int dim,int hand)
{

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  const int blockdim=blockDim.x;
  const int SZWT=blockdim;
  const int SZWM=blockdim*NVAR;
  int tid=threadIdx.x;
  int i,j,iv;
  int is,js;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];
  real maxt=0,max3=0, max1=0;
  
   int ip,jp,ipg,jpg;
  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp,kpg;
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
  __shared__ real wts[512];
  __shared__ real wms[512];



   for(ipg=0;ipg<(p->npgp[0]);ipg++)
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
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
    //set viscosities
   //if( i<((p->n[0])) && j<((p->n[1])))
   {
     #ifdef USE_SAC_3D
     wtemp2[encode3p2_hdv1(p,i+1,j+1,k+1,tmpnui)]=wtemp[encode3_hdv1(p,i,j,k,tmp6)];
     #else
     wtemp2[encode3p2_hdv1(p,i+1,j+1,0,tmpnui)]=wtemp[encode3_hdv1(p,i,j,0,tmp6)];
     #endif

   }

   }
   __syncthreads();



   for(ipg=0;ipg<(p->npgp[0]);ipg++)
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
       if(ii[0]<p->n[0] && ii[1]<(p->n[1]) && ii[2]<(p->n[2]))
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
    //set viscosities
   //if(i<((p->n[0])) && j<((p->n[1])))
   {
	
        bc_hyperdif(wtemp2, p,ii, tmpnui,dim);

   }


    }
   __syncthreads();





 
}


__global__ void hyperdifvisc1_parallel(struct params *p,real *wmod, 
     real *wd, int order, real *wtemp, real *wtemp1, real *wtemp2, int field, int dim,int hand)
{

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  const int blockdim=blockDim.x;
  const int SZWT=blockdim;
  const int SZWM=blockdim*NVAR;
  int tid=threadIdx.x;
  real maxt=0,max3=0, max1=0;
  int i,j,iv;
  int is,js;
  int index,k;
  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  real dy=p->dx[1];
  real dx=p->dx[0];


  
   int ip,jp,ipg,jpg;



  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int kp,kpg;
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


int bfac1,bfac2,bfac3;
//int bfac1=(field==rho || field>mom2)+(field>rho && field<energy);
//int bfac2= (field==rho || field>mom2);
//int bfac3=(field>rho && field<energy);
int shift=order*NVAR*dimp;
  __shared__ real wts[512];
  __shared__ real wms[512];




//init temp1 and temp2 to zero 
//the compute element initialising n[0] or n[1] element must do +1 and +2
//this is because we fit the problem geometrically to nixnj elements 
   for(ipg=0;ipg<(p->npgp[0]);ipg++)
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
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
    //set viscosities
   //if(i<((p->n[0])) && j<((p->n[1])))
   {


        for(int f=tmp1; f<=tmp8; f++)
                 wtemp[fencode3_hdv1(p,ii,f)]=0;

        for(int f=d1; f<=d3; f++)
     #ifdef USE_SAC_3D
                 wtemp1[encode3p1_hdv1(p,ii[0],ii[1],ii[2],f)]=0;
                 wtemp2[encode3p2_hdv1(p,ii[0],ii[1],ii[2],tmpnui)]=0;
     #else
                 wtemp1[encode3p1_hdv1(p,ii[0],ii[1],k,f)]=0;
                 wtemp2[encode3p2_hdv1(p,ii[0],ii[1],k,tmpnui)]=0;
     #endif

      if(i==((p->n[0])-1))
      {
        for(int f=d1; f<=d3; f++)
                 wtemp1[encode3p1_hdv1(p,ii[0]+1,ii[1],k,f)]=0;
        wtemp2[encode3p2_hdv1(p,i+1,j,k,tmpnui)]=0;
        wtemp2[encode3p2_hdv1(p,i+2,j,k,tmpnui)]=0;
      }
      if(j==((p->n[1])-1))
      {
          for(int f=d1; f<=d3; f++)
                 wtemp1[encode3p1_hdv1(p,i,j+1,k,f)]=0;
          wtemp2[encode3p2_hdv1(p,i,j+1,k,tmpnui)]=0;
          wtemp2[encode3p2_hdv1(p,i,j+2,k,tmpnui)]=0;
      }

     #ifdef USE_SAC_3D
      if(k==((p->n[2])-1))
      {
          for(int f=d1; f<=d3; f++)
                 wtemp1[encode3p1_hdv1(p,i,j,k+1,f)]=0;
          wtemp2[encode3p2_hdv1(p,i,j,k+1,tmpnui)]=0;
          wtemp2[encode3p2_hdv1(p,i,j,k+2,tmpnui)]=0;
      }

     #endif
      if(j==((p->n[1])-1)  && i==((p->n[0])-1))
      {
          for(int f=d1; f<=d3; f++)
                 wtemp1[encode3p1_hdv1(p,i+1,j+1,k,f)]=0;



          for(int di=0; di<2; di++)
             for(int dj=0; dj<2; dj++)
                   wtemp2[encode3p2_hdv1(p,i+1+di,j+1+dj,k,tmpnui)]=0;
               

      }
     #ifdef USE_SAC_3D
      if(i==((p->n[0])-1)  && k==((p->n[2])-1))
      {
          for(int f=d1; f<=d3; f++)
                 wtemp1[encode3p1_hdv1(p,i+1,j,k+1,f)]=0;
          for(int di=0; di<2; di++)
             for(int dk=0; dk<2; dk++)
                   wtemp2[encode3p2_hdv1(p,i+1+di,j,k+1+dk,tmpnui)]=0;


      }
      #endif
     #ifdef USE_SAC_3D
      if(j==((p->n[1])-1)  && k==((p->n[2])-1))
      {
          for(int f=d1; f<=d3; f++)
                 wtemp1[encode3p1_hdv1(p,i+1,j+1,k,f)]=0;

          for(int dk=0; dk<2; dk++)
             for(int dj=0; dj<2; dj++)
                   wtemp2[encode3p2_hdv1(p,i,j+1+dj,k+1+dk,tmpnui)]=0;


      }
      #endif

     #ifdef USE_SAC_3D
      if(i==((p->n[0])-1) && j==((p->n[1])-1)  && k==((p->n[2])-1))
      {
          for(int f=d1; f<=d3; f++)
                 wtemp1[encode3p1_hdv1(p,i+1,j+1,k+1,f)]=0;
       
          for(int dk=0; dk<2; dk++)
             for(int dj=0; dj<2; dj++)
               for(int di=0; di<2; di++)
                   wtemp2[encode3p2_hdv1(p,i+1+di,j+1+dj,k+1+dk,tmpnui)]=0;


      }
      #endif

   }



  }

   __syncthreads();


   for(ipg=0;ipg<(p->npgp[0]);ipg++)
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
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
    //set viscosities
   //if(i<((p->n[0])) && j<((p->n[1])))
   {

        //for(iv=0;iv<NVAR;iv++)
        //               wms[tid+iv*blockdim]=wmod[fencode_hdv1(p,i,j,iv)+shift];
        //wts[tid]=wtemp[fencode_hdv1(p,i,j,tmp6)];
        //temp value for viscosity

       //tmp6  tmpnu
#ifdef USE_SAC
        if(field==energy)
        wtemp[fencode3_hdv1(p,ii,tmp6)]=wmod[fencode3_hdv1(p,ii,energy)+shift]-0.5*((wmod[fencode3_hdv1(p,ii,b1)+shift]*wmod[fencode3_hdv1(p,ii,b1)+shift]+wmod[fencode3_hdv1(p,ii,b2)+shift]*wmod[fencode3_hdv1(p,ii,b2)+shift])+(wmod[fencode3_hdv1(p,ii,mom1)+shift]*wmod[fencode3_hdv1(p,ii,mom1)+shift]+wmod[fencode3_hdv1(p,ii,mom2)+shift]*wmod[fencode3_hdv1(p,ii,mom2)+shift])/(wmod[fencode3_hdv1(p,ii,rho)+shift]+wmod[fencode3_hdv1(p,ii,rhob)+shift] ));
        else
        {
           wtemp[fencode3_hdv1(p,ii,tmp6)]=wmod[fencode3_hdv1(p,ii,field)+shift];
	   if((field ==mom1 || field == mom2))
		wtemp[fencode3_hdv1(p,ii,tmp6)]=wmod[fencode3_hdv1(p,ii,field)+shift]/(((wmod[fencode3_hdv1(p,ii,rho)+shift] +wmod[fencode3_hdv1(p,ii,rhob)+shift])));
        }
        //wtemp2[encode3_hdv1(p,i+1,j+1,k,tmpnui)]=wtemp[fencode3_hdv1(p,ii,tmp6)];



#endif

#ifdef USE_SAC_3D
       if(field==energy)
        wtemp[fencode3_hdv1(p,ii,tmp6)]=wmod[fencode3_hdv1(p,ii,energy)+shift]-0.5*((wmod[fencode3_hdv1(p,ii,b1)+shift]*wmod[fencode3_hdv1(p,ii,b1)+shift]+wmod[fencode3_hdv1(p,ii,b2)+shift]*wmod[fencode3_hdv1(p,ii,b2)+shift]+wmod[fencode3_hdv1(p,ii,b3)+shift]*wmod[fencode3_hdv1(p,ii,b3)+shift])
+(wmod[fencode3_hdv1(p,ii,mom1)+shift]*wmod[fencode3_hdv1(p,ii,mom1)+shift]+wmod[fencode3_hdv1(p,ii,mom2)+shift]*wmod[fencode3_hdv1(p,ii,mom2)+shift]+wmod[fencode3_hdv1(p,ii,mom3)+shift]*wmod[fencode3_hdv1(p,ii,mom3)+shift])/(wmod[fencode3_hdv1(p,ii,rho)+shift]+wmod[fencode3_hdv1(p,ii,rhob)+shift] ));       
       else
       {
          wtemp[fencode3_hdv1(p,ii,tmp6)]=wmod[fencode3_hdv1(p,ii,field)+shift];
	if((field ==mom1 || field == mom2 || field == mom3))
		wtemp[fencode3_hdv1(p,ii,tmp6)]=wmod[fencode3_hdv1(p,ii,field)+shift]/(((wmod[fencode3_hdv1(p,ii,rho)+shift] +wmod[fencode3_hdv1(p,ii,rhob)+shift])));

        }
        //wtemp2[encode3_hdv1(p,i+1,j+1,k+1,tmpnui)]=wtemp[fencode3_hdv1(p,ii,tmp6)];



#endif



        wd[fencode3_hdv1(p,ii,hdnur+hand)]=0;
   }

}
   __syncthreads();




}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_hdv1(char *label)
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





int cuhyperdifvisc1(struct params **p,  struct params **d_p,   real **d_wmod,  real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2, int field, int dim,int hand)
{

  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 

// dim3 dimBlock(dimblock, 1);
 
 //   dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;

    cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);

     hyperdifvisc1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,   *d_wd, order, *d_wtemp,*d_wtemp1,*d_wtemp2, field, dim,hand);
     cudaThreadSynchronize();
     hyperdifvisc1a_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,   *d_wd, order, *d_wtemp,*d_wtemp1,*d_wtemp2, field, dim,hand);
     cudaThreadSynchronize();
     hyperdifvisc2_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,   *d_wd, order, *d_wtemp,*d_wtemp1,*d_wtemp2, field, dim,hand);
     cudaThreadSynchronize();
     hyperdifvisc3_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,   *d_wd, order, *d_wtemp,*d_wtemp1,*d_wtemp2, field, dim,hand);
     cudaThreadSynchronize();
     hyperdifvisc4_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,   *d_wd, order, *d_wtemp,*d_wtemp1,*d_wtemp2, field, dim,hand);
     cudaThreadSynchronize();
   
  /*hyperdifvisc5_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,   *d_wd, order, *d_wtemp,*d_wtemp1,*d_wtemp2, field, dim,hand);
     cudaThreadSynchronize();

    cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);

  if(hand==0)
    printf("field right hdmean hdmax %d %8.8g %8.8g \n",field, (*p)->hdmean, (*p)->hdmax);
  else
    printf("field left hdmean hdmax %d %8.8g %8.8g \n",field, (*p)->hdmean, (*p)->hdmax);*/
}







