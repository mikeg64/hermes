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
#include "../include/gradops_nshk.cuh"


__global__ void zeropadmaxdtvisc_parallel(struct params *p,   real *wmod, real *wd, int order, int dir, real *temp, int ndimp)
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
    temp[iindex]=0.0;
     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
             temp[iindex]=(   (p->maxviscoef)+wd[encode3_nshk(p,ii[0],ii[1],ii[2],nushk1+dir)]     )/((wd[encode3_nshk(p,ii[0],ii[1],ii[2],delx1+dir)])*(wd[encode3_nshk(p,ii[0],ii[1],ii[2],delx1+dir)]));//temp[iindex]/(wd[fencode3_cdf(p,ii,delx1+dir)]);

//p->dtdiffvisc=0.25/((maxtmpdt)/((wd[encode3_nshk(p,ii[0],ii[1],ii[2],delx1+dim)])*(wd[encode3_nshk(p,ii[0],ii[1],ii[2],delx1+dim)])));
}


__global__ void newreduction0computemaxdtvisc_parallel(real *cmax, real *temp,int ndimp)
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
      temp[blockIdx.x]=partialResult[0];
     }
    __syncthreads();


}



__global__ void myreduction0computemaxdtvisc_parallel(struct params *p,   real *wmod, real *wd, int order, int dir, real *temp,int ndimp,int s)
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

int tnumThreadsPerBlock = 128;
    
int numBlocks = (dimp+tnumThreadsPerBlock-1) / tnumThreadsPerBlock;
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

   if(iindex==0  && (p->dtdiffvisc<temp[0]))
      p->dtdiffvisc=0.25/temp[0];




 
}



__global__ void getdtvisc_parallel(struct params *p,real *wmod, 
     real *wd, int order, real *wtemp, real *wtemp1, real *wtemp2, int dim)
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
  real dtdiffvisc,tmpdt,maxtmpdt;
  //real g=p->g;
 //  dt=1.0;
//dt=0.05;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

  real maxt=0,max3=0, max1=0;
  
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

   
if(iindex==0)
{
  maxtmpdt=0.0;
  
 // #ifdef USE_SHOCKVISC
    for(ii[0]=0;ii[0]<((p->n[0]));ii[0]++)
      for(ii[1]=0;ii[1]<((p->n[1]));ii[1]++)
     #ifdef USE_SAC_3D
        for(ii[2]=0;ii[2]<((p->n[2]));ii[2]++)
     #endif
 // #endif
	{              
           ;//  tmpdt=(p->maxviscoef)+wd[encode3_nshk(p,ii[0],ii[1],ii[2],nushk1+dim)];

              if(tmpdt>maxtmpdt)
                    maxtmpdt=tmpdt;
              p->dtdiffvisc=0.25/((maxtmpdt)/((wd[encode3_nshk(p,ii[0],ii[1],ii[2],delx1+dim)])*(wd[encode3_nshk(p,ii[0],ii[1],ii[2],delx1+dim)])));
	}

    //p->dtdiffvisc=0.25/((maxtmpdt)/((p->dx[dim])*(p->dx[dim])));                 
//wd[encode3_nshk(p,i,j,k,nushk1+dim)]
   


 //  }
}
 __syncthreads();



 
}





__global__ void nushk2_parallel(struct params *p,real *wmod, 
     real *wd, int order, real *wtemp, real *wtemp1, real *wtemp2, int dim)
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

real cshk=0.5;




   //tmp1  tmp_nuI
 
//compute d3r and d1r
   //tmp2  d3r
    //tmp3 d1r


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
       if(ii[0]>1 && ii[1]>1 && ii[2]>1 && ii[0]<p->n[0] && ii[1]<p->n[1]  && ii[2]<p->n[2])
     #else
       if(ii[0]>1 && ii[1]>1 && ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
 
   //if(i>1 && j>1 && i<((p->n[0])) && j<((p->n[1])))       
   { 

  ;//  wd[encode3_nshk(p,i,j,k,nushk1+dim)]=cshk*(p->dx[dim])*(p->dx[dim])*wtemp[encode3_nshk(p,i,j,k,tmp3)];
  ;//  if(wtemp[encode3_nshk(p,i,j,k,tmp3)]>=0.0)
  ;//      wd[encode3_nshk(p,i,j,k,nushk1+dim)]=0.0;
  ;//  wd[encode3_nshk(p,i,j,k,nushk1+dim)]=fabs(wd[encode3_nshk(p,i,j,k,nushk1+dim)]);
   }

   __syncthreads();




}



__global__ void nushk1a_parallel(struct params *p,real *wmod, 
     real *wd, int order, real *wtemp, real *wtemp1, real *wtemp2, int dim)
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
  __shared__ real wts[512];
  __shared__ real wms[512];





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
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
    //set viscosities
   //if( i<((p->n[0])) && j<((p->n[1])))
   {
     #ifdef USE_SAC_3D
     wtemp[encode3_nshk(p,i,j,k,tmp3)]+=grad13_nshk(wtemp,p,ii,tmp1,dim);
     #else
     wtemp[encode3p2_nshk(p,i,j,0,tmp3)]+=grad13_nshk(wtemp,p,ii,tmp1,dim);
     #endif

   }

   
   __syncthreads();


 
}


__global__ void nushk1_parallel(struct params *p,real *wmod, 
     real *wd, int order, real *wtemp, real *wtemp1, real *wtemp2, int dim)
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
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
    //set viscosities
   //if(i<((p->n[0])) && j<((p->n[1])))
   {


        for(int f=tmp1; f<=tmp8; f++)
                 wtemp[fencode3_nshk(p,ii,f)]=0;


   }



  

   __syncthreads();


 

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
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
    //set viscosities
   //if(i<((p->n[0])) && j<((p->n[1])))
   {

        //for(iv=0;iv<NVAR;iv++)
        //               wms[tid+iv*blockdim]=wmod[fencode_nshk(p,i,j,iv)+shift];
        //wts[tid]=wtemp[fencode_nshk(p,i,j,tmp6)];
        //temp value for viscosity

       //tmp6  tmpnu


		wtemp[fencode3_nshk(p,ii,tmp1)]=wmod[fencode3_nshk(p,ii,mom1+dim)+shift]/(((wmod[fencode3_nshk(p,ii,rho)+shift] +wmod[fencode3_nshk(p,ii,rhob)+shift])));
      ;//   wd[fencode3_nshk(p,ii,nushk1+dim)]=0;       
        //wtemp2[encode3_nshk(p,i+1,j+1,k,tmpnui)]=wtemp[fencode3_nshk(p,ii,tmp6)];






        }
        //wtemp2[encode3_nshk(p,i+1,j+1,k+1,tmpnui)]=wtemp[fencode3_nshk(p,ii,tmp6)];





   


   __syncthreads();




}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_nshk(char *label)
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





int cunushk1(struct params **p,  struct params **d_p,   real **d_wmod,  real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2)
{

  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 

// dim3 dimBlock(dimblock, 1);
 
 //   dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;

    cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);

     for(int dir=0;dir<NDIM;dir++)
     {
         nushk1_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,   *d_wd, order, *d_wtemp,*d_wtemp1,*d_wtemp2, dir);
         cudaThreadSynchronize();
         nushk1a_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,   *d_wd, order, *d_wtemp,*d_wtemp1,*d_wtemp2, dir);
         cudaThreadSynchronize();
     }
     
     for(int dir=0;dir<NDIM;dir++)
     {
         nushk2_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,   *d_wd, order, *d_wtemp,*d_wtemp1,*d_wtemp2, dir);
         cudaThreadSynchronize();
     }


}

int cugetdtvisc1(struct params **p,  struct params **d_p,   real **d_wmod,  real **wd, real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2)
{

  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

  double *d_cmax;
  double *d_bmax;
  int i;


     real fn,fractn,in;
  int ndimp;
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

       int NTPB=512;
  int numBlocks = (ndimp+NTPB-1) / NTPB;

  int smemSize = NTPB * sizeof(double);
 double *h_cmax = (double*)malloc(numBlocks*sizeof(double));

  cudaMalloc((void**)&d_cmax, numBlocks*sizeof(double)); 
  cudaMalloc((void**)&d_bmax, numBlocks*sizeof(double)); 




// dim3 dimBlock(dimblock, 1);
 
 //   dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
 //  int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;
   
           /*for(int dim=0; dim<=(NDIM-1); dim++)
        {
        dtdiffvisc=0.25/(p->maxviscoef/((p->dx[dim])*(p->dx[dim])));
        if(dtdiffvisc>1.0e-8 && (p->dt)>dtdiffvisc )
                                      p->dt=dtdiffvisc;
        }*/

    cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);

 (*p)->dtdiffvisc=0.0000000001;


     for(int dir=0;dir<NDIM;dir++)
     {


     zeropadmaxdtvisc_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir, *d_wtemp,ndimp);
      cudaThreadSynchronize();
	cudaMemcpy(*wd, *d_wd, NDERV*dimp*sizeof(real), cudaMemcpyDeviceToHost);
	//cudaMemcpy(*d_wtemp, ((*wd)+(hdnur*dimp)), dimp*sizeof(real), cudaMemcpyHostToDevice);

  for(i=0;i<numBlocks;i++)
               h_cmax[i]=0;
  cudaMemcpy(d_bmax, h_cmax, numBlocks*sizeof(double), cudaMemcpyHostToDevice);

  newreduction0computemaxdtvisc_parallel<<<numBlocks,NTPB,smemSize>>>(d_bmax,*d_wtemp,ndimp);
  cudaThreadSynchronize();
  cudaMemcpy(h_cmax, d_bmax, numBlocks*sizeof(double), cudaMemcpyDeviceToHost);

  int oldnumBlocks,newnumBlocks;
  newnumBlocks=numBlocks;

  /*while(newnumBlocks>1)
  {
        for(i=0;i<numBlocks;i++)
                h_cmax[i]=0;
        cudaMemcpy(d_bmax, h_cmax, numBlocks*sizeof(double), cudaMemcpyHostToDevice);

       oldnumBlocks=newnumBlocks;
  	newnumBlocks = (newnumBlocks+NTPB-1) / NTPB;

  	newreduction0computemaxdtvisc_parallel<<<newnumBlocks,NTPB,smemSize>>>(d_bmax,*d_wtemp,oldnumBlocks);
       cudaThreadSynchronize();
       cudaMemcpy(h_cmax, d_bmax, newnumBlocks*sizeof(double), cudaMemcpyDeviceToHost);*/

  /*for (i=0; i<numBlocks; i++)
    {
      fprintf(stdout,"cmax# %d %f\n",i, h_cmax[i]);
    }
       fprintf(stdout,"\n");*/


 // }
  for( i=0;i<numBlocks;i++)          		
                if(h_cmax[i]>((*p)->dtdiffvisc)) ((*p)->dtdiffvisc)=h_cmax[i];
  ((*p)->dtdiffvisc)=0.25/((*p)->dtdiffvisc);
;//if(h_cmax[0]>((*p)->dtdiffvisc))
        ;//  (*p)->dtdiffvisc=h_cmax[0];






	/*int s=1;
	while(((s*=2)<=((ndimp/2)-1)) ) 
	{
	   myreduction0computemaxdtvisc_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,  *d_wd, order, dir, *d_wtemp,ndimp,s);
	   cudaThreadSynchronize();
	}*/





         //getdtvisc_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wmod,   *d_wd, order, *d_wtemp,*d_wtemp1,*d_wtemp2, dir);
         cudaThreadSynchronize();
         //cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);
         
         ;//if( ((*p)->dtdiffvisc)>1.0e-8 && ((*p)->dt)>((*p)->dtdiffvisc ))
         ;//                             (*p)->dt=(*p)->dtdiffvisc;
         cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);
         
     }
     

   free(h_cmax);
  cudaFree(d_bmax);
  cudaFree(d_cmax);


}








