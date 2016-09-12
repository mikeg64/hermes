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
#include "../include/gradops_u.cuh"


__device__ __host__
int updatestate (struct params *p, struct state *s, real *w ,int *ii, int field) {

  int status=0;
                      // atomicExch(&(p->cmax),(wd[fencode3_pre(p,ii,soundspeed)]));
                    switch(field)
                    {
                      case rho:
                    	s->rho=s->rho+(w[fencode3_u(p,ii,field)]);
		      break;
                      case mom1:
                    	s->m1=s->m1+(w[fencode3_u(p,ii,field)]);
		      break;
                      case mom2:
                    	s->m2=s->m2+(w[fencode3_u(p,ii,field)]);
		      break;
                      /*case mom3:
                    	s->m3=s->m3+(w[fencode3_u(p,ii,field)]);
		      break;*/
                      case energy:
                    	s->e=s->e+(w[fencode3_u(p,ii,field)]);
		      break;
                      case b1:
                    	s->b1=s->b1+(w[fencode3_u(p,ii,field)]);
		      break;
                      case b2:
                    	s->b2=s->b2+(w[fencode3_u(p,ii,field)]);
		      break;
                      /*case b3:
                    	s->b3=s->b3+(w[fencode3_u(p,ii,field)]);
		      break;*/
                    };
  return status;
}



__global__ void update_parallel(struct params *p, struct state *s, real *w, real *wmod)
{

   int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j,f;
  int index,k;
  __shared__ int ntot;

  int ni=p->n[0];
  int nj=p->n[1];
  real dt=p->dt;
  //real g=p->g;
  real *u,  *v,  *h;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;

   int ip,jp,ipg,jpg;
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



//int shift=order*NVAR*dimp;

  h=w+dimp*rho;
  u=w+dimp*mom1;
  v=w+dimp*mom2;


     iia[0]=ip;
     iia[1]=jp;
     i=iia[0];
     j=iia[1];
     k=0;
     #ifdef USE_SAC_3D
	   iia[2]=kp;
           k=iia[2];
           for( f=rho; f<=b3; f++)
     #else
           for( f=rho; f<=b2; f++)
     #endif
             {  
         #ifdef USE_SAC_3D
      if(i<((p->n[0])) && j<((p->n[1]))  && k<((p->n[2])))
     #else
       if(i<((p->n[0])) && j<((p->n[1])))
     #endif           
	{
            
                ;//  w[fencode3_u(p,iia,f)]=wmod[fencode3_u(p,iia,f)];
                          //   if(p->ipe==0    && f==rho)
                          //      printf("wmod,w %d %d %lg %lg\n",iia[0],iia[1],wmod[fencode3_u(p,iia,f)],w[fencode3_u(p,iia,f)]);

	}


}

__syncthreads(); 







  
}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_u(char *label)
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
int cuupdate(struct params **p, real **w, real **wmod,real **wtemp2, struct state **state,struct params **d_p, real **d_w, real **d_wmod, real ** d_wtemp2, struct state **d_state, int step)
//int cuupdate(struct params **p, real **w, real **wmod, real **wd, real **temp2, struct state **state,
//             struct params **d_p, real **d_w, real **d_wmod, real **d_wtemp2, struct state **d_state, int step)
{
  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 
    dim3 dimBlock(dimblock, 1);
 
    dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;
  // cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyHostToDevice);
cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);

//no longer necessary as w field no longer used
//just do a memcpy at end of this call
 //    update_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_state,*d_w,*d_wmod);
	    //printf("called update\n"); 
   // cudaThreadSynchronize();
//following comments removed from if def pragmas  if
//using MPI and copying all cell data to host (how slow!?)
//#ifdef USE_MPI

//#else
    if((step%((*p)->cfgsavefrequency))==0)
//#endif
    {

//following commentes removed from section if
//using MPI and copying all cell data to host (how slow!?)
/*#ifdef USE_MPI
    cudaMemcpy(*wmod, *d_w, NVAR*dimp*sizeof(real), cudaMemcpyDeviceToHost);
    #ifdef USE_SAC_3D  
           cudaMemcpy(*wtemp2, *d_wtemp2,NTEMP2*(((*p)->n[0])+2)* (((*p)->n[1])+2)* (((*p)->n[2])+2)*sizeof(real), cudaMemcpyDeviceToHost);
    #else
       cudaMemcpy(*wtemp2, *d_wtemp2,NTEMP2*(((*p)->n[0])+2)* (((*p)->n[1])+2)*sizeof(real), cudaMemcpyDeviceToHost);
    #endif

#endif */ 


 

    //cudaMemcpy(*w, *d_w, NVAR*dimp*sizeof(real), cudaMemcpyDeviceToHost);
    cudaMemcpy(*wmod, *d_wmod, 2*NVAR*dimp*sizeof(real), cudaMemcpyDeviceToHost);






   // cudaMemcpy(*w, *d_w, NVAR*dimp*sizeof(real), cudaMemcpyDeviceToHost);

    //cudaMemcpy(*wnew, *d_wd, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);

   cudaMemcpy(*state, *d_state, sizeof(struct state), cudaMemcpyDeviceToHost);
    }

//cudaMemcpy(*wnew, *d_wnew, 8*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_u, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

  //checkErrors("copy data from device");


 


}


int cuupdatehostwd(struct params **p, real **wd, real **wmod,real **wtemp2, struct state **state,struct params **d_p, real **d_wd, real **d_wmod, real ** d_wtemp2, struct state **d_state, int step)
//int cuupdate(struct params **p, real **w, real **wmod, real **wd, real **temp2, struct state **state,
//             struct params **d_p, real **d_w, real **d_wmod, real **d_wtemp2, struct state **d_state, int step)
{
  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 
    dim3 dimBlock(dimblock, 1);
 
    dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;
  // cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyHostToDevice);
cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);
     //update_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_state,*d_w,*d_wmod);
	    //printf("called update\n"); 
   // cudaThreadSynchronize();



 #ifdef USE_GPUD

	#ifdef USE_SAC_3D
	   int ndimp=((*p)->n[0])*((*p)->n[1])*((*p)->n[2]);
        #else
	   int ndimp= ((*p)->n[0])*((*p)->n[1]);
	#endif      

     real      *wt=(real *)calloc(ndimp*NDERV,sizeof(real));
 

     int shift,oshift;
     int ok1,oj1,oi1;
     int oni,onj,onk;
     int i1,j1,k1;
     int ni,nj,nk;
     real *wa=*wd;

 
     oni=((*p)->n[0])*((*p)->pnpe[0]);
     onj=((*p)->n[1])*((*p)->pnpe[1]);
     ni=((*p)->n[0]);
     nj=((*p)->n[1]);

     #ifdef USE_SAC_3D
     	onk=((*p)->n[2])*((*p)->pnpe[2]);
        nk=((*p)->n[2]);
     #endif

    cudaMemcpy(wt, *d_wd, NDERV*ndimp*sizeof(real), cudaMemcpyDeviceToHost);



     for(int ivar=0; ivar<NDERV; ivar++)
     {

		#ifdef USE_SAC_3D
		   for(k1=0; k1<nk; k1++)
		#endif
        for(j1=0; j1<nj; j1++)
        for(i1=0; i1<ni; i1++)
        {
                oi1=i1+((*p)->pipe[0]*ni);
                oj1=j1+((*p)->pipe[1]*nj);  
		#ifdef USE_SAC_3D
                         shift=(k1*ni*nj+j1*ni+i1);
                         ok1=k1+((*p)->pipe[2]*nk);

                         oshift=(ok1*oni*onj+oj1*oni+oi1);
		#else
			 shift=(j1*ni+i1);
                         oshift=(oj1*oni+oi1);
                #endif
                 //if(i1==0 && j1==0)
                 //if(ivar==0 && ((*p)->ipe)==0 && step==5)
                 // printf("called update %d %d %d %lg %lg\n",ivar,shift,oshift+oni*onj*ivar,wa[oshift+oni*onj*ivar],wt[shift+ivar*ndimp]);//, wa[oshift+oni*onj*ivar]);//,wt[shift]);
                  
                   
              wa[oshift+oni*onj*ivar]=wt[shift+ivar*ndimp];
                                              
        }
     }

       printf("here1\n");   
          free(wt);
         // free(wdt);
#else

 cudaMemcpy(*wd, *d_wd, NDERV*dimp*sizeof(real), cudaMemcpyDeviceToHost);

/*real *wad=*wd;
int iii[3];
iii[2]=0;
printf("update host wd %d\n",(*p)->ipe);
if(((*p)->ipe)==3) 
        for(iii[0]=0; iii[0]<((*p)->n[0]); iii[0]++)
          for(iii[1]=0; iii[1]<((*p)->n[1]); iii[1]++)
             {
               //if(iii[0]==0)
               printf("delx 0 %d %d %16.20f  %16.20f\n",iii[0],iii[1],wad[(fencode3_u(*p,iii,pos1))],wad[(fencode3_u(*p,iii,pos2))]);
             //printf("delx 0 %d %d %d %d\n",iii[0],iii[1],(fencode3_u(*p,iii,pos1)),(fencode3_u(*p,iii,pos2)));
              }*/




   






#endif


  //checkErrors("copy data from device");


 


}



int cuupdatedevicewd(struct params **p, real **wd, real **wmod,real **wtemp2, struct state **state,struct params **d_p, real **d_wd, real **d_wmod, real ** d_wtemp2, struct state **d_state, int step)
//int cuupdate(struct params **p, real **w, real **wmod, real **wd, real **temp2, struct state **state,
//             struct params **d_p, real **d_w, real **d_wmod, real **d_wtemp2, struct state **d_state, int step)
{
  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 
    dim3 dimBlock(dimblock, 1);
 
    dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;
  // cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyHostToDevice);
cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);
     //update_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_state,*d_w,*d_wmod);
	    //printf("called update\n"); 
   // cudaThreadSynchronize();



 #ifdef USE_GPUD

	#ifdef USE_SAC_3D
	   int ndimp=((*p)->n[0])*((*p)->n[1])*((*p)->n[2]);
        #else
	   int ndimp= ((*p)->n[0])*((*p)->n[1]);
	#endif      

     real      *wt=(real *)calloc(ndimp*NDERV,sizeof(real));
 

     int shift,oshift;
     int ok1,oj1,oi1;
     int oni,onj,onk;
     int i1,j1,k1;
     int ni,nj,nk;
     real *wa=*wd;

 
     oni=((*p)->n[0])*((*p)->pnpe[0]);
     onj=((*p)->n[1])*((*p)->pnpe[1]);
     ni=((*p)->n[0]);
     nj=((*p)->n[1]);

    cudaMemcpy(*d_wd,wt, NDERV*ndimp*sizeof(real), cudaMemcpyHostToDevice);

     #ifdef USE_SAC_3D
     	onk=((*p)->n[2])*((*p)->pnpe[2]);
        nk=((*p)->n[2]);
     #endif

    



     for(int ivar=0; ivar<NDERV; ivar++)
     {

		#ifdef USE_SAC_3D
		   for(k1=0; k1<nk; k1++)
		#endif
        for(j1=0; j1<nj; j1++)
        for(i1=0; i1<ni; i1++)
        {
                oi1=i1+((*p)->pipe[0]*ni);
                oj1=j1+((*p)->pipe[1]*nj);  
		#ifdef USE_SAC_3D
                         shift=(k1*ni*nj+j1*ni+i1);
                         ok1=k1+((*p)->pipe[2]*nk);

                         oshift=(ok1*oni*onj+oj1*oni+oi1);
		#else
			 shift=(j1*ni+i1);
                         oshift=(oj1*oni+oi1);
                #endif
                 //if(i1==0 && j1==0)
                 //if(ivar==0 && ((*p)->ipe)==0 && step==5)
                 // printf("called update %d %d %d %lg %lg\n",ivar,shift,oshift+oni*onj*ivar,wa[oshift+oni*onj*ivar],wt[shift+ivar*ndimp]);//, wa[oshift+oni*onj*ivar]);//,wt[shift]);
                  
                   
              wa[oshift+oni*onj*ivar]=wt[shift+ivar*ndimp];
                                              
        }
     }

       printf("here1\n");   
          free(wt);
         // free(wdt);
#else

    cudaMemcpy(*d_wd, *wd, NDERV*dimp*sizeof(real), cudaMemcpyHostToDevice);

#endif


  //checkErrors("copy data from device");


 


}




int cufinish(struct params **p, real **w, real **wnew, struct state **state, struct params **d_p,struct bparams **d_bp, real **d_w, real **d_wnew, real **d_wmod, real **d_dwn1, real **d_wd, struct state **d_state, real **d_wtemp, real **d_wtemp1, real **d_wtemp2)
{
  

 //cudaMemcpy(*w, *d_w, 8*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*wnew, *d_wnew, 8*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_u, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

  checkErrors_u("copy data from device");


  cudaFree(*d_p);
  cudaFree(*d_bp);
//  cudaFree(*d_state);

//  cudaFree(*d_w);
  cudaFree(*d_wnew);
 // cudaFree(*d_u);

  cudaFree(*d_wmod);
  cudaFree(*d_dwn1);
  cudaFree(*d_wd);
  cudaFree(*d_wtemp);
  cudaFree(*d_wtemp1);
  cudaFree(*d_wtemp2);
  




}

  #ifdef USE_MPI

int cufinishmgpu(struct params **p,real **w, real **wmod, real **temp2, real **gmpivisc0, real **gmpivisc1, real **gmpivisc2,   real **gmpiw0, real **gmpiwmod0,   real **gmpiw1, real **gmpiwmod1,   real **gmpiw2, real **gmpiwmod2, struct params **d_p,   real **d_w, real **d_wmod,real **d_wtemp2,    real **d_gmpivisc0,    real **d_gmpivisc1,    real **d_gmpivisc2,   real **d_gmpiw0, real **d_gmpiwmod0,   real **d_gmpiw1, real **d_gmpiwmod1,   real **d_gmpiw2, real **d_gmpiwmod2)
{
  

 //cudaMemcpy(*w, *d_w, 8*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*wnew, *d_wnew, 8*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
//cudaMemcpy(*b, *d_u, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

  //checkErrors_u("copy data from device");


  cudaFree(*d_gmpiw0);
  cudaFree(*d_gmpiwmod0);

  cudaFree(*d_gmpiw1);
  cudaFree(*d_gmpiwmod1);
#ifdef USE_SAC_3D
  cudaFree(*d_gmpiw2);
  cudaFree(*d_gmpiwmod2);
  cudaFree(*d_gmpivisc2);
#endif
  cudaFree(*d_gmpivisc0);
  cudaFree(*d_gmpivisc1);

  //free(*gmpiw0);
  //free(*gmpiwmod0);

 // free(*gmpiw1);
 // free(*gmpiwmod1);
#ifdef USE_SAC_3D
  free(*gmpiw2);
  free(*gmpiwmod2);
free(*gmpivisc2);
#endif

  free(*gmpivisc0);
free(*gmpivisc1);
  //free(*temp2);
}
#endif
