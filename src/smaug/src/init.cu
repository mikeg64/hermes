#include "../include/cudapars.h"

/////////////////////////////////////
// standard imports
/////////////////////////////////////
#include <stdio.h>
#include <math.h>

#include "../../athena.h"
#include "../../defs.h"
/////////////////////////////////////
// kernel function (CUDA device)
/////////////////////////////////////
#include "../include/gradops_i.cuh"
#include "../include/init_user_i.cuh"





/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_i(char *label)
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

int cusync(struct params **p)
{

  #ifdef USE_GPUD
     
         for(int igid=0; igid<((*p)->npe); igid++)
         {
                (*p)->ipe=igid;
                cudaSetDevice((*p)->gpid[igid]) ;
                
  #endif
  cudaThreadSynchronize();
  #ifdef USE_GPUD
                 (*p)->ipe=0;
                 cudaSetDevice((*p)->gpid[0]) ;
          }
  #endif
  return 0;
}

int cusetgpu(struct params **p)
{
  #ifdef USE_GPUD
    if(((*p)->ipe)==-1)
    {
         for(int igid=0; igid<((*p)->npe); igid++)
                (*p)->gpid[igid]=igid ;
    }
    else
      cudaSetDevice((*p)->gpid[(*p)->ipe]) ;
                
  #endif
 
  return 0;
}


int cucopyparamstogpu(MeshS *Mesh, DomainS *pD)
{
 GridS *pG=(pD->Grid);

 GPUParams *p=pG->gpuparams;
 cudaMemcpy(*(pG->gpuparams), *p, sizeof(GPUParams), cudaMemcpyHostToDevice);
  (*p)->qt=pG->time;
  (*p)->dt=pG->dt;

}

int cucopyparamsfromgpu(MeshS *Mesh, DomainS *pD)
{
 GridS *pG=(pD->Grid);

 GPUParams *p=pG->gpuparams;
 cudaMemcpy(*p, *(pG->gpuparams), sizeof(GPUParams), cudaMemcpyDeviceToHost);
  pG->time=(*p)->qt;
  pG->dt=(*p)->dt;

}

int cucopytogpu(MeshS *Mesh, DomainS *pD)
{

GridS *pG=(pD->Grid);
GPUParams *p=pG->gpuparams;
GPUData *gd=pG->gpudata;
ConsS ***U1=pG->U;

GPUparams **d_p;

Real **d_wnew
Real **d_wmod
Real **d_dwn1
Real **d_wd
Real **d_wtemp;
Real **d_wtemp1;
Real **d_wtemp2;

Real *buffer;



//copy params to gpudata


  cudaError_t code;
  int i, Nx3T, Nx2T, Nx1T;

  /* Calculate physical size of grid */
  #ifdef USE_SAC_3D 
  if (pG->Nx[2] > 1)
    Nx3T = pG->Nx[2] + 2*nghost;
  else
    Nx3T = 1;
  #endif

  if (pG->Nx[1] > 1)
    Nx2T = pG->Nx[1] + 2*nghost;
  else
    Nx2T = 1;

  if (pG->Nx[0] > 1)
    Nx1T = pG->Nx[0] + 2*nghost;
  else
    Nx1T = 1;

 buffer=(Real *)calloc(Nx1T*Nx2T,sizeof(Real));


//rho, mom1, mom2, mom3, energy, b1, b2, b3,energyb,rhob,b1b,b2b,b3b
#ifdef USE_SAC_3D


  for(k=0; k<Nx3T; k++)
{
 for(j=0; j<Nx2T; j++)
for(i=0; i<Nx1T; i++)

switch(field)
{
case 0:
     buffer[i+j*Nx1T]=U1.d[i][j][k].d;
     break;
case 1:
     buffer[i+j*Nx1T]=U1.M1[i][j][k];
     break;
case 2:
     buffer[i+j*Nx1T]=U1.M2[i][j][k];
     break;
case 3:
     buffer[i+j*Nx1T]=U1.M3[i][j][k];
     break;
case 4:
     buffer[i+j*Nx1T]=U1.E[i][j][k];
     break;
case 5:
     buffer[i+j*Nx1T]=U1.B1c[i][j][k];
     break;
case 6:
     buffer[i+j*Nx1T]=U1.B2c[i][j][k];
     break;
case 7:
     buffer[i+j*Nx1T]=U1.B3c[i][j][k];
     break;

//background variables
case 8:
     if(Mesh->nstep<=1)
     	buffer[i+j*Nx1T]=U1.Eb[i][j][k];
     break;
case 9:
     if(Mesh->nstep<=1)
     	buffer[i+j*Nx1T]=U1.db[i][j][k];
     break;

case 10:
     if(Mesh->nstep<=1)
     	buffer[i+j*Nx1T]=U1.B1cb[i][j][k];
     break;
case 11:
     if(Mesh->nstep<=1)
     	buffer[i+j*Nx1T]=U1.B2cb[i][j][k];
     break;
case 12:
     if(Mesh->nstep<=1)
     	buffer[i+j*Nx1T]=U1.B3cb[i][j][k];
     break;

}  //end of switch over fields


    code = cudaMemcpy((pG->gpudata)+field*Nx1T*Nx2T*Nx3T+k*Nx1T*Nx2T, buffer, Nx1T*NX2T, cudaMemcpyHostToDevice);

    if(code != cudaSuccess) {
      ath_error("[copy_to_gpu_mem U] error: %s\n", cudaGetErrorString(code));
}


}  //end of loop over dimensions
#endif

#ifdef USE_SAC

 for(j=0; j<Nx2T; j++)
for(i=0; i<Nx1T; i++)
{

switch(field)
{
case 0:
     buffer[field][i+j*Nx1T]=U1.d[i][j][k];
     break;
case 1:
     buffer[field][i+j*Nx1T]=U1.Mx[i][j][k];
     break;
case 2:
     buffer[field][i+j*Nx1T]=U1.My[i][j][k];
     break;
case 3:
     buffer[field][i+j*Nx1T]=U1.E[i][j][k];
     break;
case 4:
     buffer[field][i+j*Nx1T]=U1.B1c[i][j][k];
     break;
case 5:
     buffer[field][i+j*Nx1T]=U1.B2c[i][j][k];
     break;

//background variables
case 6:
     if(Mesh->nstep<=1)
     	buffer[i+j*Nx1T]=U1.Eb[i][j][k];
     break;
case 7:
     if(Mesh->nstep<=1)
     	buffer[i+j*Nx1T]=U1.db[i][j][k];
     break;

case 8:
     if(Mesh->nstep<=1)
     	buffer[i+j*Nx1T]=U1.B1cb[i][j][k];
     break;
case 9:
     if(Mesh->nstep<=1)
     	buffer[i+j*Nx1T]=U1.B2cb[i][j][k];
     break;

}//end of switch over fields
}
    code = cudaMemcpy((pG->gpudata)+field*Nx1T*Nx2T, buffer, Nx1T*NX2T, cudaMemcpyHostToDevice);

    if(code != cudaSuccess) {
      ath_error("[copy_to_gpu_mem U] error: %s\n", cudaGetErrorString(code));
}


#endif





free(buffer);



}






int cucopyfromgpu(MeshS *Mesh, DomainS *pD)
{

GridS *pG=(pD->Grid);
GPUParams *p=pG->gpuparams;
GPUData *gd=pG->gpudata;
ConsS ***U1=pG->U;

GPUparams **d_p;

Real **d_wnew
Real **d_wmod
Real **d_dwn1
Real **d_wd
Real **d_wtemp;
Real **d_wtemp1;
Real **d_wtemp2;

Real *buffer;



//copy params to gpudata


  cudaError_t code;
  int i, Nx3T, Nx2T, Nx1T;

  /* Calculate physical size of grid */
  #ifdef USE_SAC_3D 
  if (pG->Nx[2] > 1)
    Nx3T = pG->Nx[2] + 2*nghost;
  else
    Nx3T = 1;
  #endif

  if (pG->Nx[1] > 1)
    Nx2T = pG->Nx[1] + 2*nghost;
  else
    Nx2T = 1;

  if (pG->Nx[0] > 1)
    Nx1T = pG->Nx[0] + 2*nghost;
  else
    Nx1T = 1;

 buffer=(Real *)calloc(Nx1T*Nx2T,sizeof(Real));


//rho, mom1, mom2, mom3, energy, b1, b2, b3,energyb,rhob,b1b,b2b,b3b
#ifdef USE_SAC_3D



  for(k=0; k<Nx3T; k++)
{


    code = cudaMemcpy(buffer,(pG->gpudata)+field*Nx1T*Nx2T*Nx3T+k*Nx1T*Nx2T,  Nx1T*NX2T, cudaMemcpyDeviceToHost);

    if(code != cudaSuccess) {
      ath_error("[copy_to_gpu_mem U] error: %s\n", cudaGetErrorString(code));



 for(j=0; j<Nx2T; j++)
for(i=0; i<Nx1T; i++)

switch(field)
{
case 0:
     U1.d[i][j][k].d=buffer[i+j*Nx1T];
     break;
case 1:
     U1.M1[i][j][k]=buffer[i+j*Nx1T];
     break;
case 2:
     U1.M2[i][j][k]=buffer[i+j*Nx1T];
     break;
case 3:
     U1.M3[i][j][k]=buffer[i+j*Nx1T];
     break;
case 4:
     U1.E[i][j][k]=buffer[i+j*Nx1T];
     break;
case 5:
     U1.B1c[i][j][k]=buffer[i+j*Nx1T];
     break;
case 6:
     U1.B2c[i][j][k]=buffer[i+j*Nx1T];
     break;
case 7:
     U1.B3c[i][j][k]=buffer[i+j*Nx1T];
     break;

//background variables
case 8:
     if(Mesh->nstep<=1)
     	U1.Eb[i][j][k]=buffer[i+j*Nx1T];
     break;
case 9:
     if(Mesh->nstep<=1)
     	U1.db[i][j][k]=buffer[i+j*Nx1T];
     break;

case 10:
     if(Mesh->nstep<=1)
     	U1.B1cb[i][j][k]=buffer[i+j*Nx1T];
     break;
case 11:
     if(Mesh->nstep<=1)
     	U1.B2cb[i][j][k]=buffer[i+j*Nx1T];
     break;
case 12:
     if(Mesh->nstep<=1)
     	U1.B3cb[i][j][k]=buffer[i+j*Nx1T];
     break;

}  //end of switch over fields



}


}  //end of loop over dimensions
#endif

#ifdef USE_SAC

    code = cudaMemcpy(buffer,(pG->gpudata)+field*Nx1T*Nx2T, Nx1T*NX2T, cudaMemcpyDeviceToHost);

    if(code != cudaSuccess) {
      ath_error("[copy_to_gpu_mem U] error: %s\n", cudaGetErrorString(code));

 for(j=0; j<Nx2T; j++)
for(i=0; i<Nx1T; i++)
{

switch(field)
{
case 0:
     U1.d[i][j][k]=buffer[field][i+j*Nx1T];
     break;
case 1:
     U1.Mx[i][j][k]=buffer[field][i+j*Nx1T];
     break;
case 2:
     U1.My[i][j][k]=buffer[field][i+j*Nx1T];
     break;
case 3:
     U1.E[i][j][k]=buffer[field][i+j*Nx1T];
     break;
case 4:
     U1.B1c[i][j][k]=buffer[field][i+j*Nx1T];
     break;
case 5:
     U1.B2c[i][j][k]=buffer[field][i+j*Nx1T];
     break;

//background variables
case 6:
     if(Mesh->nstep<=1)
     	U1.Eb[i][j][k]=buffer[i+j*Nx1T];
     break;
case 7:
     if(Mesh->nstep<=1)
     	U1.db[i][j][k]=buffer[i+j*Nx1T];
     break;

case 8:
     if(Mesh->nstep<=1)
     	U1.B1cb[i][j][k]=buffer[i+j*Nx1T];
     break;
case 9:
     if(Mesh->nstep<=1)
     	U1.B2cb[i][j][k]=buffer[i+j*Nx1T];
     break;

}//end of switch over fields
}

}


#endif


free(buffer);






}




//int cuinit(struct params **p, struct bparams **bp, Real **wmod,Real **wnew, Real **wd, struct state **state, struct params **d_p, struct bparams **d_bp, Real **d_wnew, Real **d_wmod, Real **d_dwn1, Real **d_wd, struct state **d_state, Real **d_wtemp, Real **d_wtemp1, Real **d_wtemp2)

int cuinit(DomainS *pD)
{

GridS *pG=(pD->Grid);

/*Data and params pointer */
GPUParams **p;
GPUBParams **d_bp;
GPUData *gpudata;

cudaError_t code;
int i,dimp, Nx3T, Nx2T, Nx1T;

GPUparams **d_p;

Real *d_wnew
Real *d_wmod
Real *d_dwn1
Real *d_wd
Real *d_wtemp;
Real *d_wtemp1;
Real *d_wtemp2;




/////////////////////////////////////
  // (1) initialisations:
  //     - perform basic sanity checks
  //     - set device
  /////////////////////////////////////

  printf("in cuinit\n");


  p=&(malloc(sizeof(GPUparams)));
  pG->gpuparams= *p;



  //fill the parameter block with the athena settings


  /* Calculate physical size of grid */
  #ifdef USE_SAC_3D 
  if (pG->Nx[2] > 1)
    Nx3T = pG->Nx[2] + 2*nghost;
  else
    Nx3T = 1;
  #endif


  if (pG->Nx[1] > 1)
    Nx2T = pG->Nx[1] + 2*nghost;
  else
    Nx2T = 1;

  if (pG->Nx[0] > 1)
    Nx1T = pG->Nx[0] + 2*nghost;
  else
    Nx1T = 1;

  (*p)->n[0]=Nx1T;
  (*p)->n[1]=Nx2T;
  #ifdef USE_SAC_3D
    (*p)->n[2]=Nx3T;
  #endif
 
  dimp=(((*p)->n[0]))*(((*p)->n[1]));
   
  #ifdef USE_SAC_3D 
  	dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
  #endif  
   (*p)->rk=0;
   (*p)->mode=0;
   (*p)->ng[0]=nghost;
   (*p)->ng[1]=nghost;
  #ifdef USE_SAC_3D
   (*p)->ng[2]=nghost; 
  #endif 






  (*p)->qt=pG->time;
  (*p)->dt=pG->dt;

  (*p)->xmax[0]=pG->MaxX[0]; 
  (*p)->xmin[0]=pG->MinX[0]; 
  (*p)->dx[0]=pG->dx1; 

 (*p)->xmax[1]=pG->MaxX[1]; 
  (*p)->xmin[1]=pG->MinX[1]; 
  (*p)->dx[1]=pG->dx2; 


  #ifdef USE_SAC_3D
   (*p)->xmax[2]=pG->MaxX[2]; 
  (*p)->xmin[2]=pG->MinX[2]; 
  (*p)->dx[2]=pG->dx3; 
 
  #endif 

/**/
(*p)->chyp3=par_getd("smaug","chyp3");


for(i=0;i<NVAR;i++)
  (*p)->chyp[i]=0.0;

(*p)->chyp[rho]=par_getd("smaug","chyprho");
(*p)->chyp[energy]=par_getd("smaug","chypenergy");
(*p)->chyp[b1]=par_getd("smaug","chypb1");
(*p)->chyp[b2]=par_getd("smaug","chypb2");
(*p)->chyp[mom1]=par_getd("smaug","chypmom1");
(*p)->chyp[mom2]=par_getd("smaug","chypmom2");

  #ifdef USE_SAC_3D
   (*p)->chyp[mom3]=par_getd("smaug","chypmom3");
   (*p)->chyp[b3]=par_getd("smaug","chypb3");
  #endif 

 // (p->boundtype[ii][idir][ibound])=0;  //period=0 mpi=1 mpiperiod=2  cont=3 contcd4=4 fixed=5 symm=6 asymm=7
for(int ii=0; ii<NVAR; ii++)
{
//period=0 mpi=1 mpiperiod=2  cont=3 contcd4=4 fixed=5 symm=6 asymm=7   
(*p)->boundtype[ii][0][0])=par_getd("smaug","boundu1");  
(*p)->boundtype[ii][1][0])=par_getd("smaug","boundu2");  
(*p)->boundtype[ii][0][1])=par_getd("smaug","boundl1");  
(*p)->boundtype[ii][1][1])=par_getd("smaug","boundl2");  
}


  //stores for fixed boundaries
  cudaMalloc((void **)&d_bp,sizeof(GPUBparams));
  #ifdef USE_SAC
	  cudaMalloc((void **) & ((*d_bp)->fixed1) , nghost*nx2T*NVAR*sizeof(Real));
	  cudaMalloc((void **) & ((*d_bp)->fixed2) , nghost*nx1T*NVAR*sizeof(Real));
  #endif

  #ifdef USE_SAC_3D
	  cudaMalloc((void **) & ((*d_bp)->fixed1) , nghost*nx2T*nx3T*NVAR*sizeof(Real));
	  cudaMalloc((void **) & ((*d_bp)->fixed2) , nghost*nx1T*nx3T*NVAR*sizeof(Real));
	  cudaMalloc((void **) & ((*d_bp)->fixed3) , nghost*nx1T*nx2T*NVAR*sizeof(Real));
  #endif

  //allocate the memory for the gpudata
//run=0,scatter=1,gather,init,redistribute

	//if(((*p)->rkon)==1)
	//  cudaMalloc((void**)d_wmod, 6*NVAR*dimp*sizeof(Real));
	//else
	  cudaMalloc((void**)d_wmod, 3*NVAR*dimp*sizeof(Real));

	  cudaMalloc((void**)d_dwn1, NVAR*dimp*sizeof(Real));
	  cudaMalloc((void**)d_wd, NDERV*dimp*sizeof(Real));
	  cudaMalloc((void**)d_wtemp, NTEMP*dimp*sizeof(Real));


	  #ifdef USE_SAC
	  cudaMalloc((void**)d_wtemp1, NTEMP1*(((*p)->n[0])+1)* (((*p)->n[1])+1)*sizeof(Real));
	  cudaMalloc((void**)d_wtemp2, NTEMP2*(((*p)->n[0])+2)* (((*p)->n[1])+2)*sizeof(Real));
	  #endif
	  #ifdef USE_SAC_3D
	  cudaMalloc((void**)d_wtemp1, NTEMP1*(((*p)->n[0])+1)* (((*p)->n[1])+1)* (((*p)->n[2])+1)*sizeof(Real));
	  cudaMalloc((void**)d_wtemp2, NTEMP2*(((*p)->n[0])+2)* (((*p)->n[1])+2)* (((*p)->n[2])+2)*sizeof(Real));
	  #endif



	
	//printf("allocating %d %d %d %d\n",dimp,(*p)->n[0],(*p)->n[1],(*p)->n[2]);
	printf("allocating %d %d %d \n",dimp,(*p)->n[0],(*p)->n[1]);
        pG->gpudata=(GPUData *)malloc(sizeof(GPUData));
        gpudata=pG->gpudata;

        gpudata->d_w= *d_w;
        gpudata->d_wnew= *d_wnew;
        gpudata->d_wn1= *d_wn1;
        gpudata->d_wd= *d_wd;
        gpudata->d_wtemp= *d_wtemp;
        gpudata->d_wtemp1= *d_wtemp1;
        gpudata->d_wtemp2= *d_wtemp2;

        gpudata->d_b= *d_bp;


        //copy params to gpudata
        cudaMemcpy(*(pG->gpuparams), *p, sizeof(GPUParams), cudaMemcpyHostToDevice);
        cucopyparamstogpu(pD);
        cucopytogpu(pD);


	//printf("here1\n");






	 
	    //printf("here2\n");

	    //cudaMemcpy(*d_w, *w, NVAR*dimp*sizeof(Real), cudaMemcpyHostToDevice);
	    cudaMemcpy(*d_wmod, *wmod, 2*(1+(((*p)->rkon)==1))*NVAR*dimp*sizeof(Real), cudaMemcpyHostToDevice);
	    cudaMemcpy(*d_wd, *wd, NDERV*dimp*sizeof(Real), cudaMemcpyHostToDevice);






	//printf("here3\n");






	   // cudaMemcpy(*d_wnew, *wnew, 8*((*p)->n[0])* ((*p)->n[1])*sizeof(Real), cudaMemcpyHostToDevice);
	   // printf("here\n");
	    cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);
	    cudaMemcpy(*d_state, *state, sizeof(struct state), cudaMemcpyHostToDevice);
	    
	    dim3 dimBlock(16, 1);
	    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
	    dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
	   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;
	   

	    printf("calling initialiser\n");
	     init_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wnew, *d_wmod, *d_dwn1,  *d_wd, *d_wtemp, *d_wtemp1, *d_wtemp2);


	    printf("called initialiser\n");
	//cudaMemcpy(*w, *d_w, NVAR*dimp*sizeof(Real), cudaMemcpyDeviceToHost);
if((*p)->mode !=3)
{
	cudaMemcpy(*state, *d_state, sizeof(struct state), cudaMemcpyDeviceToHost);
        cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);
}




  return 0;



}



