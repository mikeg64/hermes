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
#include "../include/gradops_cd2.cuh"
#include "../include/dervfields_cd2.cuh"

#include "../include/usersource_cd2.cuh"

__device__ __host__
real fluxe2(real *dw, real *wd, real *w, real *wmod, struct params *p,int *ii, int dir) {

  real ddc=0;
  real fi, fim1;
  real  fip2=0, fim2=0;

  real ddcx=0,ddcy=0;

   real fluxt=0;


        #ifdef USE_SAC



      		fluxt= +wd[fencode3_cd2(p,ii,ptb)]*grad3dn_cd2(wd,wd,p,ii,vel1+dir,dir);
	fluxt += +w[fencode3_cd2(p,ii,b1b)]*w[fencode3_cd2(p,ii,b1b+dir)]*grad3dn_cd2(wd,wd,p,ii,vel1,0)+w[fencode3_cd2(p,ii,b2b)]*w[fencode3_cd2(p,ii,b1b+dir)]*grad3dn_cd2(wd,wd,p,ii,vel1+1,1);
          #endif


        #ifdef USE_SAC_3D
      		fluxt= +wd[fencode3_cd2(p,ii,ptb)]*grad3dn_cd2(wd,wd,p,ii,vel1+dir,dir);

               fluxt += +w[fencode3_cd2(p,ii,b1b)]*w[fencode3_cd2(p,ii,b1b+dir)]*grad3dn_cd2(wd,wd,p,ii,vel1,0)+w[fencode3_cd2(p,ii,b2b)]*w[fencode3_cd2(p,ii,b1b+dir)]*grad3dn_cd2(wd,wd,p,ii,vel1+1,1)+w[fencode3_cd2(p,ii,b3b)]*w[fencode3_cd2(p,ii,b1b+dir)]*grad3dn_cd2(wd,wd,p,ii,vel1+2,2);
        #endif

  return fluxt;


}



__device__ __host__
int divflux_cd2(real *dw, real *wd, real *w, struct params *p,int *ii,int field,int dir) {

  int direction;
  int status=0;
  real divflux=0;



  dw[fencode3_cd2(p,ii,field)]= grad3dn_cd2(wd,wd,p,ii,flux,dir);

 
 #ifdef USE_SAC

  //commented out to test against vac
  /*if(field==energy)
  {    
     dw[fencode3_cd2(p,ii,field)]+=fluxe2(dw, wd, w, p,ii,dir)-w[fencode3_cd2(p,ii,rho)]*((p->g[dir])*w[fencode3_cd2(p,ii,mom1+dir)]    )/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);
   }*/


 #endif
  return ( status);
}


__device__ __host__
int addenergyterms_cd2(real *dw, real *wd, real *w, real *wmod, struct params *p,int *ii,int field,int dir) {

  int direction;
  int status=0;
  real divflux=0;
  

 #if defined USE_SAC  ||  defined USE_SAC_3D

  
  if(field==energy)
  {    
     //computept3_cd2(w,wd,p,ii);
     //wmod[fencode3_cd2(p,ii,field)]+=fluxe2(dw, wd, wmod, p,ii,dir);/*+w[fencode3_cd2(p,ii,rho)]*((p->g[dir])*w[fencode3_cd2(p,ii,mom1+dir)]    )/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);*/
     

        
              		wmod[fencode3_cd2(p,ii,field)]-= +(p->dt)*wd[fencode3_cd2(p,ii,ptb)]*grad3dn_cd2(wd,wd,p,ii,vel1+dir,dir);
                      //wmod[fencode3_cd2(p,ii,field)]-= +(p->dt)*wd[fencode3_cd2(p,ii,ptb)]*grad3d_cd2(wd,p,ii,vel1+dir,dir);

                    for(int idim=0;idim<NDIM;idim++)
                         wmod[fencode3_cd2(p,ii,field)]+=(p->dt)*wmod[fencode3_cd2(p,ii,b1b+idim)]*wmod[fencode3_cd2(p,ii,b1b+dir)]*grad3dn_cd2(wd,wd,p,ii,vel1+dir,idim);
                        //wmod[fencode3_cd2(p,ii,field)]+=(p->dt)*wmod[fencode3_cd2(p,ii,b1b+idim)]*wmod[fencode3_cd2(p,ii,b1b+dir)]*grad3d_cd2(wd,p,ii,vel1+idim,idim);

		//fluxt= +(((p->gamma)-1)*w[fencode3_cd2(p,ii,energyb)]- 0.5*((p->gamma)-2)*(w[fencode3_cd2(p,ii,b1b)]*w[fencode3_cd2(p,ii,b1b)]+w[fencode3_cd2(p,ii,b2b)]*w[fencode3_cd2(p,ii,b2b)]+w[fencode3_cd2(p,ii,b3b)]*w[fencode3_cd2(p,ii,b3b)]))*grad3d_cd2(wd,p,ii,vel1+dir,dir);

               
               //flux= -(((p->gamma)-1)*w[fencode3_cd2(p,ii,energyb)]- 0.5*((p->gamma)-2)*(w[fencode3_cd2(p,ii,b1b)]*w[fencode3_cd2(p,ii,b1b)]+w[fencode3_cd2(p,ii,b2b)]*w[fencode3_cd2(p,ii,b2b)]+w[fencode3_cd2(p,ii,b3b)]*w[fencode3_cd2(p,ii,b3b)]))*grad3d_cd2(wd,p,ii,vel1+dir,dir);
              // fluxt += +w[fencode3_cd2(p,ii,b1b)]*w[fencode3_cd2(p,ii,b1b+dir)]*grad3d_cd2(wd,p,ii,vel1,0)+w[fencode3_cd2(p,ii,b2b)]*w[fencode3_cd2(p,ii,b1b+dir)]*grad3d_cd2(wd,p,ii,vel1+1,1)+w[fencode3_cd2(p,ii,b3b)]*w[fencode3_cd2(p,ii,b1b+dir)]*grad3d_cd2(wd,p,ii,vel1+2,2);

   }


 #endif
  return ( status);
}

__device__ __host__
int addgrav_cd2(real *dw, real *wd, real *w, real *wmod, struct params *p,int *ii) {

  //int direction;
  int status=0;
  int field,dir;
  //real divflux=0;
  //dw[fencode3_cd2(p,ii,field)]= grad_cd2(wd,p,ii,flux,dir);//+grad_cd2(wd,p,ii,f2,1); 


  for(field=rho;field<NVAR;field++)
  {
    switch(field)
    {
               case mom1:
               case mom2:
                    #ifdef USE_SAC_3D
                    case mom3:
                    #endif  
                         dir=field-mom1;
                         wmod[fencode3_cd2(p,ii,field)]+=(p->dt)* (p->g[dir])*w[fencode3_cd2(p,ii,rho)];

                 break;
                 
                 case energy:
                      for(dir=0; dir<NDIM; dir++)
                        wmod[fencode3_cd2(p,ii,field)]+=(p->dt)*w[fencode3_cd2(p,ii,rho)]*((p->g[dir])*w[fencode3_cd2(p,ii,mom1+dir)]    )/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);

                 break;
                 }                               
                                   
  }
  
 


  return ( status);
}


__device__ __host__
real transportflux_cd2 (real *dw, real *wd, real *w, struct params *p,int *ii,int field, int direction) {

 

  // real fluxt=0;

   //transport flux
   //use versions with velocity less ops may improve performance
        #if defined USE_SAC  || defined USE_SAC_3D
     return(w[fencode3_cd2(p,ii,mom1+direction)]*w[fencode3_cd2(p,ii,field)]/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]));
 // flux= wd[fencode3_cd2(p,ii,vel1+direction)]*w[fencode3_cd2(p,ii,field)];
        #else
     return(w[fencode3_cd2(p,ii,mom1+direction)]*w[fencode3_cd2(p,ii,field)]/w[fencode3_cd2(p,ii,rho)]);
//flux= w[fencode3_cd2(p,ii,vel1+direction)]*w[fencode3_cd2(p,ii,field)];
        #endif

 
}




__device__ __host__
real fluxb1(real *dw, real *wd, real *w, struct params *p,int *ii,int field, int direction) {

 
   real fluxt=0;

       #if defined USE_SAC  || defined USE_SAC_3D

  fluxt= -(w[fencode3_cd2(p,ii,b1+direction)]+w[fencode3_cd2(p,ii,field+(NDIM+2)+direction)])*w[fencode3_cd2(p,ii,mom1+(field-b1))]/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);

fluxt+= (w[fencode3_cd2(p,ii,field+(NDIM+2))])*w[fencode3_cd2(p,ii,mom1+direction)]/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);



         #endif

 
  return fluxt;
}



__device__ __host__
real fluxe1(real *dw, real *wd, real *w, struct params *p,int *ii, int direction) {

  real ddc=0;
  real fi, fim1;
  real  fip2=0, fim2=0;
  real ddcx=0,ddcy=0;

   real fluxt=0;

//computept3_cd2(w,wd,p,ii);

         #if defined USE_SAC


fluxt = w[fencode3_cd2(p,ii,mom1+direction)]*(wd[fencode3_cd2(p,ii,pressuret)]);


fluxt  -= w[fencode3_cd2(p,ii,b1+direction)]*(w[fencode3_cd2(p,ii,b1b)]*w[fencode3_cd2(p,ii,mom1)]+w[fencode3_cd2(p,ii,b2b)]*w[fencode3_cd2(p,ii,mom2)]);
fluxt -= w[fencode3_cd2(p,ii,b1b+direction)]*(w[fencode3_cd2(p,ii,b1)]*w[fencode3_cd2(p,ii,mom1)]+w[fencode3_cd2(p,ii,b2)]*w[fencode3_cd2(p,ii,mom2)]);
fluxt /= (w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);
fluxt += w[fencode3_cd2(p,ii,mom1+direction)]*w[fencode3_cd2(p,ii,energyb)]/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);
fluxt -=w[fencode3_cd2(p,ii,b1+direction)]*(w[fencode3_cd2(p,ii,b1)]*w[fencode3_cd2(p,ii,mom1)]+w[fencode3_cd2(p,ii,b2)]*w[fencode3_cd2(p,ii,mom2)])/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);


         #endif

#ifdef USE_SAC_3D
fluxt = w[fencode3_cd2(p,ii,mom1+direction)]*(wd[fencode3_cd2(p,ii,pressuret)]);


fluxt  -= w[fencode3_cd2(p,ii,b1+direction)]*(w[fencode3_cd2(p,ii,b1b)]*w[fencode3_cd2(p,ii,mom1)]+w[fencode3_cd2(p,ii,b2b)]*w[fencode3_cd2(p,ii,mom2)]+w[fencode3_cd2(p,ii,b3b)]*w[fencode3_cd2(p,ii,mom3)]);
fluxt -= w[fencode3_cd2(p,ii,b1b+direction)]*(w[fencode3_cd2(p,ii,b1)]*w[fencode3_cd2(p,ii,mom1)]+w[fencode3_cd2(p,ii,b2)]*w[fencode3_cd2(p,ii,mom2)]+w[fencode3_cd2(p,ii,b3)]*w[fencode3_cd2(p,ii,mom3)]);
fluxt /= (w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);
fluxt +=w[fencode3_cd2(p,ii,mom1+direction)]*w[fencode3_cd2(p,ii,energyb)]/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);
fluxt -=w[fencode3_cd2(p,ii,b1+direction)]*(w[fencode3_cd2(p,ii,b1)]*w[fencode3_cd2(p,ii,mom1)]+w[fencode3_cd2(p,ii,b2)]*w[fencode3_cd2(p,ii,mom2)]+w[fencode3_cd2(p,ii,b3)]*w[fencode3_cd2(p,ii,mom3)])/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);

#endif

  return fluxt;

}








__device__ __host__
int computefluxe(real *dw, real *wd, real *w, struct params *p,int *ii,int direction) {

  int field;//, direction;
  int status=0;
wd[fencode3_cd2(p,ii,flux)]=0.0;
         #if defined USE_SAC  || defined USE_SAC_3D
	     wd[fencode3_cd2(p,ii,flux)]= transportflux_cd2(dw,wd,w,p,ii,energy,direction)+fluxe1(dw,wd,w,p,ii,direction);
         #endif

        
  return ( status);
}

__device__ __host__
int computefluxb1 (real *dw, real *wd, real *w, struct params *p,int *ii, int field,int direction) {


  int status=0;
wd[fencode3_cd2(p,ii,flux)]=0.0;
        
      if(direction==0)
wd[fencode3_cd2(p,ii,flux)]= 0.0;
      else
 #if defined USE_SAC  || defined USE_SAC_3D  
wd[fencode3_cd2(p,ii,flux)]=  transportflux_cd2(dw,wd,w,p,ii,field,direction)-(w[fencode3_cd2(p,ii,b1+direction)]+w[fencode3_cd2(p,ii,b1b+direction)])*w[fencode3_cd2(p,ii,mom1)]/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)])+ (w[fencode3_cd2(p,ii,b1b)])*w[fencode3_cd2(p,ii,mom1+direction)]/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);//+fluxb1(dw,wd,w,p,ii,field,direction);

         #endif

  return ( status);
}

__device__ __host__
int computefluxb2 (real *dw, real *wd, real *w, struct params *p,int *ii, int field,int direction) {


  int status=0;
   wd[fencode3_cd2(p,ii,flux)]=0.0;      
      if(direction==1)
wd[fencode3_cd2(p,ii,flux)]= 0.0;
else
#if defined USE_SAC  || defined USE_SAC_3D 


wd[fencode3_cd2(p,ii,flux)]= transportflux_cd2(dw,wd,w,p,ii,field,direction)-(w[fencode3_cd2(p,ii,b1+direction)]+w[fencode3_cd2(p,ii,b1b+direction)])*w[fencode3_cd2(p,ii,mom2)]/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)])+ (w[fencode3_cd2(p,ii,b2b)])*w[fencode3_cd2(p,ii,mom1+direction)]/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);//+fluxb1(dw,wd,w,p,ii,field,direction);

         #endif


  return ( status);
}


__device__ __host__
int computefluxb3 (real *dw, real *wd, real *w, struct params *p,int *ii, int field,int direction) {

wd[fencode3_cd2(p,ii,flux)]=0.0;
  int status=0;
 #ifdef USE_SAC_3D
 

      if(direction==2)
wd[fencode3_cd2(p,ii,flux)]= 0.0;
else
wd[fencode3_cd2(p,ii,flux)]= transportflux_cd2(dw,wd,w,p,ii,field,direction)-(w[fencode3_cd2(p,ii,b1+direction)]+w[fencode3_cd2(p,ii,b1b+direction)])*w[fencode3_cd2(p,ii,mom3)]/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)])+ (w[fencode3_cd2(p,ii,b3b)])*w[fencode3_cd2(p,ii,mom1+direction)]/(w[fencode3_cd2(p,ii,rho)]+w[fencode3_cd2(p,ii,rhob)]);//+fluxb1(dw,wd,w,p,ii,field,direction);



 
  #endif
  return ( status);
}



//rho, mom1, mom2, mom3, energy, b1, b2, b3
__device__ __host__
void computeflux_cd2 (real *dw, real *wd, real *w, struct params *p,int *ii, int field,int dir) {

  //int status=0;
  switch(field)
  {
     case energy:
      computefluxe(dw,wd,w,p,ii,dir);
      
      // add the following terms for SAC
      // del((b bb+ bb b).v)+ptb del v - bb bb del v
     break;
     case b1:
      computefluxb1(dw,wd,w,p,ii,field,dir);
     break;
     case b2:
       computefluxb2(dw,wd,w,p,ii,field,dir);



     break;
#ifdef USE_SAC_3D
     case b3:
      computefluxb3(dw,wd,w,p,ii,field,dir);
     break;
#endif
  }
  //return ( status);
}



__global__ void centdiff2a_parallel(struct params *p, struct state *s, real *w, real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real dt,int f,int dir)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j,fid;
 // int index;
  int ni=p->n[0];
  int nj=p->n[1];
  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int nk=p->n[2];
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     

   fid=0;




     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif
                          #ifdef USE_SAC_3D
       				if(ii[0]<((p->n[0])-2) && ii[0]>1 && ii[1]>1 && ii[1]<((p->n[1])-2) && ii[2]>1 && ii[2]<((p->n[2])-2))
     			  #else
       				if(ii[0]<((p->n[0]))-2 && ii[0]>1  && ii[1]>1 && ii[1]<((p->n[1])-2))
     			  #endif
                                divflux_cd2(dwn1,wd,wmod+order*NVAR*dimp,p,ii,f,dir); 



__syncthreads();
                        

                         
}

__global__ void centdiff2b_parallel(struct params *p, struct state *s, real *w, real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real dt,int f,int dir)
{

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j,fid;
 
  int ni=p->n[0];
  int nj=p->n[1];
  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int nk=p->n[2];
   int kp;

   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     


   fid=0;






     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif
  real del;



                        switch(dir)
                        {
                         case 0:

                         #ifdef USE_SAC_3D
       				if(ii[0]<((p->n[0]))  && ii[1]>1 && ii[1]<((p->n[1])-2) && ii[2]>1 && ii[2]<((p->n[2])-2))
     			  #else
       				if(ii[0]<((p->n[0]))   && ii[1]>1 && ii[1]<((p->n[1])-2))
     			  #endif
 
                              wmod[fencode3_cd2(p,ii,f)+(ordero*NVAR*dimp)]=wmod[fencode3_cd2(p,ii,f)+(ordero*NVAR*dimp)]-dt*dwn1[fencode3_cd2(p,ii,f)]; 
                         break;
                         case 1:
                         #ifdef USE_SAC_3D
       				if(ii[0]>1 && ii[0]<((p->n[0])-2)  &&  ii[1]<((p->n[1])) && ii[2]>1 && ii[2]<((p->n[2])-2))
     			  #else
       				if(ii[0]>1 && ii[0]<((p->n[0])-2)   && ii[1]<((p->n[1])) )
     			  #endif
  
                              wmod[fencode3_cd2(p,ii,f)+(ordero*NVAR*dimp)]=wmod[fencode3_cd2(p,ii,f)+(ordero*NVAR*dimp)]-dt*dwn1[fencode3_cd2(p,ii,f)]; 

  
                         break;
                         #ifdef USE_SAC_3D
                         case 2:

 
      			if(ii[0]>1 && ii[0]<((p->n[0])-2)  && ii[1]>1 && ii[1]<((p->n[1])-2)  && ii[2]<((p->n[2])))
 
                              wmod[fencode3_cd2(p,ii,f)+(ordero*NVAR*dimp)]=wmod[fencode3_cd2(p,ii,f)+(ordero*NVAR*dimp)]-dt*dwn1[fencode3_cd2(p,ii,f)]; 
                         break;
                         #endif
                        }



__syncthreads(); 


}

__global__ void centdiff2ci_parallel(struct params *p, struct state *s, real *w, real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real dt,int f,int dir)
{


  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j,fid;
 
  int ni=p->n[0];
  int nj=p->n[1];
  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int nk=p->n[2];
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     


   fid=0;

   //compute pbg used in next source term


     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif



     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
   {
   //computevel3_cd2(wmod+(order*NVAR*dimp),wd,p,ii);  
   computepbg3_cd2(wmod+(ordero*NVAR*dimp),wd,p,ii);  
    }
      
    __syncthreads();





     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif


     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
                            dwn1[fencode3_cd2(p,ii,f)]=0.0;




__syncthreads();



                         
}


__global__ void centdiff2c_parallel(struct params *p, struct state *s, real *w, real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real dt,int f,int dir)
{
 

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j,fid;
 // int index;
  int ni=p->n[0];
  int nj=p->n[1];
  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int nk=p->n[2];
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     


   fid=0;

   //compute pbg used in next source term


     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif







     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif




			//if(i>1 && j >1 && i<(ni-2) && j<(nj-2))
     #ifdef USE_SAC_3D
       if(ii[0]<((p->n[0])-2) && ii[1]<((p->n[1])-2) && ii[2]<((p->n[2])-2)     && ii[0]>1    &&  ii[1]>1   && ii[2]>1   )
     #else
       if(ii[0]<(p->n[0])-2 && ii[1]<(p->n[1])-2)
     #endif
                                addenergyterms_cd2(dwn1,wd,w,wmod+ordero*NVAR*dimp,p,ii,f,dir);

    /* #if(defined(USE_SAC_3D) && defined(USE_USERSOURCE))
       //if(ii[0]<((p->n[0])-2) && ii[1]<((p->n[1])-2) && ii[2]<((p->n[2])-2)     && ii[0]>1    &&  ii[1]>1   && ii[2]>1   )
       if(ii[0]<((p->n[0])) && ii[1]<((p->n[1])) && ii[2]<((p->n[2]))    )
     #endif
     #if(defined(USE_SAC) && defined(USE_USERSOURCE))
       //if(ii[0]<(p->n[0])-2 && ii[1]<(p->n[1])-2)
      if(ii[0]<(p->n[0]) && ii[1]<(p->n[1]))
     #endif

                     #ifdef USE_USERSOURCE
                                addsourceterms2_cd2(dwn1,wd,wmod+ordero*NVAR*dimp,p,s,ii,f,dir); 
                     #endif*/





                /*if( ii[1] <(nj) && ii[0]<(ni) )
                           if(p->ipe==1 && ii[1]==125  && (p->it)==2)
                           {
                               wmod[fencode3_cd2(p,ii,rho)]=0.22113;
 			       w[fencode3_cd2(p,ii,rho)]=0.22113;
                           }*/

               /* if( ii[1] <(nj) && ii[0]<(ni) )
                           if(p->ipe==3 && ii[1]==3  && (p->it)==2)
                           {
                               wmod[fencode3_cd2(p,ii,rho)]=0.22118;
 			       w[fencode3_cd2(p,ii,rho)]=0.22118;
                           }*/

                /*if( ii[1] <(nj) && ii[0]<(ni) )
                           if(p->ipe==1 && ii[1]==127  && (p->it)==2)
                           {
                               wmod[fencode3_cd2(p,ii,rho)]=wmod[fencode_cd2(p,ii[0],ii[1]-4,rho)];
 			       w[fencode3_cd2(p,ii,rho)]= w[fencode_cd2(p,ii[0],ii[1]-4,rho)];
                           }

                if( ii[1] <(nj) && ii[0]<(ni) )
                           if(p->ipe==3 && ii[1]==0  && (p->it)==2)
                           {
                               wmod[fencode3_cd2(p,ii,rho)]=wmod[fencode_cd2(p,ii[0],ii[1]+4,rho)];
 			       w[fencode3_cd2(p,ii,rho)]= w[fencode_cd2(p,ii[0],ii[1]+4,rho)];
                           }

                if( ii[1] <(nj) && ii[0]<(ni) )
                           if(p->ipe==1 && ii[1]==126  && (p->it)==2)
                           {
                               wmod[fencode3_cd2(p,ii,rho)]=wmod[fencode_cd2(p,ii[0],ii[1]-4,rho)];
 			       w[fencode3_cd2(p,ii,rho)]= w[fencode_cd2(p,ii[0],ii[1]-4,rho)];
                           }

                if( ii[1] <(nj) && ii[0]<(ni) )
                           if(p->ipe==3 && ii[1]==1  && (p->it)==2)
                           {
                               wmod[fencode3_cd2(p,ii,rho)]=wmod[fencode_cd2(p,ii[0],ii[1]+4,rho)];
 			       w[fencode3_cd2(p,ii,rho)]= w[fencode_cd2(p,ii[0],ii[1]+4,rho)];
                           }*/



__syncthreads();



                         
}

__global__ void grav_parallel(struct params *p, struct state *s, real *w, real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real dt,int f,int dir)
{
 
  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j,fid;
  int ni=p->n[0];
  int nj=p->n[1];
  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int nk=p->n[2];
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     


   fid=0;

   //compute pbg used in next source term


     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif







     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif


     #ifdef USE_SAC_3D
       if(ii[0]<((p->n[0])-2) && ii[1]<((p->n[1])-2) && ii[2]<((p->n[2])-2)     && ii[0]>1    &&  ii[1]>1   && ii[2]>1   )
     #else
       if(ii[0]<(p->n[0])-2 && ii[1]<(p->n[1])-2)
     #endif
                                addgrav_cd2(dwn1,wd,w,wmod+ordero*NVAR*dimp,p,ii);



__syncthreads();



                         
}

__global__ void source_parallel(struct params *p, struct state *s, real *w, real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real dt)
{
  

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j,fid;
  int f,dir;

  int ni=p->n[0];
  int nj=p->n[1];
  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int nk=p->n[2];
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     


   fid=0;

   //compute pbg used in next source term


     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif







     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif


     #if(defined(USE_SAC_3D) && defined(USE_USERSOURCE))
 
       if(ii[0]<((p->n[0])) && ii[1]<((p->n[1])) && ii[2]<((p->n[2]))    )
     #endif
     #if(defined(USE_SAC) && defined(USE_USERSOURCE))

      if(ii[0]<(p->n[0]) && ii[1]<(p->n[1]))
     #endif

                     #ifdef USE_USERSOURCE
                               addsourceterms2_cd2(dwn1,wd,wmod+ordero*NVAR*dimp,p,s,ii,f,dir); 
                     #endif







	


__syncthreads();



                         
}


__global__ void centdiff2d_parallel(struct params *p, struct state *s, real *w, real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real dt,int f,int dir)
{
 

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j,fid;

  int ni=p->n[0];
  int nj=p->n[1];
  int ii[NDIM];
  int dimp=((p->n[0]))*((p->n[1]));
 #ifdef USE_SAC_3D
   int nk=p->n[2];
   int kp;
   real dz=p->dx[2];
   dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
   int ip,jp;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     


   fid=0;




     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif

			     #ifdef USE_SAC
				   if(ii[0]<ni  && ii[1]<(nj))
			     #endif
			     #ifdef USE_SAC_3D
				  if(ii[0]<ni    && ii[1]<(nj) && ii[2]<(nk))
			     #endif 
				{  
                              wmod[fencode3_cd2(p,ii,f)+(ordero*NVAR*dimp)]=wmod[fencode3_cd2(p,ii,f)+(ordero*NVAR*dimp)]-dt*dwn1[fencode3_cd2(p,ii,f)];
				  
				} 



__syncthreads(); 

                         
}






__global__ void centdiff2_parallel(struct params *p, struct state *s, real *w, real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real dt,int f,int dir)
{
 

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j,fid;
 
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

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     


   fid=0;
   



     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif


                        switch(dir)
                        {
                         case 0:
                          #ifdef USE_SAC_3D
       				if(ii[0]<p->n[0] && ii[1]>1 && ii[1]<(p->n[1]-2) && ii[2]>1 && ii[2]<(p->n[2]-2))
     			  #else
       				if(ii[0]<p->n[0] && ii[1]>1 && ii[1]<(p->n[1]-2))
     			  #endif
                         
                            computeflux_cd2(dwn1,wd,wmod+order*NVAR*dimp,p,ii,f,0); 
                         break;
                         case 1:
                          #ifdef USE_SAC_3D
       				if(ii[1]<p->n[1] && ii[0]>1 && ii[0]<(p->n[0]-2) && ii[2]>1 && ii[2]<(p->n[2]-2))
     			  #else
       				if(ii[1]<p->n[1] && ii[0]>1 && ii[0]<(p->n[0]-2))
     			  #endif
                         
                            computeflux_cd2(dwn1,wd,wmod+order*NVAR*dimp,p,ii,f,1); 
                         break;
                          #ifdef USE_SAC_3D
                         case 2:

       				if(ii[2]<p->n[2] && ii[0]>1 && ii[0]<(p->n[0]-2) && ii[1]>1 && ii[1]<(p->n[1]-2))

                         
                            computeflux_cd2(dwn1,wd,wmod+order*NVAR*dimp,p,ii,f,2); 
                         break;
                         #endif
                        }
              
 

__syncthreads();                        






                         
}


__global__ void centdiff2init_parallel(struct params *p, struct state *s, real *w, real *wmod, 
    real *dwn1, real *wd, int order, int ordero, real dt,int f,int dir)
{
 

  int iindex = blockIdx.x * blockDim.x + threadIdx.x;
  int i,j,fid;

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

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     

   fid=0;
   





     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif

     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
                        {
                            dwn1[fencode3_cd2(p,ii,f)]=0.0;

                               wd[fencode3_cd2(p,ii,flux)]=0.0;

                        }


  __syncthreads();   


                         
}


/////////////////////////////////////
// error checking routine
/////////////////////////////////////
void checkErrors_cd2(char *label)
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




int cucentdiff2(struct params **p, struct params **d_p, struct state **d_s, real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd, int order,int ordero, real dt, int field,int dir)
{
 int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 
 
   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;
   //  cudaMemcpy(*w, *d_w, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
  
    cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);

     centdiff2init_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_s,*d_w, *d_wmod, *d_dwn1,  *d_wd, order,ordero,dt,field,dir);
     cudaThreadSynchronize();

     centdiff2_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_s,*d_w, *d_wmod, *d_dwn1,  *d_wd, order,ordero,dt,field,dir);
     cudaThreadSynchronize();



     centdiff2a_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_s, *d_w, *d_wmod, *d_dwn1,  *d_wd, order,ordero,dt,field,dir);
     cudaThreadSynchronize();

     centdiff2b_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_s, *d_w, *d_wmod, *d_dwn1,  *d_wd, order,ordero,dt,field,dir);
     cudaThreadSynchronize();

 

     centdiff2ci_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_s, *d_w, *d_wmod, *d_dwn1,  *d_wd, order,ordero,dt,field,dir);
     cudaThreadSynchronize();


    centdiff2c_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_s, *d_w, *d_wmod, *d_dwn1,  *d_wd, order,ordero,dt,field,dir);
     cudaThreadSynchronize();


    //cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);
    //printf("source params %G %f %f\n",(*p)->test, (*p)->chyp[0] , (*p)->chyp[1]);
    //printf("source params %G \n",(*p)->test);


     //centdiff2d_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_s, *d_w, *d_wmod, *d_dwn1,  *d_wd, order,ordero,dt,field,dir);
     //cudaThreadSynchronize();


     // cudaMemcpy(*w, *d_w, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
     //cudaMemcpy(*wnew, *d_wnew, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
     //cudaMemcpy(*b, *d_b, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

     //checkErrors("copy data from device");

}

int cugrav(struct params **p, struct params **d_p, struct state **d_s, real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd, int order,int ordero, real dt)
{
 int dimp=(((*p)->n[0]))*(((*p)->n[1]));

  int field=rho;
  int dir=0;   
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 

   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;
  
    cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);


     grav_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_s, *d_w, *d_wmod, *d_dwn1,  *d_wd, order,ordero,dt,field,dir);
     cudaThreadSynchronize();


    //cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);
    //printf("source params %G %f %f %G\n",(*p)->test, (*p)->chyp[0] , (*p)->chyp[1] , (*p)->chyp[2]);



     //centdiff2d_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_s, *d_w, *d_wmod, *d_dwn1,  *d_wd, order,ordero,dt,field,dir);
     //cudaThreadSynchronize();


     // cudaMemcpy(*w, *d_w, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
     //cudaMemcpy(*wnew, *d_wnew, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
     //cudaMemcpy(*b, *d_b, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

     //checkErrors("copy data from device");

}

int cusource(struct params **p, struct params **d_p, struct state **d_s, real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd, int order,int ordero, real dt)
{
 int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   int field=rho;
  int dir=0;     
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif 
  
   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;
   //  cudaMemcpy(*w, *d_w, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
   // if(order==0)
    cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);


     //centdiff2ci_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_s, *d_w, *d_wmod, *d_dwn1,  *d_wd, order,ordero,dt,field,dir);
     //cudaThreadSynchronize();


     source_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_s, *d_w, *d_wmod, *d_dwn1,  *d_wd, order,ordero,dt);
     cudaThreadSynchronize();


    cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);
    //printf("vx vy e %8.16G %8.16G %8.16G\n", (*p)->chyp[0] , (*p)->chyp[1] ,(*p)->test);



     //centdiff2d_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_s, *d_w, *d_wmod, *d_dwn1,  *d_wd, order,ordero,dt,field,dir);
     //cudaThreadSynchronize();


     // cudaMemcpy(*w, *d_w, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
     //cudaMemcpy(*wnew, *d_wnew, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
     //cudaMemcpy(*b, *d_b, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

     //checkErrors("copy data from device");

}
