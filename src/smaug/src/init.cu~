#include "../include/cudapars.h"
#include "../include/iotypes.h"
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
#include "../include/gradops_i.cuh"
#include "../include/init_user_i.cuh"


//*d_p,*d_w, *d_wnew, *d_wmod, *d_dwn1,  *d_wd

__global__ void init_parallel(struct params *p, real *wnew, real *wmod, 
    real *dwn1, real *wd, real *wtemp, real *wtemp1, real *wtemp2)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  // int i = blockIdx.x * blockDim.x + threadIdx.x;
  // int j = blockIdx.y * blockDim.y + threadIdx.y;

 int iindex = blockIdx.x * blockDim.x + threadIdx.x;
 // int index,k;
int ni=p->n[0];
  int nj=p->n[1];
#ifdef USE_SAC_3D
  int nk=p->n[2];
#endif


// Block index
    int bx = blockIdx.x;
   // int by = blockIdx.y;
    // Thread index
    int tx = threadIdx.x;
   // int ty = threadIdx.y;
    
  real *u,  *v,  *h;

   int ord;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;


  int i,j;
  int ip,jp;
  int ii[NDIM];
   int dimp=((p->n[0]))*((p->n[1]));

   
 #ifdef USE_SAC_3D
   int kp;
  dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
/*   int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni/((p->npgp[1])*(p->npgp[0])));
   jp=(iindex-(kp*(nj*ni/((p->npgp[1])*(p->npgp[0])))))/(ni/(p->npgp[0]));
   ip=iindex-(kp*nj*ni/((p->npgp[1])*(p->npgp[0])))-(jp*(ni/(p->npgp[0])));
#else
    jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));
#endif */ 

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     

   

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
		//b[i+j*(p->n[0])]=0;

                 //Define b	

 


	//apply this special condition
	//initiate alfven wave propagtion 
	//if no initial config read

	    /*for(int f=0; f<NVAR; f++)
            { 		         
                          for(ord=0;ord<(2+3*(p->rkon==1));ord++)
                              wmod[fencode3_i(p,ii,f)+ord*NVAR*dimp]=0;
	    }*/


//	 __syncthreads();

			}

        	
	 __syncthreads();


    /* #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
     
               for(int f=vel1; f<NDERV; f++)
                    wd[fencode3_i(p,ii,f)]=0.0;
     

 __syncthreads(); */



     #ifdef USE_SAC_3D
      // if((p->readini==0) && ii[0]>1 && ii[1]>1  && ii[2]>1 && ii[0]<(p->n[0])-1 && ii[1]<(p->n[1])-1 && ii[2]<(p->n[2])-1)
         if((p->readini==0) && ii[0]<(p->n[0]) && ii[1]<(p->n[1])   && ii[2]<(p->n[2])) 
     #else
      // if((p->readini==0) && ii[0]>2 && ii[1]>2 && ii[0]<(p->n[0])-3 && ii[1]<(p->n[1])-3)  //this form for OZT test???? 
     
     
     //if((p->readini==0) && ii[0]>1 && ii[1]>1  && ii[0]<(p->n[0])-1 && ii[1]<(p->n[1])-1)  //this form for OZT test???? 
        if((p->readini==0) && ii[0]<(p->n[0]) && ii[1]<(p->n[1]))  //this form for BW test  //still issue here
     #endif
	{


            #ifdef ADIABHYDRO
		    if(i> (((p->n[0])/2)-2) && i<(((p->n[0])/2)+2) && j>(((p->n[1])/2)-2) && j<(((p->n[1])/2)+2) ) 
				;//w[fencode3_i(p,ii,rho)]=1.3;
            #else
                   // init_alftest (real *w, struct params *p,int i, int j)
                   // init_alftest(w,p,i,j);
                   // init_ozttest (real *w, struct params *p,int i, int j)
                   // init_ozttest(w,p,i,j);
                   // init_bwtest(w,p,i,j);

	           //default values for positions these may be updated by the initialisation routines
                   wd[fencode3_i(p,ii,delx1)]=(p->dx[0]);
		   wd[fencode3_i(p,ii,delx2)]=(p->dx[1]);
                   wd[fencode3_i(p,ii,pos1)]=(p->xmin[0])+ii[0]*(p->dx[0]);
		   wd[fencode3_i(p,ii,pos2)]=(p->xmin[1])+ii[1]*(p->dx[1]);
                 #ifdef USE_SAC_3D
		   wd[fencode3_i(p,ii,pos3)]=(p->xmin[2])+ii[2]*(p->dx[2]);
                   wd[fencode3_i(p,ii,delx3)]=(p->dx[2]);
                 #endif

                   //init_user_i(w,p,ii);  //initilise using w field

                   //commented out because spicule problem
                   //constructed on host
               if(p->mode!=3)
                   init_user_i(wmod,wd,p,ii);
           #endif

	

        }
	
	 __syncthreads();


       





     #ifdef USE_SAC_3D
       if(ii[0]<p->n[0] && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
	{
        /*for(int f=energyb; f<NVAR; f++)
             if(f != rhob)
                      w[fencode3_i(p,ii,f)]=0.0;*/
        //w[fencode3_i(p,ii,b2b)]=w[fencode3_i(p,ii,b3b)];
        for(int f=rho; f<NVAR; f++)
        {               
                  //wmod[fencode3_i(p,ii,f)]=w[fencode3_i(p,ii,f)];
                  //wmod[  (((3*(1+(p->rkon)))-1)*NVAR*dimp)+fencode3_i(p,ii,f)]=w[fencode3_i(p,ii,f)];              
                  dwn1[fencode3_i(p,ii,f)]=0;

                  //initial value of ord changed to 1 ensure have correct background fields set
                  for(ord=1;ord<(2+3*(p->rkon==1));ord++)
                  {
                              //only the wmod field is used w now redundant
                              wmod[fencode3_i(p,ii,f)+ord*NVAR*dimp]=wmod[fencode3_i(p,ii,f)];

                              //original version using w
                              //wmod[fencode3_i(p,ii,f)+ord*NVAR*dimp]=w[fencode3_i(p,ii,f)];
                              //wmod[fencode3_i(p,ii,b2b)+ord*NVAR*dimp]=w[fencode3_i(p,ii,b3b)];
                  }
  /*        int dir=0;
 for(int field=rho;field<=rho ; field++)
if( ii[0]<4 && (p->ipe)==0  && ((p)->it)==1 && ( isnan(wmod[fencode3_i(p,ii,field)])|| wmod[fencode3_i(p,ii,field)]==0 ))
        { 
    				printf("nant %d %d %d %d %lg %lg \n",ii[0],ii[1],field,dir, wmod[fencode3_i(p,ii,rho)],wmod[fencode3_i(p,ii,field)] );
}*/
                  
        }

        for(int f=tmp1; f<NTEMP; f++)
                 wtemp[fencode3_i(p,ii,f)]=0;


	/*for(int field=rho;field<=rho ; field++)
	if(  (p->ipe)==0  && (  wmod[fencode3_i(p,ii,field)]==0 ))
		{ 
	    				printf("nanti %d %d %d %d %lg %lg \n",ii[0],ii[1],field,0, wmod[fencode3_i(p,ii,rho)],wmod[fencode3_i(p,ii,field)+dimp*NVAR] );
	}*/


}

 __syncthreads();



}

__global__ void updatemod_parallel(struct params *p, real *w, real *wnew, real *wmod, 
    real *dwn1, real *wd, real *wtemp, real *wtemp1, real *wtemp2)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  // int i = blockIdx.x * blockDim.x + threadIdx.x;
  // int j = blockIdx.y * blockDim.y + threadIdx.y;

 int iindex = blockIdx.x * blockDim.x + threadIdx.x;
 // int index,k;
int ni=p->n[0];
  int nj=p->n[1];
#ifdef USE_SAC_3D
  int nk=p->n[2];
#endif


// Block index
    int bx = blockIdx.x;
   // int by = blockIdx.y;
    // Thread index
    int tx = threadIdx.x;
   // int ty = threadIdx.y;
    
  real *u,  *v,  *h;

   int ord;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;


  int i,j;
  int ip,jp;
  int ii[NDIM];
   int dimp=((p->n[0]))*((p->n[1]));

   
 #ifdef USE_SAC_3D
   int kp;
  dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
/*   int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni/((p->npgp[1])*(p->npgp[0])));
   jp=(iindex-(kp*(nj*ni/((p->npgp[1])*(p->npgp[0])))))/(ni/(p->npgp[0]));
   ip=iindex-(kp*nj*ni/((p->npgp[1])*(p->npgp[0])))-(jp*(ni/(p->npgp[0])));
#else
    jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));
#endif */ 

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     

   

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
        for(int f=rho; f<NVAR; f++)
        {               
                  for(ord=1;ord<(2+3*(p->rkon==1));ord++)
                  {
                              //wmod[fencode3_i(p,ii,f)+ord*NVAR*dimp]=w[fencode3_i(p,ii,f)];
                              wmod[fencode3_i(p,ii,f)+ord*NVAR*dimp]=wmod[fencode3_i(p,ii,f)];

                            // if(p->ipe==0    && f==rho)
                            //    printf("wmod,w %d %d %lg %lg\n",ii[0],ii[1],wmod[fencode3_i(p,ii,f)+ord*NVAR*dimp],w[fencode3_i(p,ii,f)]);
 
                  }
          int dir=0;


 //for(int field=rho;field<=rho ; field++)
//if( /*ii[0]<4 &&*/ (p->ipe)==0  && /*((p)->it)==1 &&*/ (/* isnan(wmod[fencode3_i(p,ii,field)])||*/ wmod[fencode3_i(p,ii,field)]==0 ))
//        { 
//    				printf("nant %d %d %d %d %lg %lg \n",ii[0],ii[1],field,dir, wmod[fencode3_i(p,ii,rho)],wmod[fencode3_i(p,ii,field)] );
//}
                  
        }




}

 __syncthreads();



}


 //initialise grid on the gpu
 //we currently don't do this to avoid use of additional memory on GPU
//set up a temporary grid

__global__ void gridsetup_parallel(struct params *p, real *w, real *wnew, real *wmod, 
    real *dwn1, real *wd, real *wtemp, real *wtemp1, real *wtemp2, int dir)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  // int i = blockIdx.x * blockDim.x + threadIdx.x;
  // int j = blockIdx.y * blockDim.y + threadIdx.y;

 int iindex = blockIdx.x * blockDim.x + threadIdx.x;
 // int index,k;
int ni=p->n[0];
  int nj=p->n[1];
#ifdef USE_SAC_3D
  int nk=p->n[2];
#endif


// Block index
    int bx = blockIdx.x;
   // int by = blockIdx.y;
    // Thread index
    int tx = threadIdx.x;
   // int ty = threadIdx.y;
    
  real *u,  *v,  *h;

   int ord;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;


  int i,j;
  int ip,jp,kp;
  int ii[NDIM];
   int dimp=((p->n[0]))*((p->n[1]));
   kp=0;
   
 #ifdef USE_SAC_3D
 
  dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
/*   int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni/((p->npgp[1])*(p->npgp[0])));
   jp=(iindex-(kp*(nj*ni/((p->npgp[1])*(p->npgp[0])))))/(ni/(p->npgp[0]));
   ip=iindex-(kp*nj*ni/((p->npgp[1])*(p->npgp[0])))-(jp*(ni/(p->npgp[0])));
#else
    jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));
#endif */ 

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     

   

     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif


     #ifdef USE_SAC_3D
       if(ii[0]>0 && ii[0]<(p->n[0]-1) && ii[1]<p->n[1] && ii[2]<p->n[2])
     #else
       if(ii[0]<p->n[0] && ii[1]<p->n[1])
     #endif
     {


        switch(dir)
        {

    case 0:
           wtemp2[encode3p2_i(p,ip+1,jp+1,kp+1,tmpnui)]=wd[fencode3_i(p,ii,pos1)];
    break;
    case 1:
           wtemp2[encode3p2_i(p,ip+1,jp+1,kp+1,tmpnui1)]=wd[fencode3_i(p,ii,pos2)];
    break;
    #ifdef USE_SAC_3D
           case 2:
                        wtemp2[encode3p2_i(p,ip+1,jp+1,kp+1,tmpnui2)]=wd[fencode3_i(p,ii,pos3)];
           break;
     #endif
           }
     }


        	
	 __syncthreads();




       





}




 //initialise grid on the gpu
 //we currently don't do this to avoid use of additional memory on GPU
//calculate the dx values

__global__ void setupdx_parallel(struct params *p, real *w, real *wnew, real *wmod, 
    real *dwn1, real *wd, real *wtemp, real *wtemp1, real *wtemp2, int dir)
{
  // compute the global index in the vector from
  // the number of the current block, blockIdx,
  // the number of threads per block, blockDim,
  // and the number of the current thread within the block, threadIdx
  // int i = blockIdx.x * blockDim.x + threadIdx.x;
  // int j = blockIdx.y * blockDim.y + threadIdx.y;

 int iindex = blockIdx.x * blockDim.x + threadIdx.x;
 // int index,k;
int ni=p->n[0];
  int nj=p->n[1];
#ifdef USE_SAC_3D
  int nk=p->n[2];
#endif


// Block index
    int bx = blockIdx.x;
   // int by = blockIdx.y;
    // Thread index
    int tx = threadIdx.x;
   // int ty = threadIdx.y;
    
  real *u,  *v,  *h;

   int ord;
//enum vars rho, mom1, mom2, mom3, energy, b1, b2, b3;


  int i,j;
  int ip,jp,kp;
  int ii[NDIM];
   int dimp=((p->n[0]))*((p->n[1]));

   
 #ifdef USE_SAC_3D
 
  dimp=((p->n[0]))*((p->n[1]))*((p->n[2]));
#endif  
/*   int ip,jp,ipg,jpg;

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni/((p->npgp[1])*(p->npgp[0])));
   jp=(iindex-(kp*(nj*ni/((p->npgp[1])*(p->npgp[0])))))/(ni/(p->npgp[0]));
   ip=iindex-(kp*nj*ni/((p->npgp[1])*(p->npgp[0])))-(jp*(ni/(p->npgp[0])));
#else
    jp=iindex/(ni/(p->npgp[0]));
   ip=iindex-(jp*(ni/(p->npgp[0])));
#endif */ 

  #ifdef USE_SAC_3D
   kp=iindex/(nj*ni);
   jp=(iindex-(kp*(nj*ni)))/ni;
   ip=iindex-(kp*nj*ni)-(jp*ni);
#else
    jp=iindex/ni;
   ip=iindex-(jp*ni);
#endif     

   

     ii[0]=ip;
     ii[1]=jp;
     #ifdef USE_SAC_3D
	   ii[2]=kp;
     #endif

   //calculate the dx values


	    switch(dir)
	    {
		     case 0:
		     #ifdef USE_SAC_3D
		       if( ii[0]>0 && ii[0]<(p->n[0])+1 && ii[1]>0 &&  ii[1]<(p->n[1])+1 && ii[2]>0 &&  ii[2]<(p->n[2])+1)
		     #else
		       if( ii[0]>0 && ii[0]<(p->n[0])+1  && ii[1]>0 && ii[1]<(p->n[1])+1)
		     #endif
	                wd[fencode3_i(p,ii,delx1)]=0.5*(wtemp2[encode3p2_i(p,ip+1,jp,kp,tmpnui)]-wtemp2[encode3p2_i(p,ip-1,jp,kp,tmpnui)]);
		     break;
	
		     case 1:
		     #ifdef USE_SAC_3D
		       if(ii[0]>0 && ii[0]<(p->n[0])+1 && ii[1]>0 &&  ii[1]<(p->n[1])+1 && ii[2]>0 &&  ii[2]<(p->n[2])+1)
		     #else
		       if(ii[0]>0 && ii[0]<(p->n[0])+1 && ii[1]>0 && ii[1]<(p->n[1])+1)
		     #endif
			wd[fencode3_i(p,ii,delx2)]=0.5*(wtemp2[encode3p2_i(p,ip,jp+1,kp,tmpnui)]-wtemp2[encode3p2_i(p,ip,jp-1,kp,tmpnui)]);
		     break;
		         
		     #ifdef USE_SAC_3D
		     case 2:

		       if(ii[0]>0 && ii[0]<(p->n[0])+1 && ii[1]>0 && ii[1]<(p->n[1])+1 && ii[2]>0 && ii[2]<(p->n[2])+1)
			wd[fencode3_i(p,ii,delx3)]=0.5*(wtemp2[encode3p2_i(p,ip,jp,kp+1,tmpnui)]-wtemp2[encode3p2_i(p,ip,jp,kp-1,tmpnui)]);
		     break;			
		     #endif
	     }
     
        	
	 __syncthreads();







       





}

 //initialise grid on the gpu
 //we currently don't do this to avoid use of additional memory on GPU
//intialise temporrary matrix needs t be completed
__global__ void zerotempv_parallel(struct params *p, real *w, real *wnew, real *wmod, 
real *dwn1,  real *wd, real *wtemp, real *wtemp1, real *wtemp2,  int dir)
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
//int shift=order*NVAR*dimp;




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


        for(int f=d1; f<=d3; f++)
     #ifdef USE_SAC_3D
                 wtemp2[encode3p2_i(p,ii[0],ii[1],ii[2],tmpnui)]=0;
     #else
                 wtemp2[encode3p2_i(p,ii[0],ii[1],k,tmpnui)]=0;
     #endif

      if(i==((p->n[0])-1))
      {
        wtemp2[encode3p2_i(p,i+1,j,k,tmpnui)]=0;
        wtemp2[encode3p2_i(p,i+2,j,k,tmpnui)]=0;
      }
      if(j==((p->n[1])-1))
      {
          wtemp2[encode3p2_i(p,i,j+1,k,tmpnui)]=0;
          wtemp2[encode3p2_i(p,i,j+2,k,tmpnui)]=0;
      }

     #ifdef USE_SAC_3D
      if(k==((p->n[2])-1))
      {
          wtemp2[encode3p2_i(p,i,j,k+1,tmpnui)]=0;
          wtemp2[encode3p2_i(p,i,j,k+2,tmpnui)]=0;
      }

     #endif
      if(j==((p->n[1])-1)  && i==((p->n[0])-1))
      {
          for(int di=0; di<2; di++)
             for(int dj=0; dj<2; dj++)
                   wtemp2[encode3p2_i(p,i+1+di,j+1+dj,k,tmpnui)]=0;
      }
     #ifdef USE_SAC_3D
      if(i==((p->n[0])-1)  && k==((p->n[2])-1))
      {
          for(int di=0; di<2; di++)
             for(int dk=0; dk<2; dk++)
                   wtemp2[encode3p2_i(p,i+1+di,j,k+1+dk,tmpnui)]=0;
      }
      #endif

    

     #ifdef USE_SAC_3D
      if(j==((p->n[1])-1)  && k==((p->n[2])-1))
      {
          for(int dk=0; dk<2; dk++)
             for(int dj=0; dj<2; dj++)
                   wtemp2[encode3p2_i(p,i,j+1+dj,k+1+dk,tmpnui)]=0;
      }
      #endif

     #ifdef USE_SAC_3D
      if(i==((p->n[0])-1) && j==((p->n[1])-1)  && k==((p->n[2])-1))
      {
          for(int dk=0; dk<2; dk++)
             for(int dj=0; dj<2; dj++)
               for(int di=0; di<2; di++)
                   wtemp2[encode3p2_i(p,i+1+di,j+1+dj,k+1+dk,tmpnui)]=0;
      }
      #endif

   }

}



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

int cuinit(struct params **p, struct bparams **bp, real **wmod,real **wnew, real **wd, struct state **state, struct params **d_p, struct bparams **d_bp, real **d_wnew, real **d_wmod, real **d_dwn1, real **d_wd, struct state **d_state, real **d_wtemp, real **d_wtemp1, real **d_wtemp2)
{



/////////////////////////////////////
  // (1) initialisations:
  //     - perform basic sanity checks
  //     - set device
  /////////////////////////////////////
  int deviceCount;
  int dir;
 /* cudaGetDeviceCount(&deviceCount);
   
 // if (deviceCount == 0)
 // {
 //   fprintf(stderr, "Sorry, no CUDA device fount");
 //   return 1;
//  }

  #ifdef USE_MPI
     int lipe=(*p)->ipe;
     int gpugid=lipe/4;
     selectedDevice=lipe-4*gpugid;
  #endif
  if (selectedDevice >= deviceCount)
  {
    fprintf(stderr, "Choose device ID between 0 and %d\n", deviceCount-1);
    return 1;
  }


        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, selectedDevice);
        if (deviceProp.major < 1) {
            fprintf(stderr, "gpuDeviceInit(): GPU device does not support CUDA.\n");
            exit(-1);                                                  \
        }

        cudaSetDevice(selectedDevice) ;
        printf("> gpuDeviceInit() CUDA device [%d]: %s %s\n", selectedDevice, deviceProp.name, getenv("HOSTNAME"));



  cudaSetDevice(selectedDevice);
  printf("device count %d selected %d\n", deviceCount,selectedDevice);
  checkErrors_i("initialisations");*/
  
	// Build empty u, v, b matrices

  printf("in cuinit\n");
 // real *adb;
  real *adw, *adwnew;
  struct params *adp;
  struct bparams *adbp;
  struct state *ads;


 
  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   
 #ifdef USE_SAC_3D
   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif  
if((*p)->mode != 3)
{
	if(((*p)->rkon)==1)
	  cudaMalloc((void**)d_wmod, 6*NVAR*dimp*sizeof(real));
	else
	  cudaMalloc((void**)d_wmod, 3*NVAR*dimp*sizeof(real));

	  cudaMalloc((void**)d_dwn1, NVAR*dimp*sizeof(real));
	  cudaMalloc((void**)d_wd, NDERV*dimp*sizeof(real));
	  cudaMalloc((void**)d_wtemp, NTEMP*dimp*sizeof(real));


	  #ifdef USE_SAC
	  cudaMalloc((void**)d_wtemp1, NTEMP1*(((*p)->n[0])+1)* (((*p)->n[1])+1)*sizeof(real));
	  cudaMalloc((void**)d_wtemp2, NTEMP2*(((*p)->n[0])+2)* (((*p)->n[1])+2)*sizeof(real));
	  #endif
	  #ifdef USE_SAC_3D
	  cudaMalloc((void**)d_wtemp1, NTEMP1*(((*p)->n[0])+1)* (((*p)->n[1])+1)* (((*p)->n[2])+1)*sizeof(real));
	  cudaMalloc((void**)d_wtemp2, NTEMP2*(((*p)->n[0])+2)* (((*p)->n[1])+2)* (((*p)->n[2])+2)*sizeof(real));
	  #endif

	  //cudaMalloc((void**)&adw, NVAR*dimp*sizeof(real));
	  //cudaMalloc((void**)&adwnew, NVAR*dimp*sizeof(real));

	  cudaMalloc((void**)&adbp, sizeof(struct bparams));
	  cudaMalloc((void**)&adp, sizeof(struct params));
	  cudaMalloc((void**)&ads, sizeof(struct state));
	 // checkErrors_i("memory allocation");

	printf("ni is %d\n",(*p)->n[1]);

	   // *d_b=adb;
	    *d_bp=adbp;
	    *d_p=adp;
	    //*d_w=adw;
	    //*d_wnew=adwnew;
	    *d_state=ads;

	     
	//printf("allocating %d %d %d %d\n",dimp,(*p)->n[0],(*p)->n[1],(*p)->n[2]);
	printf("allocating %d %d %d \n",dimp,(*p)->n[0],(*p)->n[1]);




	//printf("here1\n");






	 
	    //printf("here2\n");

	    //cudaMemcpy(*d_w, *w, NVAR*dimp*sizeof(real), cudaMemcpyHostToDevice);
	    cudaMemcpy(*d_wmod, *wmod, 2*(1+(((*p)->rkon)==1))*NVAR*dimp*sizeof(real), cudaMemcpyHostToDevice);
	    cudaMemcpy(*d_wd, *wd, NDERV*dimp*sizeof(real), cudaMemcpyHostToDevice);






	//printf("here3\n");






	   // cudaMemcpy(*d_wnew, *wnew, 8*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyHostToDevice);
	   // printf("here\n");
	    cudaMemcpy(*d_p, *p, sizeof(struct params), cudaMemcpyHostToDevice);
	    cudaMemcpy(*d_state, *state, sizeof(struct state), cudaMemcpyHostToDevice);
	    
	    dim3 dimBlock(16, 1);
	    //dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
	    dim3 dimGrid(((*p)->n[0])/dimBlock.x,((*p)->n[1])/dimBlock.y);
	   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;
	   

	    printf("calling initialiser\n");
	     //init_parallel(struct params *p, real *b, real *u, real *v, real *h)
	    // init_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_b,*d_u,*d_v,*d_h);
	    // init_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w, *d_wnew, *d_b);
	     init_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p, *d_wnew, *d_wmod, *d_dwn1,  *d_wd, *d_wtemp, *d_wtemp1, *d_wtemp2);


}//end of if(p->mode !=3)

     //cudaThreadSynchronize();
     

//checkErrors_i("memory allocation");

     //copy data back to cpu so we can compute and update the grid (on the cpu)

 
    //cudaMemcpy(*w, *d_w, NVAR*dimp*sizeof(real), cudaMemcpyDeviceToHost);

/*if((*p)->mode==3)
{
  
  int ii[3];
  ii[0]=0;
  ii[1]=0;
  ii[2]=0;
 init_user_i(*wmod,*wd,*p,ii);
}*/







    //cudaMemcpy(*w, *d_w, NVAR*dimp*sizeof(real), cudaMemcpyDeviceToHost);
    //setup the grid and dx values here


    //cudaMemcpy(*d_w, *w, NVAR*dimp*sizeof(real), cudaMemcpyHostToDevice);


 //initialise grid on the gpu
 //we currently don't do this to avoid use of additional memory on GPU
 /*for(dir=0; dir<NDIM; dir++)
 {
     zerotempv_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w, *d_wnew, *d_wmod, *d_dwn1,  *d_wd, *d_wtemp, *d_wtemp1, *d_wtemp2,dir);
     cudaThreadSynchronize();     
     gridsetup_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w, *d_wnew, *d_wmod, *d_dwn1,  *d_wd, *d_wtemp, *d_wtemp1, *d_wtemp2,dir);
     cudaThreadSynchronize();
     setupdx_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w, *d_wnew, *d_wmod, *d_dwn1,  *d_wd, *d_wtemp, *d_wtemp1, *d_wtemp2,dir);
     cudaThreadSynchronize();
  }*/

	    printf("called initialiser\n");
	//cudaMemcpy(*w, *d_w, NVAR*dimp*sizeof(real), cudaMemcpyDeviceToHost);
if((*p)->mode !=3)
{
	cudaMemcpy(*state, *d_state, sizeof(struct state), cudaMemcpyDeviceToHost);
        cudaMemcpy(*p, *d_p, sizeof(struct params), cudaMemcpyDeviceToHost);
}


//checkErrors_i("memory allocation");checkErrors_i("memory allocation");



	//cudaMemcpy(*wnew, *d_wnew, NVAR*((*p)->n[0])* ((*p)->n[1])*sizeof(real), cudaMemcpyDeviceToHost);
	//cudaMemcpy(*b, *d_b, (((*p)->n[0])* ((*p)->n[1]))*sizeof(real), cudaMemcpyDeviceToHost);

        // printf("mod times step %f %f\n",(*p)->dt, ((*wnew)[10+16*((*p)->n[0])+((*p)->n[0])*((*p)->n[1])*b1]));



  return 0;



}




int cuupdatemod(struct params **p, struct bparams **bp,real **w, real **wnew, real **wd, struct state **state, struct params **d_p, struct bparams **d_bp,real **d_w, real **d_wnew, real **d_wmod, real **d_dwn1, real **d_wd, struct state **d_state, real **d_wtemp, real **d_wtemp1, real **d_wtemp2)
{
  int deviceCount;
  int dir;
 
  printf("in cuinit\n");
 // real *adb;
  real *adw, *adwnew;
  struct params *adp;
  struct bparams *adbp;
  struct state *ads;

 
 
  int dimp=(((*p)->n[0]))*(((*p)->n[1]));

   
 #ifdef USE_SAC_3D   
  dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif  

   int numBlocks = (dimp+numThreadsPerBlock-1) / numThreadsPerBlock;
   
    printf("calling updatemod\n");
     //init_parallel(struct params *p, real *b, real *u, real *v, real *h)
    // init_parallel<<<dimGrid,dimBlock>>>(*d_p,*d_b,*d_u,*d_v,*d_h);
    // init_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w, *d_wnew, *d_b);
     updatemod_parallel<<<numBlocks, numThreadsPerBlock>>>(*d_p,*d_w, *d_wnew, *d_wmod, *d_dwn1,  *d_wd, *d_wtemp, *d_wtemp1, *d_wtemp2);
     //cudaThreadSynchronize();
     


 


  return 0;
}




/*! Cartesian or polar grid. Determine x at the boundaries.
! Determine often needed combinations of x, such as dx or dvolume.
! Determine variables for axial symmetry
!
! ixe          - edge coordinate of the grid touching the boundary region
! ixf          - coordinate inside of ixe
! qx           - x with an extended index range for calculation of dx   */

int initgrid(struct params **p,   struct state **state, real **wd, struct params **d_p,  real **d_dwn1, real **d_wd, struct state **d_state, real **d_wtemp, real **d_wtemp1, real **d_wtemp2)
{
    real *ttemp2;
    int ii[NDIM];
    int ii1[3],ii2[3],ix;
    int ip,jp,kp,kpo;
    int dir,dir1,dir2;
    int ixmin,ixmax,ixe,ixf;
    real *wda=*wd;
    //real *wa=*wmod;
 int dimp=(((*p)->n[0]))*(((*p)->n[1]));

/*if(((*p)->ipe)==2)
      {
checkErrors_i("initgrid memory allocation");
}*/
 #ifdef USE_SAC_3D
 
   dimp=(((*p)->n[0]))*(((*p)->n[1]))*(((*p)->n[2]));
#endif      
    kp=0;
    //printf("called initgrid %d\n",(*p)->ipe);
    

    for(int i=0;i<3;i++)
    {
       ii1[i]=0;
       ii2[i]=0;
    }
    #ifdef USE_SAC
    ttemp2=(real *) malloc( (NTEMP2+2)*(((*p)->n[0])+2)* (((*p)->n[1])+2)*sizeof(real));
    #endif
    #ifdef USE_SAC_3D
    ttemp2=(real *)malloc((NTEMP2+2)*(((*p)->n[0])+2)* (((*p)->n[1])+2)* (((*p)->n[2])+2)*sizeof(real));
    #endif
    
     //cudaMemcpy(*wmod, *d_wmod, NVAR*dimp*sizeof(real), cudaMemcpyDeviceToHost);
     cudaMemcpy(*wd, *d_wd, NDERV*dimp*sizeof(real), cudaMemcpyDeviceToHost);
     for(dir=0;dir<NDIM;dir++)
     for(ii[0]=0; ii[0]<((*p)->n[0])+2; ii[0]++)
     for(ii[1]=0; ii[1]<((*p)->n[1])+2; ii[1]++)
     		     #ifdef USE_SAC_3D
                   for(ii[2]=0; ii[2]<((*p)->n[2])+2; ii[2]++)
                 #endif
                 {
                        ip=ii[0];
                        jp=ii[1];
         		     #ifdef USE_SAC_3D
                       kp=ii[2];
                     #endif                   
                       
	    switch(dir)
	    {
		     case 0:
	                ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui)]= 0;
		     break;
	
		     case 1:
			 ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui1)]= 0;
		     break;
		         
		     #ifdef USE_SAC_3D
		     case 2:
			 ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui2)]= 0;
		     break;			
		     #endif
	     }
      }	
 

     kp=1;
     kpo=0;
     for(dir=0;dir<NDIM;dir++)
        for(ii[0]=1; ii[0]<((*p)->n[0])+1; ii[0]++)
           for(ii[1]=1; ii[1]<((*p)->n[1])+1; ii[1]++)
		#ifdef USE_SAC_3D
		   for(ii[2]=1; ii[2]<((*p)->n[2])+1; ii[2]++)
		#endif
                {
                        ip=ii[0];
                        jp=ii[1];
         		     #ifdef USE_SAC_3D
                       kp=ii[2];
                       kpo=kp;
                     #endif                   
                       
	    switch(dir)
	    {
		     case 0:
	                ttemp2[encode3p2_i(*p,ip,jp,kpo,tmpnui)]= (wda[encode3_i(*p,ip-1,jp-1,kp-1,pos1)]);
		     break;
	
		     case 1:
			 ttemp2[encode3p2_i(*p,ip,jp,kpo,tmpnui1)]= (wda[(encode3_i(*p,ip-1,jp-1,kp-1,pos2))]);
		     break;
		         
		     #ifdef USE_SAC_3D
		     case 2:
			 ttemp2[encode3p2_i(*p,ip,jp,kpo,tmpnui2)]= (wda[(encode3_i(*p,ip-1,jp-1,kp-1,pos3))]);
		     break;			
		     #endif
	     }


      }	


   /* if((*p)->ipe==3   ) 
        for(ii[1]=1; ii[1]<((*p)->n[1])+1; ii[1]++)                                                            
        for(ii[0]=1; ii[0]<((*p)->n[0])+1; ii[0]++)    
           {
                      ip=ii[0];
                        jp=ii[1];
                         printf("ii0, ii1 %d %d %16.20f %16.20f\n",ip,jp, ttemp2[encode3p2_i(*p,ip,jp,kpo,tmpnui)],ttemp2[(encode3p2_i(*p,ip,jp,kpo,tmpnui1))]);

            }*/


  	
   	//update grid edges
     kp=0;
     for(dir=0;dir<NDIM;dir++)
     {
                
                       
	    switch(dir)
	    {
		     case 0:
                       ixmax=((*p)->n[0])+1;//ixGmax1+1; 
                       ixmin=((*p)->n[0])-1;//ixmin1=ixGmax1+1                      

                      #ifdef USE_MULTIGPU
			if(((*p)->fullgridini)==1    ||   ((*p)->mpiupperb[dir])==1) ixmin=((*p)->n[0])+1;//ixGmax1+1;
                      #endif

                       ixe=ixmin-1; 
                       ixf=ixe-1;


                       //upper layers
			     for(dir1=0;dir1<NDIM;dir1++)
			     {
				     for(ii[0]=ixmin; ii[0]<=ixmax; ii[0]++)
				     for(ii[1]=0; ii[1]<((*p)->n[1])+2; ii[1]++)
				     		 #ifdef USE_SAC_3D
						   for(ii[2]=0; ii[2]<((*p)->n[2])+2; ii[2]++)
						 #endif
						 {
				                        ix=ii[0];
                                                        ip=ii[0];
							jp=ii[1];
					 		     #ifdef USE_SAC_3D
						       kp=ii[2];
						     #endif  
                                                       for(dir2=0;dir2<NDIM;dir2++)
                                                       {
                                                         ii1[dir2]=ii[dir2];
                                                         ii2[dir2]=ii[dir2];
                                                       }
                                                       ii1[0]=ixe;
                                                       ii2[0]=ixf; 

 

                                                       ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui+dir1)]=(1+abs(ixe-ix))* (ttemp2[(fencode3p2_i(*p,ii1,tmpnui+dir1))])-(abs(ixe-ix))* (ttemp2[(fencode3p2_i(*p,ii2,tmpnui+dir1))]);
						      //ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui+dir1)]=(1+abs(ixe-ix))* (wda[fencode3_i(*p,ii1,pos1+dir1)]);
						  }

				}



                      //lower layers

                       ixmin=0;//ixmin1=ixGmin1-1;
                       ixmax=2;//ixmax1=ixGmin1-1 

                     #ifdef USE_MULTIGPU
			if(((*p)->fullgridini)==1    ||  ((*p)->mpilowerb[dir])==1) ixmax=0;
                      #endif

                  
                       ixe=ixmax+1; 
                       ixf=ixe+1;

			     for(dir1=0;dir1<NDIM;dir1++)
			     {
				     for(ii[0]=ixmin; ii[0]<=ixmax; ii[0]++)
				     for(ii[1]=0; ii[1]<((*p)->n[1])+2; ii[1]++)
				     		 #ifdef USE_SAC_3D
						   for(ii[2]=0; ii[2]<((*p)->n[2])+2; ii[2]++)
						 #endif
						 {
							ix=ip=ii[0];
							jp=ii[1];
					 		     #ifdef USE_SAC_3D
						       kp=ii[2];
						     #endif  
                                                       for(dir2=0;dir2<NDIM;dir2++)
                                                       {
                                                         ii1[dir2]=ii[dir2];
                                                         ii2[dir2]=ii[dir2];
                                                       }
                                                       ii1[0]=ixe;
                                                       ii2[0]=ixf;




 
    ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui+dir1)]=(1+abs(ixe-ix))* (ttemp2[(fencode3p2_i(*p,ii1,tmpnui+dir1))])-(abs(ixe-ix))* (ttemp2[(fencode3p2_i(*p,ii2,tmpnui+dir1))]);
// ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui+dir1)]= (ttemp2[(fencode3p2_i(*p,ii1,tmpnui+dir1))])+ (ttemp2[(fencode3p2_i(*p,ii2,tmpnui+dir1))]);
   // qx(ix,ixmin2:ixmax2,jdim)=(1+abs(ixe-ix))*qx(ixe,ixmin2:ixmax2,jdim)- abs(ixe-ix) *qx(ixf,ixmin2:ixmax2,jdim)


//		if((*p)->ipe==0   && ii[1]==0) 
//                         printf("ii0, ii1 %d %d %16.20f %16.20f %d %d %d %d %d  %d %d %d %d\n",ip,jp, ttemp2[encode3p2_i(*p,ip,jp,kpo,tmpnui)],ttemp2[(encode3p2_i(*p,ip,jp,kpo,tmpnui1))],dir1,ixe,ixf,ixmin,ixmax,ii1[0],ii1[1],ii2[0],ii2[1]);





						  }

				}
		     break;
	
		     case 1:
                       ixmax=((*p)->n[1])+1;//ixGmax1+1; 
                       ixmin=((*p)->n[1])-1;//ixmin1=ixGmax1+1                      

                      #ifdef USE_MULTIGPU
			if(((*p)->fullgridini)==1    ||  ((*p)->mpiupperb[dir])==1) ixmin=((*p)->n[1])+1;//ixGmax1+1;
                      #endif
                     
                       ixe=ixmin-1; 
                       ixf=ixe-1;


                       //upper layers
			     for(dir1=0;dir1<NDIM;dir1++)
			     {
                 for(ii[0]=0; ii[0]<((*p)->n[0])+2; ii[0]++)
				     for(ii[1]=ixmin; ii[1]<=ixmax; ii[1]++)
				     
				     		 #ifdef USE_SAC_3D
						   for(ii[2]=0; ii[2]<((*p)->n[2])+2; ii[2]++)
						 #endif
						 {
							ip=ii[0];
							ix=jp=ii[1];
					 		     #ifdef USE_SAC_3D
						       kp=ii[2];
						     #endif  
                                                       for(dir2=0;dir2<NDIM;dir2++)
                                                       {
                                                         ii1[dir2]=ii[dir2];
                                                         ii2[dir2]=ii[dir2];
                                                       }
                                                       ii1[1]=ixe;
                                                       ii2[1]=ixf; 





						       ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui+dir1)]=(1+abs(ixe-ix))* (ttemp2[(fencode3p2_i(*p,ii1,tmpnui+dir1))])-(abs(ixe-ix))* (ttemp2[(fencode3p2_i(*p,ii2,tmpnui+dir1))]);
						      //ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui+dir1)]=(1+abs(ixe-ix))* (wda[fencode3_i(*p,ii1,pos1+dir1)]);
						  }

				}



                      //lower layers

                       ixmin=0;//ixmin1=ixGmin1-1;
                       ixmax=2;//ixmax1=ixGmin1-1 

                     #ifdef USE_MULTIGPU
			if(((*p)->fullgridini)==1    ||  ((*p)->mpilowerb[dir])==1) ixmax=0;
                      #endif
                
                       ixe=ixmax+1; 
                       ixf=ixe+1;

			     for(dir1=0;dir1<NDIM;dir1++)
			     {
			         for(ii[0]=0; ii[0]<((*p)->n[0])+2; ii[0]++)	
				     for(ii[1]=ixmin; ii[1]<=ixmax; ii[1]++)
				     		 #ifdef USE_SAC_3D
						   for(ii[2]=0; ii[2]<((*p)->n[2])+2; ii[2]++)
						 #endif
						 {
							ip=ii[0];
							ix=jp=ii[1];
					 		     #ifdef USE_SAC_3D
						       kp=ii[2];
						     #endif  
                                                       for(dir2=0;dir2<NDIM;dir2++)
                                                       {
                                                         ii1[dir2]=ii[dir2];
                                                         ii2[dir2]=ii[dir2];
                                                       }
                                                       ii1[1]=ixe;
                                                       ii2[1]=ixf; 



                                                        

						       ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui+dir1)]=(1+abs(ixe-ix))* (ttemp2[(fencode3p2_i(*p,ii1,tmpnui+dir1))])-(abs(ixe-ix))* (ttemp2[(fencode3p2_i(*p,ii2,tmpnui+dir1))]);
//write(*,*) jdim,ixe,ixf,ix,ixmin1,ixmax1,ixmin2,ixmax2, qx(ixmin1:ixmax1,ix,jdim),qx(ixmin1:ixmax1,&
//                   ixe,jdim),qx(ixmin1:ixmax1,ixf,jdim)
//if((*p)->ipe==0   && ii[0]==0)
//                                                             printf("ixe, ix %d %d %d %d %d %d %16.20f %16.20f %16.20f\n",dir1,ixe,ixf,ix,ixmin,ixmax, ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui+dir1)],ttemp2[(fencode3p2_i(*p,ii1,tmpnui+dir1))],ttemp2[(fencode3p2_i(*p,ii2,tmpnui+dir1))]);
//if((*p)->ipe==0   && ii[0]==0)
//                                                             printf("ixe, ix %d %d %d %d %d %d %16.20f %16.20f %16.20f\n",dir1,ixe,ixf,ix,ixmin,ixmax, ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui+dir1)],ttemp2[(fencode3p2_i(*p,ii1,tmpnui+dir1))],ttemp2[(fencode3p2_i(*p,ii2,tmpnui+dir1))]);


						  }

				}




		     break;
		         
		     #ifdef USE_SAC_3D
		     case 2:
                       ixmax=((*p)->n[2])+1;//ixGmax1+1; 
                       ixmin=((*p)->n[2])-1;//ixmin1=ixGmax1+1                      

                      #ifdef USE_MULTIGPU
			if(((*p)->mpiupperb[dir])==1) ixmin=((*p)->n[2])+1;//ixGmax1+1;
                      #endif
                  
                       ixe=ixmin-1; 
                       ixf=ixe-1;


                       //upper layers
			     for(dir1=0;dir1<NDIM;dir1++)
			     {
                 for(ii[0]=0; ii[0]<((*p)->n[0])+2; ii[0]++)
                 for(ii[1]=0; ii[1]<((*p)->n[1])+2; ii[1]++)
				     
				     		 #ifdef USE_SAC_3D
						  
			        for(ii[2]=ixmin; ii[2]<=ixmax; ii[2]++)
						 #endif
						 {
							ip=ii[0];
							jp=ii[1];
					 		     #ifdef USE_SAC_3D
						       ix=kp=ii[2];
						     #endif  
                                                       for(dir2=0;dir2<NDIM;dir2++)
                                                       {
                                                         ii1[dir2]=ii[dir2];
                                                         ii2[dir2]=ii[dir2];
                                                       }
                                                       ii1[2]=ixe;
                                                       ii2[2]=ixf; 
						       ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui+dir1)]=(1+abs(ixe-ix))* (ttemp2[(fencode3p2_i(*p,ii1,tmpnui+dir1))])-(abs(ixe-ix))* (ttemp2[(fencode3p2_i(*p,ii2,tmpnui+dir1))]);
						      //ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui+dir1)]=(1+abs(ixe-ix))* (wda[fencode3_i(*p,ii1,pos1+dir1)]);
						  }

				}



                      //lower layers
                      //lower layers

                       ixmin=0;//ixmin1=ixGmin1-1;
                       ixmax=2;//ixmax1=ixGmin1-1 

                     #ifdef USE_MULTIGPU
			if(((*p)->fullgridini)==1    ||  ((*p)->mpilowerb[dir])==1) ixmax=0;
                      #endif
                   
                       ixe=ixmax+1; 
                       ixf=ixe+1;

			     for(dir1=0;dir1<NDIM;dir1++)
			     {
			         for(ii[0]=0; ii[0]<((*p)->n[0])+2; ii[0]++)
                     for(ii[1]=0; ii[1]<((*p)->n[1])+2; ii[1]++)	
				     
				     		 #ifdef USE_SAC_3D
						   
						    for(ii[2]=ixmin; ii[2]<=ixmax; ii[2]++)
						 #endif
						 {
							ip=ii[0];
							jp=ii[1];
					 		     #ifdef USE_SAC_3D
						       ix=kp=ii[2];
						     #endif  
                                                       for(dir2=0;dir2<NDIM;dir2++)
                                                       {
                                                         ii1[dir2]=ii[dir2];
                                                         ii2[dir2]=ii[dir2];
                                                       }
                                                       ii1[2]=ixe;
                                                       ii2[2]=ixf; 
						       ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui+dir1)]=(1+abs(ixe-ix))* (ttemp2[(fencode3p2_i(*p,ii1,tmpnui+dir1))])-(abs(ixe-ix))* (ttemp2[(fencode3p2_i(*p,ii2,tmpnui+dir1))]);
						  }

				}



		     break;			
		     #endif
	     }
      }	


kp=0;

     for(dir=0;dir<NDIM;dir++)
        for(ii[0]=0; ii[0]<((*p)->n[0]); ii[0]++)
           for(ii[1]=0; ii[1]<((*p)->n[1]); ii[1]++)
		#ifdef USE_SAC_3D
		   for(ii[2]=0; ii[2]<((*p)->n[2]); ii[2]++)
		#endif
                {
                        ip=ii[0]+1;
                        jp=ii[1]+1;
         		     #ifdef USE_SAC_3D
                       kp=ii[2]+1;
                     #endif                   
                       
	    switch(dir)
	    {
		     case 0:
	                 (wda[fencode3_i(*p,ii,pos1)])=ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui)];
                      //  if(ip==1)
                      //  printf("delx 0 %d %d %16.20f  %16.20f \n",ii[0],ii[1],wda[(encode3_i(*p,ip-1,jp-1,kp,delx1))],wda[(encode3_i(*p,ip-1,jp-1,kp,delx2))]);
		     break;
	
		     case 1:
			  (wda[(fencode3_i(*p,ii,pos2))])=ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui1)];
//if(ip==1)
                 //       printf("delx 1 %d %d %16.20f  %16.20f \n",ii[0],ii[1],wda[(encode3_i(*p,ip-1,jp-1,kp,delx1))],wda[(encode3_i(*p,ip-1,jp-1,kp,delx2))]);

		     break;
		         
		     #ifdef USE_SAC_3D
		     case 2:
			  (wda[(fencode3_i(*p,ii,pos3))])=ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui2)];
		     break;			
		     #endif
	     }
      }	



   	//calculate dx
  kp=0;
  kpo=0;

                   for(dir=0;dir<NDIM;dir++)
                 {

  for(ii[0]=1; ii[0]<((*p)->n[0])+1; ii[0]++)
     for(ii[1]=1; ii[1]<((*p)->n[1])+1; ii[1]++)
  //for(ii[0]=0; ii[0]<((*p)->n[0]); ii[0]++)
  //   for(ii[1]=0; ii[1]<((*p)->n[1]); ii[1]++)

     		     #ifdef USE_SAC_3D
                   for(ii[2]=1; ii[2]<((*p)->n[2])+1; ii[2]++)
                 #endif
{

                        ip=ii[0];
                        jp=ii[1];
         		     #ifdef USE_SAC_3D
                       
                       kp=ii[2];
                        kpo=kp-1;
                     #endif                   
                       
	    switch(dir)
	    {
		     case 0:
	               // (wda[(encode3_i(*p,ip-1,jp-1,kpo,delx1))])=/*(*p)->dx[0];//*/0.5*(ttemp2[encode3p2_i(*p,ip+1,jp,kp,tmpnui)]-ttemp2[encode3p2_i(*p,ip-1,jp,kp,tmpnui)]);
                  (wda[(encode3_i(*p,ip-1,jp-1,kpo,delx1))])=/*(*p)->dx[0];//*/0.5*(ttemp2[encode3p2_i(*p,ip+1,jp,kp,tmpnui)]-ttemp2[encode3p2_i(*p,ip-1,jp,kp,tmpnui)]);

		if(wda[(encode3_i(*p,ip-1,jp-1,kpo,delx1))]==0) wda[(encode3_i(*p,ip-1,jp-1,kpo,delx1))]=(*p)->dx[0];
	              //  if(ip==128  && jp==128 && kp==128)
                      //  printf("delx 0 %d %d %d %16.20f  %16.20f   %16.20f \n",ii[0]-1,ii[1]-1,ii[2]-1,wda[(encode3_i(*p,ip-1,jp-1,kp-1,delx1))],wda[(encode3_i(*p,ip-1,jp-1,kp-1,delx2))],wda[(encode3_i(*p,ip-1,jp-1,kp-1,delx3))]);
		     break;
	
		     case 1:
			(wda[(encode3_i(*p,ip-1,jp-1,kpo,delx2))])=/*(*p)->dx[1];//*/0.5*(ttemp2[encode3p2_i(*p,ip,jp+1,kp,tmpnui1)]-ttemp2[encode3p2_i(*p,ip,jp-1,kp,tmpnui1)]);
		if(wda[(encode3_i(*p,ip-1,jp-1,kpo,delx2))]==0) wda[(encode3_i(*p,ip-1,jp-1,kpo,delx2))]=(*p)->dx[1];
	               // if(ip==128  && jp==128 && kp==128)
                       //   printf("delx 1 %d %d %d %16.20f  %16.20f   %16.20f \n",ii[0]-1,ii[1]-1,ii[2]-1,wda[(encode3_i(*p,ip-1,jp-1,kp-1,delx1))],wda[(encode3_i(*p,ip-1,jp-1,kp-1,delx2))],wda[(encode3_i(*p,ip-1,jp-1,kp-1,delx3))]);

		        //printf("delx2 %d %d %g ",ii[0],ii[1],wda[(fencode3_i(*p,ii,delx2))]);
		     break;
		         
		     #ifdef USE_SAC_3D
		     case 2:
			(wda[(encode3_i(*p,ip-1,jp-1,kpo,delx3))])=0.5*(ttemp2[encode3p2_i(*p,ip,jp,kp+1,tmpnui2)]-ttemp2[encode3p2_i(*p,ip,jp,kp-1,tmpnui2)]);
		if(wda[(encode3_i(*p,ip-1,jp-1,kpo,delx3))]==0) wda[(encode3_i(*p,ip-1,jp-1,kpo,delx3))]=(*p)->dx[2];
	              //  if(ip==128  && jp==128 && kp==128)
                      //  printf("delx 2 %d %d %d %16.20f  %16.20f   %16.20f \n",ii[0]-1,ii[1]-1,ii[2]-1,wda[(encode3_i(*p,ip-1,jp-1,kp-1,delx1))],wda[(encode3_i(*p,ip-1,jp-1,kp-1,delx2))],wda[(encode3_i(*p,ip-1,jp-1,kp-1,delx3))]);

		     break;			
		     #endif
	     }
      }
  printf("\n");
}


printf("dx=%g dy=%g\n",(*p)->dx[0], (*p)->dx[1] );




     kp=0;

   // if((*p)->ipe==3)
    // for(dir=0;dir<NDIM;dir++)
//for(ii[1]=0; ii[1]<((*p)->n[1])+2; ii[1]++)
        //for(ii[0]=0; ii[0]<((*p)->n[0])+2; ii[0]++)
           
            // {

              //          ip=ii[0];
              //          jp=ii[1];
                //if(ii[0]==0)
              //  printf("delx %d %d %16.20f  %16.20f  %16.20f  %16.20f \n",ii[0],ii[1],wda[(fencode3_i(*p,ii,pos1))],wda[(fencode3_i(*p,ii,pos2))],wda[(fencode3_i(*p,ii,delx1))],wda[(fencode3_i(*p,ii,delx2))]);
//printf("ttemp2 %d %d %16.20f  %16.20f  \n",ii[0],ii[1],ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui)],ttemp2[encode3p2_i(*p,ip,jp,kp,tmpnui1)]);
            //  }

  //  cudaMemcpy(*d_w, *w, NVAR*dimp*sizeof(real), cudaMemcpyHostToDevice);
  //  cudaMemcpy(*d_wd, *wd, NDERV*dimp*sizeof(real), cudaMemcpyHostToDevice);

 free(ttemp2);

 

    //cudaMemcpy(*d_wmod, *wmod, NVAR*dimp*sizeof(real), cudaMemcpyHostToDevice);
    cudaMemcpy(*d_wd, *wd, NDERV*dimp*sizeof(real), cudaMemcpyHostToDevice);


    //  
     



  

   
  return 0;



}


