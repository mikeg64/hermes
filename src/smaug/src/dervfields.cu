/*Operators for derived fields*/


//These operators for the derived fields are routines which may be used by the kernel functions
//They are used as follows during make the field MODID is replaced by a unique identifier
//for the particular cuda source file 
//For example the file centdiff1.cu has identifier cd1
//so that computej3_MODID becomes computej3_cd1

//The make routine copies the resulting file to a new file called dervfields_cd1.cuh
//This file is then included using the line #include "../include/dervfields_cd1.cuh"
//in centdiff1.cu

//The routines in centdiff1.cu must call these routines with _MODID replaced by _cd1



__device__ __host__
void computej3_MODID(real *wmod,real *wd,struct params *p,int *ii)
{


 /* wd[fencode3_MODID(p,ii,current1)]=(grad3d_MODID(wmod,p,ii,b3,1))/(p->mu);
  wd[fencode3_MODID(p,ii,current2)]=(grad3d_MODID(wmod,p,ii,b3,0))/(p->mu);
  wd[fencode3_MODID(p,ii,current3)]=(grad3d_MODID(wmod,p,ii,b2,0)-grad3d_MODID(wmod,p,ii,b1,1))/(p->mu);*/
  
          #ifdef USE_SAC
	 /* wd[fencode3_MODID(p,ii,current1)]+=(grad3d_MODID(wmod,p,ii,b3b,1))/(p->mu);
	  wd[fencode3_MODID(p,ii,current2)]+=(grad3d_MODID(wmod,p,ii,b3b,0))/(p->mu);
	  wd[fencode3_MODID(p,ii,current3)]+=(grad3d_MODID(wmod,p,ii,b2b,0)-grad3d_MODID(wmod,p,ii,b1b,1))/(p->mu);*/
         #endif

          #ifdef USE_SAC_3D

         /* wd[fencode3_MODID(p,ii,current1)]-=(  (grad3d_MODID(wmod,p,ii,b2b,2))+ (grad3d_MODID(wmod,p,ii,b2,2)) )/(p->mu)
          wd[fencode3_MODID(p,ii,current2)]+=(  (grad3d_MODID(wmod,p,ii,b1b,2))+ (grad3d_MODID(wmod,p,ii,b1,2)) )/(p->mu)*/

	 /* wd[fencode3_MODID(p,ii,current1)]+=(grad3d_MODID(wmod,p,ii,b3b,1))/(p->mu);
	  wd[fencode3_MODID(p,ii,current2)]+=(grad3d_MODID(wmod,p,ii,b3b,0))/(p->mu);
	  wd[fencode3_MODID(p,ii,current3)]+=(grad3d_MODID(wmod,p,ii,b2b,0)-grad3d_MODID(wmod,p,ii,b1b,1))/(p->mu);*/
         #endif

}

__device__ __host__
void computebdotv3_MODID(real *wmod,real *wd,struct params *p,int *ii)
{
        #ifdef USE_SAC


wd[fencode3_MODID(p,ii,bdotv)]=((wmod[fencode3_MODID(p,ii,b1)]+wmod[fencode3_MODID(p,ii,b1b)])*wmod[fencode3_MODID(p,ii,mom1)]+(wmod[fencode3_MODID(p,ii,b2)]+wmod[fencode3_MODID(p,ii,b2b)])*wmod[fencode3_MODID(p,ii,mom2)])/(wmod[fencode3_MODID(p,ii,rho)]+wmod[fencode3_MODID(p,ii,rhob)]);

         #endif
        #ifdef USE_SAC_3D


wd[fencode3_MODID(p,ii,bdotv)]=((wmod[fencode3_MODID(p,ii,b1)]+wmod[fencode3_MODID(p,ii,b1b)])*wmod[fencode3_MODID(p,ii,mom1)]+(wmod[fencode3_MODID(p,ii,b2)]+wmod[fencode3_MODID(p,ii,b2b)])*wmod[fencode3_MODID(p,ii,mom2)]+(wmod[fencode3_MODID(p,ii,b3)]+wmod[fencode3_MODID(p,ii,b3b)])*wmod[fencode3_MODID(p,ii,mom3)])/(wmod[fencode3_MODID(p,ii,rho)]+wmod[fencode3_MODID(p,ii,rhob)]);

         #endif
 // return ( status);
}

__device__ __host__
void computedivb3_MODID(real *wmod,real *wd,struct params *p,int *ii)
{
      //#ifdef USE_SAC

		wd[fencode3_MODID(p,ii,divb)]=grad3d_MODID(wmod,p,ii,b1,0)+grad3d_MODID(wmod,p,ii,b2,1);
		wd[fencode3_MODID(p,ii,divb)]+=grad3d_MODID(wmod,p,ii,b1b,0)+grad3d_MODID(wmod,p,ii,b2b,1);
        //#endif
        #ifdef USE_SAC_3D
		wd[fencode3_MODID(p,ii,divb)]=grad3d_MODID(wmod,p,ii,b2,2);
		wd[fencode3_MODID(p,ii,divb)]+=grad3d_MODID(wmod,p,ii,b2b,2);		
         #endif
 // return ( status);
}

__device__ __host__
void computevel3_MODID(real *wmod,real *wd,struct params *p,int *ii)
{

        #ifdef USE_SAC_3D
		wd[fencode3_MODID(p,ii,vel1)]=wmod[fencode3_MODID(p,ii,mom1)]/(wmod[fencode3_MODID(p,ii,rho)]+wmod[fencode3_MODID(p,ii,rhob)]);
                wd[fencode3_MODID(p,ii,vel2)]=wmod[fencode3_MODID(p,ii,mom2)]/(wmod[fencode3_MODID(p,ii,rho)]+wmod[fencode3_MODID(p,ii,rhob)]);
                wd[fencode3_MODID(p,ii,vel3)]=wmod[fencode3_MODID(p,ii,mom3)]/(wmod[fencode3_MODID(p,ii,rho)]+wmod[fencode3_MODID(p,ii,rhob)]);
        #endif

        #ifdef USE_SAC
		wd[fencode3_MODID(p,ii,vel1)]=wmod[fencode3_MODID(p,ii,mom1)]/(wmod[fencode3_MODID(p,ii,rho)]+wmod[fencode3_MODID(p,ii,rhob)]);
                wd[fencode3_MODID(p,ii,vel2)]=wmod[fencode3_MODID(p,ii,mom2)]/(wmod[fencode3_MODID(p,ii,rho)]+wmod[fencode3_MODID(p,ii,rhob)]);
        #endif
       #ifdef ADIABHYDRO
		wd[fencode3_MODID(p,ii,vel1)]=wmod[fencode3_MODID(p,ii,mom1)]/(wmod[fencode3_MODID(p,ii,rho)]);
                wd[fencode3_MODID(p,ii,vel2)]=wmod[fencode3_MODID(p,ii,mom2)]/(wmod[fencode3_MODID(p,ii,rho)]);

         #endif
 // return ( status);
}







__device__ __host__
void computept3_MODID(real *wmod,real *wd,struct params *p,int *ii)
{
 // int status=0;

#ifdef ADIABHYDRO

/*below used for adiabatic hydrodynamics*/
 wd[fencode3_MODID(p,ii,pressuret)]=(p->adiab)*pow(wmod[fencode3_MODID(p,ii,rho)],p->gamma);

#elif defined(USE_SAC)
 //wmod[fencode3_MODID(p,ii,b1b)]=0;
// wmod[fencode3_MODID(p,ii,b2b)]=0;
wd[fencode3_MODID(p,ii,pressuret)]=((p->gamma)-2.0)*((wmod[fencode3_MODID(p,ii,b1)]*wmod[fencode3_MODID(p,ii,b1b)]+wmod[fencode3_MODID(p,ii,b2)]*wmod[fencode3_MODID(p,ii,b2b)])+0.5*(wmod[fencode3_MODID(p,ii,b1)]*wmod[fencode3_MODID(p,ii,b1)]+wmod[fencode3_MODID(p,ii,b2)]*wmod[fencode3_MODID(p,ii,b2)]));
wd[fencode3_MODID(p,ii,pressuret)]=((p->gamma)-1.0)*( wmod[fencode3_MODID(p,ii,energy)]-0.5*(wmod[fencode3_MODID(p,ii,mom1)]*wmod[fencode3_MODID(p,ii,mom1)]+wmod[fencode3_MODID(p,ii,mom2)]*wmod[fencode3_MODID(p,ii,mom2)])/(wmod[fencode3_MODID(p,ii,rho)]+wmod[fencode3_MODID(p,ii,rhob)]))-wd[fencode3_MODID(p,ii,pressuret)];



#elif defined(USE_SAC_3D)

wd[fencode3_MODID(p,ii,pressuret)]=((p->gamma)-2.0)*((wmod[fencode3_MODID(p,ii,b1)]*wmod[fencode3_MODID(p,ii,b1b)]+wmod[fencode3_MODID(p,ii,b2)]*wmod[fencode3_MODID(p,ii,b2b)]+wmod[fencode3_MODID(p,ii,b3)]*wmod[fencode3_MODID(p,ii,b3b)])+0.5*(wmod[fencode3_MODID(p,ii,b1)]*wmod[fencode3_MODID(p,ii,b1)]+wmod[fencode3_MODID(p,ii,b2)]*wmod[fencode3_MODID(p,ii,b2)]+wmod[fencode3_MODID(p,ii,b3)]*wmod[fencode3_MODID(p,ii,b3)]));




//wd[fencode3_MODID(p,ii,pressuret)]=0.0;
//wd[fencode3_MODID(p,ii,pressuret)]=((p->gamma)-1.0)*( wmod[fencode3_MODID(p,ii,energy)]);
//wd[fencode3_MODID(p,ii,pressuret)]=((p->gamma)-1.0)*( wmod[fencode3_MODID(p,ii,energy)]-0.5*(wmod[fencode3_MODID(p,ii,mom1)]*wmod[fencode3_MODID(p,ii,mom1)]+wmod[fencode3_MODID(p,ii,mom2)]*wmod[fencode3_MODID(p,ii,mom2)]+wmod[fencode3_MODID(p,ii,mom3)]*wmod[fencode3_MODID(p,ii,mom3)])/(wmod[fencode3_MODID(p,ii,rho)]+wmod[fencode3_MODID(p,ii,rhob)]))-wd[fencode3_MODID(p,ii,pressuret)];



wd[fencode3_MODID(p,ii,pressuret)]=((p->gamma)-1.0)*(wmod[fencode3_MODID(p,ii,energy)] -0.5*(wmod[fencode3_MODID(p,ii,mom1)]*wmod[fencode3_MODID(p,ii,mom1)]+wmod[fencode3_MODID(p,ii,mom2)]*wmod[fencode3_MODID(p,ii,mom2)]+wmod[fencode3_MODID(p,ii,mom3)]*wmod[fencode3_MODID(p,ii,mom3)])/(wmod[fencode3_MODID(p,ii,rho)]+wmod[fencode3_MODID(p,ii,rhob)]))-wd[fencode3_MODID(p,ii,pressuret)];

#else

wd[fencode3_MODID(p,ii,pressuret)]=  ((p->gamma)-1.0)*wmod[fencode3_MODID(p,ii,energy)]+(1.0-0.5*(p->gamma))*(wmod[fencode3_MODID(p,ii,b1)]*wmod[fencode3_MODID(p,ii,b1)]+wmod[fencode3_MODID(p,ii,b2)]*wmod[fencode3_MODID(p,ii,b2)])+0.5*(1.0-(p->gamma))*(wmod[fencode3_MODID(p,ii,mom1)]*wmod[fencode3_MODID(p,ii,mom1)]+wmod[fencode3_MODID(p,ii,mom2)]*wmod[fencode3_MODID(p,ii,mom2)])/wmod[fencode3_MODID(p,ii,rho)];

#endif



  //if(wd[fencode3_MODID(p,ii,pressuret)]<0)
              //wd[fencode3_MODID(p,ii,pressuret)]=1.0e-10;
	//      wd[fencode3_MODID(p,ii,pressuret)]=0.01;


 // return ( status);
}


__device__ __host__
void computepbg3_MODID(real *wmod,real *wd,struct params *p,int *ii)
{
 // int status=0;

#ifdef ADIABHYDRO

/*below used for adiabatic hydrodynamics*/
;// wd[fencode3_MODID(p,ii,pressuret)]=(p->adiab)*pow(wmod[fencode3_MODID(p,ii,rho)],p->gamma);

#elif defined(USE_SAC)
 //wmod[fencode3_MODID(p,ii,b1b)]=0;
// wmod[fencode3_MODID(p,ii,b2b)]=0;



 wd[fencode3_MODID(p,ii,ptb)]=  ((p->gamma)-1)*wmod[fencode3_MODID(p,ii,energyb)]- 0.5*((p->gamma)-2)*(wmod[fencode3_MODID(p,ii,b1b)]*wmod[fencode3_MODID(p,ii,b1b)]+wmod[fencode3_MODID(p,ii,b2b)]*wmod[fencode3_MODID(p,ii,b2b)]) ;

#elif defined(USE_SAC_3D)



 wd[fencode3_MODID(p,ii,ptb)]=  ((p->gamma)-1)*wmod[fencode3_MODID(p,ii,energyb)]- 0.5*((p->gamma)-2)*(wmod[fencode3_MODID(p,ii,b1b)]*wmod[fencode3_MODID(p,ii,b1b)]+wmod[fencode3_MODID(p,ii,b2b)]*wmod[fencode3_MODID(p,ii,b2b)]+wmod[fencode3_MODID(p,ii,b3b)]*wmod[fencode3_MODID(p,ii,b3b)]) ;


#endif



  //if(wd[fencode3_MODID(p,ii,pressuret)]<0)
              //wd[fencode3_MODID(p,ii,pressuret)]=1.0e-10;
	//      wd[fencode3_MODID(p,ii,pressuret)]=0.01;


 // return ( status);
}





__device__ __host__
void computepk3_MODID(real *wmod,real *wd,struct params *p,int *ii)
{
  //int status=0;

#ifdef ADIABHYDRO

/*below used for adiabatic hydrodynamics*/
wd[fencode3_MODID(p,ii,pressurek)]=(p->adiab)*pow(wmod[fencode3_MODID(p,ii,rho)],p->gamma);
wd[fencode3_MODID(p,ii,vel1)]=wmod[fencode3_MODID(p,ii,mom1)]/(wmod[fencode3_MODID(p,ii,rho)]);
wd[fencode3_MODID(p,ii,vel2)]=wmod[fencode3_MODID(p,ii,mom2)]/(wmod[fencode3_MODID(p,ii,rho)]);

#elif defined(USE_SAC)
 wd[fencode3_MODID(p,ii,pressurek)]= 0.5*((wmod[fencode3_MODID(p,ii,mom1)]*wmod[fencode3_MODID(p,ii,mom1)]+wmod[fencode3_MODID(p,ii,mom2)]*wmod[fencode3_MODID(p,ii,mom2)])/(wmod[fencode3_MODID(p,ii,rho)]+wmod[fencode3_MODID(p,ii,rhob)]));

wd[fencode3_MODID(p,ii,pressurek)]=wd[fencode3_MODID(p,ii,pressurek)]+(0.5*(wmod[fencode3_MODID(p,ii,b1)]*wmod[fencode3_MODID(p,ii,b1)]+wmod[fencode3_MODID(p,ii,b2)]*wmod[fencode3_MODID(p,ii,b2)]) +(wmod[fencode3_MODID(p,ii,b1)]*wmod[fencode3_MODID(p,ii,b1b)]+wmod[fencode3_MODID(p,ii,b2)]*wmod[fencode3_MODID(p,ii,b2b)]) );

wd[fencode3_MODID(p,ii,pressurek)]=((p->gamma)-1)*wmod[fencode3_MODID(p,ii,energy)]-wd[fencode3_MODID(p,ii,pressurek)];


#elif defined(USE_SAC_3D)

 wd[fencode3_MODID(p,ii,pressurek)]= 0.5*((wmod[fencode3_MODID(p,ii,mom1)]*wmod[fencode3_MODID(p,ii,mom1)]+wmod[fencode3_MODID(p,ii,mom2)]*wmod[fencode3_MODID(p,ii,mom2)]+wmod[fencode3_MODID(p,ii,mom3)]*wmod[fencode3_MODID(p,ii,mom3)])/(wmod[fencode3_MODID(p,ii,rho)]+wmod[fencode3_MODID(p,ii,rhob)]));

wd[fencode3_MODID(p,ii,pressurek)]=wd[fencode3_MODID(p,ii,pressurek)]+(0.5*(wmod[fencode3_MODID(p,ii,b1)]*wmod[fencode3_MODID(p,ii,b1)]+wmod[fencode3_MODID(p,ii,b2)]*wmod[fencode3_MODID(p,ii,b2)]+wmod[fencode3_MODID(p,ii,b3)]*wmod[fencode3_MODID(p,ii,b3)]) +(wmod[fencode3_MODID(p,ii,b1)]*wmod[fencode3_MODID(p,ii,b1b)]+wmod[fencode3_MODID(p,ii,b2)]*wmod[fencode3_MODID(p,ii,b2b)]+wmod[fencode3_MODID(p,ii,b3)]*wmod[fencode3_MODID(p,ii,b3b)]) );

wd[fencode3_MODID(p,ii,pressurek)]=((p->gamma)-1)*wmod[fencode3_MODID(p,ii,energy)]-wd[fencode3_MODID(p,ii,pressurek)];


#else

 wd[fencode3_MODID(p,ii,pressurek)]=((p->gamma)-1)*(wmod[fencode3_MODID(p,ii,energy)]- 0.5*(wmod[fencode3_MODID(p,ii,mom1)]*wmod[fencode3_MODID(p,ii,mom1)]+wmod[fencode3_MODID(p,ii,mom2)]*wmod[fencode3_MODID(p,ii,mom2)])/wmod[fencode3_MODID(p,ii,rho)]-0.5*(wmod[fencode3_MODID(p,ii,b1)]*wmod[fencode3_MODID(p,ii,b1)]+wmod[fencode3_MODID(p,ii,b2)]*wmod[fencode3_MODID(p,ii,b2)]) );


#endif



  //if(wd[fencode3_MODID(p,ii,pressurek)]<0)
  //           wd[fencode3_MODID(p,ii,pressurek)]=0.001;
  //return ( status);
}

__device__ __host__
void computec3_MODID(real *wmod,real *wd,struct params *p,int *ii,int dir)
{

 real cfasti,pk; 
#ifdef ADIABHYDRO
/*below used for adiabatic hydrodynamics*/
  wd[fencode3_MODID(p,ii,soundspeed)]=sqrt((p->adiab)/wmod[fencode3_MODID(p,ii,rho)]);
#elif defined(USE_SAC)

pk=((p->gamma)-1)*(wmod[fencode3_MODID(p,ii,energy)]
- 0.5*((wmod[fencode3_MODID(p,ii,mom1)]*wmod[fencode3_MODID(p,ii,mom1)]+wmod[fencode3_MODID(p,ii,mom2)]*wmod[fencode3_MODID(p,ii,mom2)])/(wmod[fencode3_MODID(p,ii,rho)]+wmod[fencode3_MODID(p,ii,rhob)]))-0.5*(wmod[fencode3_MODID(p,ii,b1)]*wmod[fencode3_MODID(p,ii,b1)]+wmod[fencode3_MODID(p,ii,b2)]*wmod[fencode3_MODID(p,ii,b2)]) -(wmod[fencode3_MODID(p,ii,b1b)]*wmod[fencode3_MODID(p,ii,b1)]+wmod[fencode3_MODID(p,ii,b2b)]*wmod[fencode3_MODID(p,ii,b2)]) );

wd[fencode3_MODID(p,ii,soundspeed)]=(((p->gamma))
*(pk+(((p->gamma))-1)*(
wmod[fencode3_MODID(p,ii,energyb)] -0.5*(wmod[fencode3_MODID(p,ii,b1b)]*wmod[fencode3_MODID(p,ii,b1b)]+wmod[fencode3_MODID(p,ii,b2b)]*wmod[fencode3_MODID(p,ii,b2b)])))
/(wmod[fencode3_MODID(p,ii,rho)]+wmod[fencode3_MODID(p,ii,rhob)]));


wd[fencode3_MODID(p,ii,cfast)]=( ((wmod[fencode3_MODID(p,ii,b1)]*wmod[fencode3_MODID(p,ii,b1)]+wmod[fencode3_MODID(p,ii,b2)]*wmod[fencode3_MODID(p,ii,b2)]) + (wmod[fencode3_MODID(p,ii,b1b)]*wmod[fencode3_MODID(p,ii,b1b)]+wmod[fencode3_MODID(p,ii,b2b)]*wmod[fencode3_MODID(p,ii,b2b)]) +2.0*(wmod[fencode3_MODID(p,ii,b1b)]*wmod[fencode3_MODID(p,ii,b1)]+wmod[fencode3_MODID(p,ii,b2b)]*wmod[fencode3_MODID(p,ii,b2)]))/(wmod[fencode3_MODID(p,ii,rho)]+wmod[fencode3_MODID(p,ii,rhob)]))+(wd[fencode3_MODID(p,ii,soundspeed)]);

cfasti=0.5*(
wd[fencode3_MODID(p,ii,cfast)]
+sqrt(wd[fencode3_MODID(p,ii,cfast)]*wd[fencode3_MODID(p,ii,cfast)]
-4.0*wd[fencode3_MODID(p,ii,soundspeed)]*((wmod[fencode3_MODID(p,ii,b1b+dir)]+wmod[fencode3_MODID(p,ii,b1+dir)])*(wmod[fencode3_MODID(p,ii,b1b+dir)]+wmod[fencode3_MODID(p,ii,b1+dir)]))
/(wmod[fencode3_MODID(p,ii,rho)]+wmod[fencode3_MODID(p,ii,rhob)])));

wd[fencode3_MODID(p,ii,cfast)]=sqrt(cfasti)+sacdabs_MODID(wmod[fencode3_MODID(p,ii,mom1+dir)]/(wmod[fencode3_MODID(p,ii,rho)]+wmod[fencode3_MODID(p,ii,rhob)]));
wd[fencode3_MODID(p,ii,soundspeed)]=sqrt(wd[fencode3_MODID(p,ii,soundspeed)]);
//wd[fencode3_MODID(p,ii,cfast)]=( ((wmod[fencode3_MODID(p,ii,b1)]*wmod[fencode3_MODID(p,ii,b1)]+wmod[fencode3_MODID(p,ii,b2)]*wmod[fencode3_MODID(p,ii,b2)]) + (wmod[fencode3_MODID(p,ii,b1b)]*wmod[fencode3_MODID(p,ii,b1b)]+wmod[fencode3_MODID(p,ii,b2b)]*wmod[fencode3_MODID(p,ii,b2b)]) +2.0*(wmod[fencode3_MODID(p,ii,b1b)]*wmod[fencode3_MODID(p,ii,b1)]+wmod[fencode3_MODID(p,ii,b2b)]*wmod[fencode3_MODID(p,ii,b2)]))/(wmod[fencode3_MODID(p,ii,rho)]+wmod[fencode3_MODID(p,ii,rhob)]));

//wd[fencode3_MODID(p,ii,cfast)]=cfasti;

#elif defined(USE_SAC_3D)


pk=((p->gamma)-1)*(wmod[fencode3_MODID(p,ii,energy)]
- 0.5*((wmod[fencode3_MODID(p,ii,mom1)]*wmod[fencode3_MODID(p,ii,mom1)]+wmod[fencode3_MODID(p,ii,mom2)]*wmod[fencode3_MODID(p,ii,mom2)]+wmod[fencode3_MODID(p,ii,mom3)]*wmod[fencode3_MODID(p,ii,mom3)])/(wmod[fencode3_MODID(p,ii,rho)]+wmod[fencode3_MODID(p,ii,rhob)]))-0.5*(wmod[fencode3_MODID(p,ii,b1)]*wmod[fencode3_MODID(p,ii,b1)]+wmod[fencode3_MODID(p,ii,b2)]*wmod[fencode3_MODID(p,ii,b2)]+wmod[fencode3_MODID(p,ii,b3)]*wmod[fencode3_MODID(p,ii,b3)]) -(wmod[fencode3_MODID(p,ii,b1b)]*wmod[fencode3_MODID(p,ii,b1)]+wmod[fencode3_MODID(p,ii,b2b)]*wmod[fencode3_MODID(p,ii,b2)]+wmod[fencode3_MODID(p,ii,b3b)]*wmod[fencode3_MODID(p,ii,b3)]) );

wd[fencode3_MODID(p,ii,soundspeed)]=(((p->gamma))
*(pk+(((p->gamma))-1)*(
wmod[fencode3_MODID(p,ii,energyb)] -0.5*(wmod[fencode3_MODID(p,ii,b1b)]*wmod[fencode3_MODID(p,ii,b1b)]+wmod[fencode3_MODID(p,ii,b2b)]*wmod[fencode3_MODID(p,ii,b2b)]+wmod[fencode3_MODID(p,ii,b3b)]*wmod[fencode3_MODID(p,ii,b3b)])))
/(wmod[fencode3_MODID(p,ii,rho)]+wmod[fencode3_MODID(p,ii,rhob)]));


wd[fencode3_MODID(p,ii,cfast)]=( ((wmod[fencode3_MODID(p,ii,b1)]*wmod[fencode3_MODID(p,ii,b1)]+wmod[fencode3_MODID(p,ii,b2)]*wmod[fencode3_MODID(p,ii,b2)]+wmod[fencode3_MODID(p,ii,b3)]*wmod[fencode3_MODID(p,ii,b3)]) + (wmod[fencode3_MODID(p,ii,b1b)]*wmod[fencode3_MODID(p,ii,b1b)]+wmod[fencode3_MODID(p,ii,b2b)]*wmod[fencode3_MODID(p,ii,b2b)]+wmod[fencode3_MODID(p,ii,b3b)]*wmod[fencode3_MODID(p,ii,b3b)]) +2.0*(wmod[fencode3_MODID(p,ii,b1b)]*wmod[fencode3_MODID(p,ii,b1)]+wmod[fencode3_MODID(p,ii,b2b)]*wmod[fencode3_MODID(p,ii,b2)]+wmod[fencode3_MODID(p,ii,b3b)]*wmod[fencode3_MODID(p,ii,b3)]))/(wmod[fencode3_MODID(p,ii,rho)]+wmod[fencode3_MODID(p,ii,rhob)]))+(wd[fencode3_MODID(p,ii,soundspeed)]);

cfasti=0.5*(
wd[fencode3_MODID(p,ii,cfast)]
+sqrt(wd[fencode3_MODID(p,ii,cfast)]*wd[fencode3_MODID(p,ii,cfast)]
-4.0*wd[fencode3_MODID(p,ii,soundspeed)]*((wmod[fencode3_MODID(p,ii,b1b+dir)]+wmod[fencode3_MODID(p,ii,b1+dir)])*(wmod[fencode3_MODID(p,ii,b1b+dir)]+wmod[fencode3_MODID(p,ii,b1+dir)]))
/(wmod[fencode3_MODID(p,ii,rho)]+wmod[fencode3_MODID(p,ii,rhob)])));

wd[fencode3_MODID(p,ii,cfast)]=sqrt(cfasti)+sacdabs_MODID(wmod[fencode3_MODID(p,ii,mom1+dir)]/(wmod[fencode3_MODID(p,ii,rho)]+wmod[fencode3_MODID(p,ii,rhob)]));
wd[fencode3_MODID(p,ii,soundspeed)]=sqrt(wd[fencode3_MODID(p,ii,soundspeed)]);



#else
wd[fencode3_MODID(p,ii,soundspeed)]=sqrt(((p->gamma))*wd[fencode3_MODID(p,ii,pressuret)]/wmod[fencode3_MODID(p,ii,rho)]);


wd[fencode3_MODID(p,ii,cfast)]=sqrt(((wmod[fencode3_MODID(p,ii,b1)]*wmod[fencode3_MODID(p,ii,b1)]+wmod[fencode3_MODID(p,ii,b2)]*wmod[fencode3_MODID(p,ii,b2)])/wmod[fencode3_MODID(p,ii,rho)])+(wd[fencode3_MODID(p,ii,soundspeed)]*wd[fencode3_MODID(p,ii,soundspeed)]));

#endif



  
}
//uptohere so far thursday  24th march
__device__ __host__
void computecmax3_MODID(real *wmod,real *wd,struct params *p,int *ii)
{
 //p->cmax=0.02;
#ifdef ADIABHYDRO
       if(wd[fencode3_MODID(p,ii,soundspeed)]>(p->cmax))
                    // atomicExch(&(p->cmax),(wd[fencode3_MODID(p,ii,soundspeed)]));
                    p->cmax=(wd[fencode3_MODID(p,ii,soundspeed)]);
#else
       if(wd[fencode3_MODID(p,ii,soundspeed)]>(p->cmax))
                    // atomicExch(&(p->cmax),(wd[fencode3_MODID(p,ii,soundspeed)]));
                    p->cmax=(wd[fencode3_MODID(p,ii,soundspeed)]);
       if(wd[fencode3_MODID(p,ii,cfast)]>(p->cmax))
                     //atomicExch(&(p->cmax),(wd[fencode3_MODID(p,ii,cfast)]));
                    p->cmax=(wd[fencode3_MODID(p,ii,cfast)]);
#endif

}


