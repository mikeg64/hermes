
__device__ __host__
real usersource2a_MODID(real *dw, real *wd, real *w, struct params *p,int *ii, int dir) {

  real ddc=0;
  real fi, fim1;
  real  fip2=0, fim2=0;
 // real ddc1;
  real ddcx=0,ddcy=0;

   real source=0;


        #if defined USE_SAC  || defined USE_SAC_3D
//computept3_MODID(w,wd,p,ii);



  /*    		source= -wd[fencode3_MODID(p,ii,ptb)]*grad3d_MODID(wd,p,ii,vel1+dir,dir);
               source += +w[fencode3_MODID(p,ii,b1b)]*w[fencode3_MODID(p,ii,b1b+dir)]*grad3d_MODID(wd,p,ii,vel1,0)+w[fencode3_MODID(p,ii,b2b)]*w[fencode3_MODID(p,ii,b1b+dir)]*grad3d_MODID(wd,p,ii,vel1+1,1);*/
         #endif


        #if defined USE_SAC_3D
            //   source += +w[fencode3_MODID(p,ii,b3b)]*w[fencode3_MODID(p,ii,b1b+dir)]*grad3d_MODID(wd,p,ii,vel3,2);
        #endif

  return source;


 
}

__device__ __host__
real usersource1a_MODID(real *dw, real *wd, real *w, struct params *p,int *ii, int dir) {

  real ddc=0;
  real fi, fim1;
  real  fip2=0, fim2=0;
 // real ddc1;
  real ddcx=0,ddcy=0;

   real source=0;


        #if defined USE_SAC  || defined USE_SAC_3D
//computept3_MODID(w,wd,p,ii);



      	/*	source= -wd[fencode3_MODID(p,ii,ptb)]*grad3d_MODID(wd,p,ii,vel1+dir,dir);
               source += +w[fencode3_MODID(p,ii,b1b)]*w[fencode3_MODID(p,ii,b1b+dir)]*grad3d_MODID(wd,p,ii,vel1,0)+w[fencode3_MODID(p,ii,b2b)]*w[fencode3_MODID(p,ii,b1b+dir)]*grad3d_MODID(wd,p,ii,vel1+1,1);*/
         #endif


        #if defined USE_SAC_3D
              // source += +w[fencode3_MODID(p,ii,b3b)]*w[fencode3_MODID(p,ii,b1b+dir)]*grad3d_MODID(wd,p,ii,vel3,2);
        #endif

  return source;


}



__device__ __host__
int addsourceterms2_MODID(real *dw, real *wd, real *w, struct params *p, struct state *s,int *ii,int field,int dir) {

  int direction;
  int status=0;
  real divsource=0;
  //dw[fencode3_MODID(p,ii,field)]= grad_MODID(wd,p,ii,source,dir);//+grad_MODID(wd,p,ii,f2,1); 


 #if defined USE_SAC  ||  defined USE_SAC_3D

  
 /* if(field==energy)
  {    
     computept3_MODID(w,wd,p,ii);
     dw[fencode3_MODID(p,ii,field)]=usersource2a_MODID(dw, wd, w, p,ii,dir)+w[fencode3_MODID(p,ii,rho)]*((p->g[dir])*w[fencode3_MODID(p,ii,mom1+dir)]    )/(w[fencode3_MODID(p,ii,rho)]+w[fencode3_MODID(p,ii,rhob)]);
   }*/


 #endif
  return ( status);
}

__device__ __host__
int addsourceterms1_MODID(real *dw, real *wd, real *w, struct params *p, struct state *s,int *ii,int field,int dir) {

  int direction;
  int status=0;
  real divsource=0;
  //dw[fencode3_MODID(p,ii,field)]= grad_MODID(wd,p,ii,source,dir);//+grad_MODID(wd,p,ii,f2,1); 


 #if defined USE_SAC  ||  defined USE_SAC_3D

  
  /*if(field==energy)
  {    
     computept3_MODID(w,wd,p,ii);
     dw[fencode3_MODID(p,ii,field)]=usersource2a_MODID(dw, wd, w, p,ii,dir)+w[fencode3_MODID(p,ii,rho)]*((p->g[dir])*w[fencode3_MODID(p,ii,mom1+dir)]    )/(w[fencode3_MODID(p,ii,rho)]+w[fencode3_MODID(p,ii,rhob)]);
   }*/


 #endif
  return ( status);
}

