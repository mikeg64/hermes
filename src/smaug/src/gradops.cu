/* differential operators and boundary condition*/


//These operators for the differential operators
//and the boundary condition routines and may be used by the kernel functions
//They are used as follows during make the field MODID is replaced by a unique identifier
//for the particular cuda source file 
//For example the file centdiff1.cu has identifier cd1
//so that dimproduct_MODID becomes dimproduct_cd1

//The make routine copies the resulting file to a new file called gradops_cd1.cuh
//This file is then included using the line #include "../include/gradops_cd1.cuh"
//in centdiff1.cu

//The routines in centdiff1.cu must call these routines with _MODID replaced by _cd1





__device__ __host__
int dimproduct_MODID (struct params *dp) {

  int tot=1;
  for(int i=0;i<NDIM;i++)
    tot*=dp->n[i];
  return tot; 
}






__device__ __host__
int fencode_MODID (struct params *dp,int ix, int iy, int field) {


    return ( (iy * ((dp)->n[0]) + ix)+(field*((dp)->n[0])*((dp)->n[1])));

}


__device__ __host__
int fencode3_MODID (struct params *dp,int *ii, int field) {


#ifdef USE_SAC_3D
   return (ii[2]*((dp)->n[0])*((dp)->n[1])  + ii[1] * ((dp)->n[0]) + ii[0]+(field*((dp)->n[0])*((dp)->n[1])*((dp)->n[2])));
#else
   return ( ii[1] * ((dp)->n[0]) + ii[0]+(field*((dp)->n[0])*((dp)->n[1])));
#endif

}

__device__ __host__
int encode3p1_MODID (struct params *dp,int ix, int iy, int iz, int field) {


  #ifdef USE_SAC_3D
    return ( (iz*(((dp)->n[0])+1)*(((dp)->n[1])+1)  + iy * (((dp)->n[0])+1) + ix)+(field*(((dp)->n[0])+1)*(((dp)->n[1])+1)*(((dp)->n[2])+1)));
  #else
    return ( (iy * (((dp)->n[0])+1) + ix)+(field*(((dp)->n[0])+1)*(((dp)->n[1])+1)));
  #endif
}




__device__ __host__
int encode3p2_MODID (struct params *dp,int ix, int iy, int iz, int field) {


  #ifdef USE_SAC_3D
    return ( (iz*(((dp)->n[0])+2)*(((dp)->n[1])+2)  + iy * (((dp)->n[0])+2) + ix)+(field*(((dp)->n[0])+2)*(((dp)->n[1])+2)*(((dp)->n[2])+2)));
  #else
    return ( (iy * (((dp)->n[0])+2) + ix)+(field*(((dp)->n[0])+2)*(((dp)->n[1])+2)));
  #endif
}

__device__ __host__
int fencode3p2_MODID (struct params *dp,int *ii, int field) {

  return(encode3p2_MODID(dp,ii[0],ii[1],ii[2],field));
}


__device__ __host__
int encode3_MODID (struct params *dp,int ix, int iy, int iz, int field) {


  #ifdef USE_SAC_3D
    return ( (iz*((dp)->n[0])*((dp)->n[1])  + iy * ((dp)->n[0]) + ix)+(field*((dp)->n[0])*((dp)->n[1])*((dp)->n[2])));
  #else
    return ( (iy * ((dp)->n[0]) + ix)+(field*((dp)->n[0])*((dp)->n[1])));
  #endif
}

__device__ __host__
int encodefixed13_MODID (struct params *dp,int ix, int iy, int iz, int field) {
  #ifdef USE_SAC_3D
    return ( (ix*((dp)->n[1])*((dp)->n[1])  + iy * ((dp)->n[2]) + iz)+(field*4*((dp)->n[1])*((dp)->n[2])));
  #else
    return ( (ix * ((dp)->n[1]) + iy)+(field*4*((dp)->n[1])));
  #endif
}

__device__ __host__
int encodefixed23_MODID (struct params *dp,int ix, int iy, int iz, int field) {
  #ifdef USE_SAC_3D
    return ( (iy*((dp)->n[0])*((dp)->n[2])  + ix * ((dp)->n[0]) + iz)+(4*field*((dp)->n[0])*((dp)->n[2])));
  #else
    return ( (  iy * ((dp)->n[0]) + ix)+(4*field*((dp)->n[0])));
  #endif
}

__device__ __host__
int encodefixed33_MODID (struct params *dp,int ix, int iy, int iz, int field) {
  #ifdef USE_SAC_3D
    return ( ( iz*((dp)->n[0])*((dp)->n[1])  + iy * ((dp)->n[0]) + ix)+(4*field*((dp)->n[0])*((dp)->n[1])));
  #endif
}








__device__ __host__
real grad3d_MODID(real *wmod,struct params *p,int *ii,int field,int dir)
{


 real grad=0;

 
 

 switch(dir)
 {
   case 0:
 
#ifdef USE_SAC_3D
  #ifdef USE_DORDER3
 if(ii[0]>2 && ii[0]<((p->n[0])-3) )
  grad=(  ( ((3*wmod[encode3_MODID(p,ii[0]+1,ii[1],ii[2],field)]-3*wmod[encode3_MODID(p,ii[0]-1,ii[1],ii[2],field)]+3.0*(wmod[encode3_MODID(p,ii[0]-2,ii[1],ii[2],field)]-wmod[encode3_MODID(p,ii[0]+2,ii[1],ii[2],field)])/5.0-(wmod[encode3_MODID(p,ii[0]-3,ii[1],ii[2],field)]-wmod[encode3_MODID(p,ii[0]+3,ii[1],ii[2],field)])/15.0)/2.0))/(2.0*(p->dx[0]))    );
 else 
  #endif
if(ii[0]>1 && ii[0]<((p->n[0])-2) )
 grad=(  ( ((8*wmod[encode3_MODID(p,ii[0]+1,ii[1],ii[2],field)]-8*wmod[encode3_MODID(p,ii[0]-1,ii[1],ii[2],field)]+wmod[encode3_MODID(p,ii[0]-2,ii[1],ii[2],field)]-wmod[encode3_MODID(p,ii[0]+2,ii[1],ii[2],field)])/6.0))/(2.0*(p->dx[0]))    );

   if((ii[0]==(p->n[0])-3) || (ii[0]==(p->n[0])-4)  && ii[1]>1   && ii[1]<(p->n[1])-2 && ii[2]>1   && ii[2]<(p->n[2])-2  )
       grad=0;
   else if(ii[0]==2 || ii[0]==3  && ii[1]>1   && ii[1]<(p->n[1])-2 && ii[2]>1   && ii[2]<(p->n[2])-2  )
       grad=0;
#else

  #ifdef USE_DORDER3
if(ii[0]>2 && ii[0]<((p->n[0])-3) )
 grad=(  ( ((3*wmod[encode3_MODID(p,ii[0]+1,ii[1],0,field)]-3*wmod[encode3_MODID(p,ii[0]-1,ii[1],0,field)]+3.0*(wmod[encode3_MODID(p,ii[0]-2,ii[1],0,field)]-wmod[encode3_MODID(p,ii[0]+2,ii[1],0,field)])/5.0-(wmod[encode3_MODID(p,ii[0]-3,ii[1],0,field)]-wmod[encode3_MODID(p,ii[0]+3,ii[1],0,field)])/15.0)/2.0))/(2.0*(p->dx[0]))    );
 else 
  #endif
if(ii[0]>1 && ii[0]<((p->n[0])-2) )
 grad=(  ( ((8*wmod[encode3_MODID(p,ii[0]+1,ii[1],0,field)]-8*wmod[encode3_MODID(p,ii[0]-1,ii[1],0,field)]+wmod[encode3_MODID(p,ii[0]-2,ii[1],0,field)]-wmod[encode3_MODID(p,ii[0]+2,ii[1],0,field)])/6.0))/(2.0*(p->dx[0]))    );

   if((ii[0]==(p->n[0])-3) || (ii[0]==(p->n[0])-4)  && ii[1]>1   && ii[1]<(p->n[1])-2  )
       grad=0;
   else if(ii[0]==2 || ii[0]==3  && ii[1]>1   && ii[1]<(p->n[1])-2  )
       grad=0;
#endif



   break;

   case 1:

#ifdef USE_SAC_3D

  #ifdef USE_DORDER3
 if(ii[1]>2 && ii[1]<((p->n[1])-3) )
  grad=(  ( ((3*wmod[encode3_MODID(p,ii[0],ii[1]+1,ii[2],field)]-3*wmod[encode3_MODID(p,ii[0],ii[1]-1,ii[2],field)]+3.0*(wmod[encode3_MODID(p,ii[0],ii[1]-2,ii[2],field)]-wmod[encode3_MODID(p,ii[0],ii[1]+2,ii[2],field)])/5.0-(wmod[encode3_MODID(p,ii[0],ii[1]-3,ii[2],field)]-wmod[encode3_MODID(p,ii[0],ii[1]+3,ii[2],field)])/15.0)/2.0))/(2.0*(p->dx[1]))    );
 else 
#endif
if( ii[1] >1 &&  ii[1]<((p->n[1])-2))
	grad=(  ( ((8*wmod[encode3_MODID(p,ii[0],ii[1]+1,ii[2],field)]-8*wmod[encode3_MODID(p,ii[0],ii[1]-1,ii[2],field)]+wmod[encode3_MODID(p,ii[0],ii[1]-2,ii[2],field)]-wmod[encode3_MODID(p,ii[0],ii[1]+2,ii[2],field)])/6.0))/(2.0*(p->dx[1]))    );

   if((ii[1]==(p->n[1])-3) || (ii[1]==(p->n[1])-4)  && ii[0]>1   && ii[0]<(p->n[0])-2  && ii[2]>1   && ii[2]<(p->n[2])-2  )
       grad=0;
   else if(ii[1]==2 || ii[1]==3  && ii[0]>1   && ii[0]<(p->n[0])-2  && ii[2]>1   && ii[2]<(p->n[2])-2  )
       grad=0;
#else

  #ifdef USE_DORDER3
if(ii[1]>2 && ii[1]<((p->n[1])-3) )
 grad=(  ( ((3*wmod[encode3_MODID(p,ii[0],ii[1]+1,0,field)]-3*wmod[encode3_MODID(p,ii[0],ii[1]-1,0,field)]+3.0*(wmod[encode3_MODID(p,ii[0],ii[1]-2,0,field)]-wmod[encode3_MODID(p,ii[0],ii[1]+2,0,field)])/5.0-(wmod[encode3_MODID(p,ii[0],ii[1]-3,0,field)]-wmod[encode3_MODID(p,ii[0],ii[1]+3,0,field)])/15.0)/2.0))/(2.0*(p->dx[1]))    );
else  
#endif
if( ii[1] >1 &&  ii[1]<((p->n[1])-2))
	grad=(  ( ((8*wmod[encode3_MODID(p,ii[0],ii[1]+1,0,field)]-8*wmod[encode3_MODID(p,ii[0],ii[1]-1,0,field)]+wmod[encode3_MODID(p,ii[0],ii[1]-2,0,field)]-wmod[encode3_MODID(p,ii[0],ii[1]+2,0,field)])/6.0))/(2.0*(p->dx[1]))    );

   if((ii[1]==(p->n[1])-3) || (ii[1]==(p->n[1])-4)  && ii[0]>1   && ii[0]<(p->n[0])-2  )
       grad=0;
   else if(ii[1]==2 || ii[1]==3  && ii[0]>1   && ii[0]<(p->n[0])-2  )
       grad=0;
#endif
   break;


   case 2:

#ifdef USE_SAC_3D
  #ifdef USE_DORDER3
 if(ii[2]>2 && ii[2]<((p->n[2])-3) )
  grad=(  ( ((3*wmod[encode3_MODID(p,ii[0],ii[1],ii[2]+1,field)]-3*wmod[encode3_MODID(p,ii[0],ii[1],ii[2]-1,field)]+3.0*(wmod[encode3_MODID(p,ii[0],ii[1],ii[2]-2,field)]-wmod[encode3_MODID(p,ii[0],ii[1],ii[2]+2,field)])/5.0-(wmod[encode3_MODID(p,ii[0],ii[1],ii[2]-3,field)]-wmod[encode3_MODID(p,ii[0],ii[1],ii[2]+3,field)])/15.0)/2.0))/(2.0*(p->dx[2]))    );
 else 
#endif
if( ii[2] >1 &&  ii[2]<((p->n[2])-2))
	grad=(  ( ((8*wmod[encode3_MODID(p,ii[0],ii[1],ii[2]+1,field)]-8*wmod[encode3_MODID(p,ii[0],ii[1],ii[2]-1,field)]+wmod[encode3_MODID(p,ii[0],ii[1],ii[2]-2,field)]-wmod[encode3_MODID(p,ii[0],ii[1],ii[2]+2,field)])/6.0))/(2.0*(p->dx[2]))    );

   if((ii[2]==(p->n[2])-3) || (ii[2]==(p->n[2])-4)  && ii[0]>1   && ii[0]<(p->n[0])-2 && ii[1]>1   && ii[1]<(p->n[1])-2  )
       grad=0;
   else if(ii[2]==2 || ii[2]==3  && ii[0]>1   && ii[0]<(p->n[0])-2 && ii[1]>1   && ii[1]<(p->n[1])-2  )
       grad=0;
#endif
   break;

}



 return grad;


}

__device__ __host__
real grad1l3_MODID(real *wmod,struct params *p,int *ii,int field,int dir)
{
 real grad=0;
   int i,j,k;
   i=ii[0];
   j=ii[1];
   k=0;
   #ifdef USE_SAC_3D
    k=ii[2];
   #endif


 if((dir == 0) && i>0 && i<((p->n[0])))
 {
    grad=(  ( wmod[encode3_MODID(p,i,j,k,field)]-wmod[encode3_MODID(p,i-1,j,k,field)]) /((p->dx[0]))    );

   #ifdef USE_SAC_3D
	   if((i==(p->n[0])-2) || (i==(p->n[0])-3)  && j>0   && j<(p->n[1])-1 && k>1   && k<(p->n[2])-1 )
	       grad=0;
	   else if(i==1 || i==2  && j>0   && j<(p->n[1])-1 && k>0   && k<(p->n[2])-1 )
	       grad=0;
   #else
	   if((i==(p->n[0])-2) || (i==(p->n[0])-3)  && j>0   && j<(p->n[1])-1  )
	       grad=0;
	   else if(i==1 || i==2  && j>0   && j<(p->n[1])-1  )
	       grad=0;
   #endif
 }
 else if((dir == 1)    && j>0 && j<((p->n[1])))
 {
    grad=(  ( wmod[encode3_MODID(p,i,j,k,field)]-wmod[encode3_MODID(p,i,j-1,k,field)])/((p->dx[1]))    );
   #ifdef USE_SAC_3D
	   if((j==(p->n[1])-2) || (j==(p->n[1])-3)  && i>0   && i<(p->n[0])-1 && k>0   && k<(p->n[2])-1 )
	       grad=0;
	   else if(j==1 || j==2  && i>0   && i<(p->n[0])-1 && k>0   && k<(p->n[2])-1 )
	       grad=0;
   #else
	   if((j==(p->n[1])-2) || (j==(p->n[1])-3)  && i>0   && i<(p->n[0])-1  )
	       grad=0;
	   else if(j==1 || j==2  && i>0   && i<(p->n[0])-1  )
	       grad=0;
   #endif


  }
   #ifdef USE_SAC_3D
 else if((dir == 2)    && k>0 && k<((p->n[2])))
 {
    grad=(  ( wmod[encode3_MODID(p,i,j,k,field)]-wmod[encode3_MODID(p,i,j,k-1,field)])/((p->dx[2]))    );

   if((k==(p->n[2])-2) || (k==(p->n[2])-3)  && i>0   && i<(p->n[0])-1  && j>0   && j<(p->n[1])-1 )
       grad=0;
   else if(k==1 || k==2  && i>0   && i<(p->n[0])-1  && j>0   && j<(p->n[1])-1 )
       grad=0;


  }
  #endif
 return grad;

}

__device__ __host__
real grad1r3_MODID(real *wmod,struct params *p,int *ii,int field,int dir)
{
  real grad=0;
   int i,j,k;
   i=ii[0];
   j=ii[1];
   k=0;
   #ifdef USE_SAC_3D
    k=ii[2];
   #endif


 if((dir == 0) && /*i>0 &&*/ i<((p->n[0])-1))
 {
    grad=(  ( wmod[encode3_MODID(p,i+1,j,k,field)]-wmod[encode3_MODID(p,i,j,k,field)]) /((p->dx[0]))    );

   #ifdef USE_SAC_3D
	   if((i==(p->n[0])-2) || (i==(p->n[0])-3)  && j>0   && j<(p->n[1])-1 && k>0   && k<(p->n[2])-1 )
	       grad=0;
	   else if(i==1 || i==2  && j>0   && j<(p->n[1])-1 && k>0   && k<(p->n[2])-1 )
	       grad=0;
   #else
	   if((i==(p->n[0])-2) || (i==(p->n[0])-3)  && j>0   && j<(p->n[1])-1  )
	       grad=0;
	   else if(i==1 || i==2  && j>0   && j<(p->n[1])-1  )
	       grad=0;
   #endif
 }
 else if((dir == 1)    /*&& j>0*/ && j<((p->n[1])-1))
 {
    grad=(  ( wmod[encode3_MODID(p,i,j+1,k,field)]-wmod[encode3_MODID(p,i,j,k,field)])/((p->dx[1]))    );
   #ifdef USE_SAC_3D
	   if((j==(p->n[1])-2) || (j==(p->n[1])-3)  && i>0   && i<(p->n[0])-1 && k>0   && k<(p->n[2])-1 )
	       grad=0;
	   else if(j==1 || j==2  && i>0   && i<(p->n[0])-1 && k>0   && k<(p->n[2])-1 )
	       grad=0;
   #else
	   if((j==(p->n[1])-2) || (j==(p->n[1])-3)  && i>0   && i<(p->n[0])-1  )
	       grad=0;
	   else if(j==1 || j==2  && i>0   && i<(p->n[0])-1  )
	       grad=0;
   #endif


  }
   #ifdef USE_SAC_3D
 else if((dir == 2)    /*&& k>0*/ && k<((p->n[2])-1))
 {
    grad=(  ( wmod[encode3_MODID(p,i,j,k+1,field)]-wmod[encode3_MODID(p,i,j,k,field)])/((p->dx[2]))    );

   if((k==(p->n[2])-2) || (k==(p->n[2])-3)  && i>0   && i<(p->n[0])-1  && j>0   && j<(p->n[1])-1 )
       grad=0;
   else if(k==1 || k==2  && i>0   && i<(p->n[0])-1  && j>0   && j<(p->n[1])-1 )
       grad=0;


  }
  #endif
 return grad;
}



__device__ __host__
real grad13_MODID(real *wmod,struct params *p,int *ii,int field,int dir)
{
  real grad=0;
   int i,j,k;
   i=ii[0];
   j=ii[1];
   k=0;
   #ifdef USE_SAC_3D
    k=ii[2];
   #endif


 if((dir == 0) && i>0 && i<((p->n[0])-1))
 {
    grad=(  ( wmod[encode3_MODID(p,i+1,j,k,field)]-wmod[encode3_MODID(p,i-1,j,k,field)]) /((p->dx[0]))/2.0    );

   #ifdef USE_SAC_3D
	   if((i==(p->n[0])-2) || (i==(p->n[0])-3)  && j>0   && j<(p->n[1])-1 && k>0   && k<(p->n[2])-1 )
	       grad=0;
	   else if(i==1 || i==2  && j>0   && j<(p->n[1])-1 && k>0   && k<(p->n[2])-1 )
	       grad=0;
   #else
	   if((i==(p->n[0])-2) || (i==(p->n[0])-3)  && j>0   && j<(p->n[1])-1  )
	       grad=0;
	   else if(i==1 || i==2  && j>0   && j<(p->n[1])-1  )
	       grad=0;
   #endif
 }
 else if((dir == 1)    && j>0 && j<((p->n[1])-1))
 {
    grad=(  ( wmod[encode3_MODID(p,i,j+1,k,field)]-wmod[encode3_MODID(p,i,j-1,k,field)])/((p->dx[1]))/2.0    );
   #ifdef USE_SAC_3D
	   if((j==(p->n[1])-2) || (j==(p->n[1])-3)  && i>0   && i<(p->n[0])-1 && k>0   && k<(p->n[2])-1 )
	       grad=0;
	   else if(j==1 || j==2  && i>0   && i<(p->n[0])-1 && k>0   && k<(p->n[2])-1 )
	       grad=0;
   #else
	   if((j==(p->n[1])-2) || (j==(p->n[1])-3)  && i>0   && i<(p->n[0])-1  )
	       grad=0;
	   else if(j==1 || j==2  && i>0   && i<(p->n[0])-1  )
	       grad=0;
   #endif


  }
   #ifdef USE_SAC_3D
 else if((dir == 2)    && k>0 && k<((p->n[2])-1))
 {
    grad=(  ( wmod[encode3_MODID(p,i,j,k+1,field)]-wmod[encode3_MODID(p,i,j,k-1,field)])/((p->dx[2]))/2.0    );

   if((k==(p->n[2])-2) || (k==(p->n[2])-3)  && i>0   && i<(p->n[0])-1  && j>0   && j<(p->n[1])-1 )
       grad=0;
   else if(k==1 || k==2  && i>0   && i<(p->n[0])-1  && j>1   && j<(p->n[1])-1 )
       grad=0;


  }
  #endif
 return grad;
}




/*****************************************************/




__device__ __host__
real grad3dn_MODID(real *wmod, real *wd,struct params *p,int *ii,int field,int dir)
{


 real grad=0;

 

 switch(dir)
 {
   case 0:
 
#ifdef USE_SAC_3D
  #ifdef USE_DORDER3
 if(ii[0]>2 && ii[0]<((p->n[0])-3) )
  grad=(  ( ((3*wmod[encode3_MODID(p,ii[0]+1,ii[1],ii[2],field)]-3*wmod[encode3_MODID(p,ii[0]-1,ii[1],ii[2],field)]+3.0*(wmod[encode3_MODID(p,ii[0]-2,ii[1],ii[2],field)]-wmod[encode3_MODID(p,ii[0]+2,ii[1],ii[2],field)])/5.0-(wmod[encode3_MODID(p,ii[0]-3,ii[1],ii[2],field)]-wmod[encode3_MODID(p,ii[0]+3,ii[1],ii[2],field)])/15.0)/2.0))/(2.0*(wd[fencode3_MODID(p,ii,delx1)]))    );

  #else
if(ii[0]>1 && ii[0]<((p->n[0])-2) )
 grad=(  ( ((8*wmod[encode3_MODID(p,ii[0]+1,ii[1],ii[2],field)]-8*wmod[encode3_MODID(p,ii[0]-1,ii[1],ii[2],field)]+wmod[encode3_MODID(p,ii[0]-2,ii[1],ii[2],field)]-wmod[encode3_MODID(p,ii[0]+2,ii[1],ii[2],field)])/6.0))/(2.0*(wd[fencode3_MODID(p,ii,delx1)]))    );
 #endif

#ifdef USE_MPI
if(p->boundtype[field][dir][0] !=1  )
  if(p->mpiupperb[dir]==1  )
#else
if(p->boundtype[field][dir][0] !=0  )
#endif
{

  ;//for OZT test using MPI use this directive further clarification needed
  #ifndef USE_MPI
   if((ii[0]==(p->n[0])-3) || (ii[0]==(p->n[0])-4)  && ii[1]>1   && ii[1]<(p->n[1])-2 && ii[2]>1   && ii[2]<(p->n[2])-2  )
       grad=0;
   else if(ii[0]==2 || ii[0]==3  && ii[1]>1   && ii[1]<(p->n[1])-2 && ii[2]>1   && ii[2]<(p->n[2])-2  )
       grad=0;
  #endif
}
#else  

  #ifdef USE_DORDER3
if(ii[0]>2 && ii[0]<((p->n[0])-3) )
 grad=(  ( ((3*wmod[encode3_MODID(p,ii[0]+1,ii[1],0,field)]-3*wmod[encode3_MODID(p,ii[0]-1,ii[1],0,field)]+3.0*(wmod[encode3_MODID(p,ii[0]-2,ii[1],0,field)]-wmod[encode3_MODID(p,ii[0]+2,ii[1],0,field)])/5.0-(wmod[encode3_MODID(p,ii[0]-3,ii[1],0,field)]-wmod[encode3_MODID(p,ii[0]+3,ii[1],0,field)])/15.0)/2.0))/(2.0*(wd[fencode3_MODID(p,ii,delx1)]))    );

  #else
if(ii[0]>1 && ii[0]<((p->n[0])-2) )
 grad=(  ( ((8*wmod[encode3_MODID(p,ii[0]+1,ii[1],0,field)]-8*wmod[encode3_MODID(p,ii[0]-1,ii[1],0,field)]+wmod[encode3_MODID(p,ii[0]-2,ii[1],0,field)]-wmod[encode3_MODID(p,ii[0]+2,ii[1],0,field)])/6.0))/(2.0*(wd[fencode3_MODID(p,ii,delx1)]))    );
 #endif
#ifdef USE_MPI
if(p->boundtype[field][dir][0] !=1  )
  if(p->mpiupperb[dir]==1  )
#else
if(p->boundtype[field][dir][0] !=0  )
#endif
{

  ;//for OZT test using MPI use this directive further clarification needed
  #ifndef USE_MPI
   if((ii[0]==(p->n[0])-3) || (ii[0]==(p->n[0])-4)  && ii[1]>1   && ii[1]<(p->n[1])-2  )
       grad=0;
   else if(ii[0]==2 || ii[0]==3  && ii[1]>1   && ii[1]<(p->n[1])-2  )
       grad=0;
  #endif
}
#endif



   break;

   case 1:

#ifdef USE_SAC_3D

  #ifdef USE_DORDER3
 if(ii[1]>2 && ii[1]<((p->n[1])-3) )
  grad=(  ( ((3*wmod[encode3_MODID(p,ii[0],ii[1]+1,ii[2],field)]-3*wmod[encode3_MODID(p,ii[0],ii[1]-1,ii[2],field)]+3.0*(wmod[encode3_MODID(p,ii[0],ii[1]-2,ii[2],field)]-wmod[encode3_MODID(p,ii[0],ii[1]+2,ii[2],field)])/5.0-(wmod[encode3_MODID(p,ii[0],ii[1]-3,ii[2],field)]-wmod[encode3_MODID(p,ii[0],ii[1]+3,ii[2],field)])/15.0)/2.0))/(2.0*(wd[fencode3_MODID(p,ii,delx2)]))    );

#else
if( ii[1] >1 &&  ii[1]<((p->n[1])-2))
	grad=(  ( ((8*wmod[encode3_MODID(p,ii[0],ii[1]+1,ii[2],field)]-8*wmod[encode3_MODID(p,ii[0],ii[1]-1,ii[2],field)]+wmod[encode3_MODID(p,ii[0],ii[1]-2,ii[2],field)]-wmod[encode3_MODID(p,ii[0],ii[1]+2,ii[2],field)])/6.0))/(2.0*(wd[fencode3_MODID(p,ii,delx2)]))    );
 #endif
#ifdef USE_MPI
if(p->boundtype[field][dir][0] !=1  )
  if(p->mpiupperb[dir]==1  )
#else
if(p->boundtype[field][dir][0] !=0  )
#endif
{
  ;//for OZT test using MPI use this directive further clarification needed
  #ifndef USE_MPI
   if((ii[1]==(p->n[1])-3) || (ii[1]==(p->n[1])-4)  && ii[0]>1   && ii[0]<(p->n[0])-2  && ii[2]>1   && ii[2]<(p->n[2])-2  )
       grad=0;
   else if(ii[1]==2 || ii[1]==3  && ii[0]>1   && ii[0]<(p->n[0])-2  && ii[2]>1   && ii[2]<(p->n[2])-2  )
       grad=0;
  #endif
}
#else

  #ifdef USE_DORDER3
if(ii[1]>2 && ii[1]<((p->n[1])-3) )
 grad=(  ( ((3*wmod[encode3_MODID(p,ii[0],ii[1]+1,0,field)]-3*wmod[encode3_MODID(p,ii[0],ii[1]-1,0,field)]+3.0*(wmod[encode3_MODID(p,ii[0],ii[1]-2,0,field)]-wmod[encode3_MODID(p,ii[0],ii[1]+2,0,field)])/5.0-(wmod[encode3_MODID(p,ii[0],ii[1]-3,0,field)]-wmod[encode3_MODID(p,ii[0],ii[1]+3,0,field)])/15.0)/2.0))/(2.0*(wd[fencode3_MODID(p,ii,delx2)]))    );

#endif
if( ii[1] >1 &&  ii[1]<((p->n[1])-2))
	grad=(  ( ((8*wmod[encode3_MODID(p,ii[0],ii[1]+1,0,field)]-8*wmod[encode3_MODID(p,ii[0],ii[1]-1,0,field)]+wmod[encode3_MODID(p,ii[0],ii[1]-2,0,field)]-wmod[encode3_MODID(p,ii[0],ii[1]+2,0,field)])/6.0))/(2.0*(wd[fencode3_MODID(p,ii,delx2)]))    );

#ifdef USE_MPI
if(p->boundtype[field][dir][0] !=1  )
  if(p->mpiupperb[dir]==1  )
#else
if(p->boundtype[field][dir][0] !=0  )
#endif
{

  ;//for OZT test using MPI use this directive further clarification needed
  #ifndef USE_MPI
   if((ii[1]==(p->n[1])-3) || (ii[1]==(p->n[1])-4)  && ii[0]>1   && ii[0]<(p->n[0])-2  )
       grad=0;
   else if(ii[1]==2 || ii[1]==3  && ii[0]>1   && ii[0]<(p->n[0])-2  )
       grad=0;
  #endif
}
#endif
   break;


   case 2:

#ifdef USE_SAC_3D
  #ifdef USE_DORDER3
 if(ii[2]>2 && ii[2]<((p->n[2])-3) )
  grad=(  ( ((3*wmod[encode3_MODID(p,ii[0],ii[1],ii[2]+1,field)]-3*wmod[encode3_MODID(p,ii[0],ii[1],ii[2]-1,field)]+3.0*(wmod[encode3_MODID(p,ii[0],ii[1],ii[2]-2,field)]-wmod[encode3_MODID(p,ii[0],ii[1],ii[2]+2,field)])/5.0-(wmod[encode3_MODID(p,ii[0],ii[1],ii[2]-3,field)]-wmod[encode3_MODID(p,ii[0],ii[1],ii[2]+3,field)])/15.0)/2.0))/(2.0*(wd[fencode3_MODID(p,ii,delx3)]))    );

#else
if( ii[2] >1 &&  ii[2]<((p->n[2])-2))
	grad=(  ( ((8*wmod[encode3_MODID(p,ii[0],ii[1],ii[2]+1,field)]-8*wmod[encode3_MODID(p,ii[0],ii[1],ii[2]-1,field)]+wmod[encode3_MODID(p,ii[0],ii[1],ii[2]-2,field)]-wmod[encode3_MODID(p,ii[0],ii[1],ii[2]+2,field)])/6.0))/(2.0*(wd[fencode3_MODID(p,ii,delx3)]))    );
#endif

#ifdef USE_MPI
if(p->boundtype[field][dir][0] !=1  )
  if(p->mpiupperb[dir]==1  )
#else
if(p->boundtype[field][dir][0] !=0  )
#endif
{

  ;//for OZT test using MPI use this directive further clarification needed
  #ifndef USE_MPI
   if((ii[2]==(p->n[2])-3) || (ii[2]==(p->n[2])-4)  && ii[0]>1   && ii[0]<(p->n[0])-2 && ii[1]>1   && ii[1]<(p->n[1])-2  )
       grad=0;
   else if(ii[2]==2 || ii[2]==3  && ii[0]>1   && ii[0]<(p->n[0])-2 && ii[1]>1   && ii[1]<(p->n[1])-2  )
       grad=0;
  #endif
}
#endif
   break;

}



 return grad;


}




__device__ __host__
real grad1l3n_MODID(real *wmod, real *wd,struct params *p,int *ii,int field,int dir)
{
 real grad=0;
   int i,j,k;
   i=ii[0];
   j=ii[1];
   k=0;
   #ifdef USE_SAC_3D
    k=ii[2];
   #endif


 if((dir == 0) && i>0 && i<((p->n[0])))
 {
    grad=(  ( wmod[encode3_MODID(p,i,j,k,field)]-wmod[encode3_MODID(p,i-1,j,k,field)]) /((wd[fencode3_MODID(p,ii,delx1)]))    );

   if(p->boundtype[field][dir][0] !=0)
	{
	   #ifdef USE_SAC_3D
		   if((i==(p->n[0])-2) || (i==(p->n[0])-3)  && j>0   && j<(p->n[1])-1 && k>1   && k<(p->n[2])-1 )
		       grad=0;
		   else if(i==1 || i==2  && j>0   && j<(p->n[1])-1 && k>0   && k<(p->n[2])-1 )
		       grad=0;
	   #else
		   if((i==(p->n[0])-2) || (i==(p->n[0])-3)  && j>0   && j<(p->n[1])-1  )
		       grad=0;
		   else if(i==1 || i==2  && j>0   && j<(p->n[1])-1  )
		       grad=0;
	   #endif
	}
 }
 else if((dir == 1)    && j>0 && j<((p->n[1])))
 {
    grad=(  ( wmod[encode3_MODID(p,i,j,k,field)]-wmod[encode3_MODID(p,i,j-1,k,field)])/((wd[fencode3_MODID(p,ii,delx2)]))    );

  if(p->boundtype[field][dir][0] !=0)
  {
   #ifdef USE_SAC_3D
	   if((j==(p->n[1])-2) || (j==(p->n[1])-3)  && i>0   && i<(p->n[0])-1 && k>0   && k<(p->n[2])-1 )
	       grad=0;
	   else if(j==1 || j==2  && i>0   && i<(p->n[0])-1 && k>0   && k<(p->n[2])-1 )
	       grad=0;
   #else
	   if((j==(p->n[1])-2) || (j==(p->n[1])-3)  && i>0   && i<(p->n[0])-1  )
	       grad=0;
	   else if(j==1 || j==2  && i>0   && i<(p->n[0])-1  )
	       grad=0;
   #endif
   }

  }
   #ifdef USE_SAC_3D
 else if((dir == 2)    && k>0 && k<((p->n[2])))
 {
    grad=(  ( wmod[encode3_MODID(p,i,j,k,field)]-wmod[encode3_MODID(p,i,j,k-1,field)])/((wd[fencode3_MODID(p,ii,delx3)]))    );

 if(p->boundtype[field][dir][0] !=0)
 {
   if((k==(p->n[2])-2) || (k==(p->n[2])-3)  && i>0   && i<(p->n[0])-1  && j>0   && j<(p->n[1])-1 )
       grad=0;
   else if(k==1 || k==2  && i>0   && i<(p->n[0])-1  && j>0   && j<(p->n[1])-1 )
       grad=0;
 }


  }
  #endif
 return grad;

}

__device__ __host__
real grad1r3n_MODID(real *wmod, real *wd,struct params *p,int *ii,int field,int dir)
{
  real grad=0;
   int i,j,k;
   i=ii[0];
   j=ii[1];
   k=0;
   #ifdef USE_SAC_3D
    k=ii[2];
   #endif


 if((dir == 0) && /*i>0 &&*/ i<((p->n[0])-1))
 {
    grad=(  ( wmod[encode3_MODID(p,i+1,j,k,field)]-wmod[encode3_MODID(p,i,j,k,field)]) /((wd[fencode3_MODID(p,ii,delx1)]))    );


   if(p->boundtype[field][dir][0] !=0)
   {
   #ifdef USE_SAC_3D
	   if((i==(p->n[0])-2) || (i==(p->n[0])-3)  && j>0   && j<(p->n[1])-1 && k>0   && k<(p->n[2])-1 )
	       grad=0;
	   else if(i==1 || i==2  && j>0   && j<(p->n[1])-1 && k>0   && k<(p->n[2])-1 )
	       grad=0;
   #else
	   if((i==(p->n[0])-2) || (i==(p->n[0])-3)  && j>0   && j<(p->n[1])-1  )
	       grad=0;
	   else if(i==1 || i==2  && j>0   && j<(p->n[1])-1  )
	       grad=0;
   #endif
   }
 }
 else if((dir == 1)    /*&& j>0*/ && j<((p->n[1])-1))
 {
    grad=(  ( wmod[encode3_MODID(p,i,j+1,k,field)]-wmod[encode3_MODID(p,i,j,k,field)])/((wd[fencode3_MODID(p,ii,delx2)]))    );


  if(p->boundtype[field][dir][0] !=0)
  {
   #ifdef USE_SAC_3D
	   if((j==(p->n[1])-2) || (j==(p->n[1])-3)  && i>0   && i<(p->n[0])-1 && k>0   && k<(p->n[2])-1 )
	       grad=0;
	   else if(j==1 || j==2  && i>0   && i<(p->n[0])-1 && k>0   && k<(p->n[2])-1 )
	       grad=0;
   #else
	   if((j==(p->n[1])-2) || (j==(p->n[1])-3)  && i>0   && i<(p->n[0])-1  )
	       grad=0;
	   else if(j==1 || j==2  && i>0   && i<(p->n[0])-1  )
	       grad=0;
   #endif
   }


  }
   #ifdef USE_SAC_3D
 else if((dir == 2)    /*&& k>0*/ && k<((p->n[2])-1))
 {
    grad=(  ( wmod[encode3_MODID(p,i,j,k+1,field)]-wmod[encode3_MODID(p,i,j,k,field)])/((wd[fencode3_MODID(p,ii,delx3)]))    );

if(p->boundtype[field][dir][0] !=0)
 {
   if((k==(p->n[2])-2) || (k==(p->n[2])-3)  && i>0   && i<(p->n[0])-1  && j>0   && j<(p->n[1])-1 )
       grad=0;
   else if(k==1 || k==2  && i>0   && i<(p->n[0])-1  && j>0   && j<(p->n[1])-1 )
       grad=0;
  }

  }
  #endif
 return grad;
}



__device__ __host__
real grad13n_MODID(real *wmod, real *wd,struct params *p,int *ii,int field,int dir)
{
  real grad=0;
   int i,j,k;
   i=ii[0];
   j=ii[1];
   k=0;
   #ifdef USE_SAC_3D
    k=ii[2];
   #endif


 if((dir == 0) && i>0 && i<((p->n[0])-1))
 {
    grad=(  ( wmod[encode3_MODID(p,i+1,j,k,field)]-wmod[encode3_MODID(p,i-1,j,k,field)]) /((wd[fencode3_MODID(p,ii,delx1)]))/2.0    );


	if(p->boundtype[field][dir][0] !=0)
	{
	   #ifdef USE_SAC_3D
		   if((i==(p->n[0])-2) || (i==(p->n[0])-3)  && j>0   && j<(p->n[1])-1 && k>0   && k<(p->n[2])-1 )
		       grad=0;
		   else if(i==1 || i==2  && j>0   && j<(p->n[1])-1 && k>0   && k<(p->n[2])-1 )
		       grad=0;
	   #else
		   if((i==(p->n[0])-2) || (i==(p->n[0])-3)  && j>0   && j<(p->n[1])-1  )
		       grad=0;
		   else if(i==1 || i==2  && j>0   && j<(p->n[1])-1  )
		       grad=0;
	   #endif
	}
 }
 else if((dir == 1)    && j>0 && j<((p->n[1])-1))
 {
    grad=(  ( wmod[encode3_MODID(p,i,j+1,k,field)]-wmod[encode3_MODID(p,i,j-1,k,field)])/((wd[fencode3_MODID(p,ii,delx2)]))/2.0    );

	if(p->boundtype[field][dir][0] !=0)
	{
	   #ifdef USE_SAC_3D
		   if((j==(p->n[1])-2) || (j==(p->n[1])-3)  && i>0   && i<(p->n[0])-1 && k>0   && k<(p->n[2])-1 )
		       grad=0;
		   else if(j==1 || j==2  && i>0   && i<(p->n[0])-1 && k>0   && k<(p->n[2])-1 )
		       grad=0;
	   #else
		   if((j==(p->n[1])-2) || (j==(p->n[1])-3)  && i>0   && i<(p->n[0])-1  )
		       grad=0;
		   else if(j==1 || j==2  && i>0   && i<(p->n[0])-1  )
		       grad=0;
	   #endif
	}

  }
   #ifdef USE_SAC_3D
 else if((dir == 2)    && k>0 && k<((p->n[2])-1))
 {
    grad=(  ( wmod[encode3_MODID(p,i,j,k+1,field)]-wmod[encode3_MODID(p,i,j,k-1,field)])/((wd[fencode3_MODID(p,ii,delx3)]))/2.0    );

	if(p->boundtype[field][dir][0] !=0)
	{
	   if((k==(p->n[2])-2) || (k==(p->n[2])-3)  && i>0   && i<(p->n[0])-1  && j>0   && j<(p->n[1])-1 )
	       grad=0;
	   else if(k==1 || k==2  && i>0   && i<(p->n[0])-1  && j>1   && j<(p->n[1])-1 )
	       grad=0;
	}

  }
  #endif
 return grad;
}








__device__ __host__
void bc3_cont_MODID(real *wt, struct params *p,int *ii, int f) {

   
int i,j,k;
i=ii[0];
j=ii[1];
k=0;
        #ifdef USE_SAC_3D
          k=ii[2];

                /*if(i<2 && j<2  && k<2)
                {
                 if(i==j==k )
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,2,2,2,f)];
                }
                else if(i<2 && j>((p->n[1])-3)  &&  k<2)
                {
                  if(i==(j-(p->n[1]))  && k==(j-(p->n[1])))                  
                     wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,2,(p->n[1])-3,2,f)];                     
                }
                else if(i>((p->n[0])-3) && j<2 && k<2)
                {
                  if(j==(i-(p->n[0]))  && k==(i-(p->n[0])))                  
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,((p->n[0])-3),2,2,f)];                  
                }
                else if(i>((p->n[0])-3) && j>((p->n[1])-3) && k<2)
                {
                  if(i==j  && k==(i-(p->n[0])) )                  
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,((p->n[0])-3),((p->n[1])-3),2,f)];                                                  
                }
                else if(i<2 && j<2  && k>((p->n[2])-3))
                {
                 if(i==j && k==(i-(p->n[0]))  )
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,2,2,(p->n[2])-3,f)];
                } 
                else if(i>((p->n[0])-3) && j<2  && k>((p->n[2])-3))
                {
                 if(i==k && j==(i-(p->n[0]))  )
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-3,2,(p->n[2])-3,f)];
                }
                else if(i<2 && j>((p->n[1])-3)  && k>((p->n[2])-3) )
                {
                 if(j==k && i==(j-(p->n[1]))  )
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,2,(p->n[1])-3,(p->n[2])-3,f)];
                } 
                else if(i>((p->n[0])-3) && j>((p->n[1])-3)  && k>((p->n[2])-3) )
                {
                 if(i==j==k  )
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-3,(p->n[1])-3,(p->n[2])-3,f)];
                }                     
                else*/ if(i==0 || i==1  && ((p->boundtype[f][0][0])==3))                
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,2,j,k,f)];              
                else if((i==((p->n[0])-1)) || (i==((p->n[0])-2)) && ((p->boundtype[f][0][0])==3))                
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-3,j,k,f)];                            
                else if(j==0 || j==1 && ((p->boundtype[f][1][0])==3))                
                   wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,2,k,f)];                    
                else if((j==((p->n[1])-1)) || (j==((p->n[1])-2))  && ((p->boundtype[f][1][0])==3))                
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,(p->n[1])-3,k,f)];
                else if(k==0 || k==1  && ((p->boundtype[f][2][0])==3))                
                   wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,2,k,f)];                    
                else if((k==((p->n[2])-1)) || (k==((p->n[2])-2))  && ((p->boundtype[f][2][0])==3))                
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,(p->n[2])-3,f)];

        #else
             /*if(i<2 && j<2)
                {
                  if(i==j)
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,2,j,f)];
                  else                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,2,f)];                  
                }
                else if(i<2 && j>((p->n[1])-3))
                {
                  if(i==(j-(p->n[1])))                  
                     wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,2,j,f)];                     
                  else                  
                     wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,((p->n[1])-3),f)];                   
                }
                else if(i>((p->n[0])-3) && j<2)
                {
                  if((i-(p->n[0]))==j)                  
                    //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(i-3),j,f)];
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,((p->n[0])-3),j,f)];                  
                  else                  
                   // wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,j+2,f)];
                   wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,2,f)];                        
                }
                else if(i>((p->n[0])-3) && j>((p->n[1])-3))
                {
                  if(i==j)                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(i-3),j,f)];                   
                  else                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(j-3),f)];                  
                }                       
                else*/ if(i==0 || i==1 && ((p->boundtype[f][0][0])==3))                
   
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,2,j,f)];              
                else if((i==((p->n[0])-1)) || (i==((p->n[0])-2))  && ((p->boundtype[f][0][0])==3))                
                  //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(i-3),j,f)];    
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-3,j,f)];                            
                else if(j==0 || j==1  && ((p->boundtype[f][1][0])==3))                
                  //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,j+2,f)]; 
                   wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,2,f)];                    
                else if((j==((p->n[1])-1)) || (j==((p->n[1])-2))  && ((p->boundtype[f][1][0])==3))                
                  //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(j-3),f)];
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-3,f)];
                
         #endif



}

__device__ __host__
void bc3_cont_dir_MODID(real *wt, struct params *p,int *ii, int f, int dir) {

   
int i,j,k;
i=ii[0];
j=ii[1];
k=0;
        #ifdef USE_SAC_3D
          k=ii[2];

               /* if(i<2 && j<2  && k<2)
                {
                 if(i==j==k )
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,2,2,2,f)];
                }
                else if(i<2 && j>((p->n[1])-3)  &&  k<2)
                {
                  if(i==(j-(p->n[1]))  && k==(j-(p->n[1])))                  
                     wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,2,(p->n[1])-3,2,f)];                     
                }
                else if(i>((p->n[0])-3) && j<2 && k<2)
                {
                  if(j==(i-(p->n[0]))  && k==(i-(p->n[0])))                  
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,((p->n[0])-3),2,2,f)];                  
                }
                else if(i>((p->n[0])-3) && j>((p->n[1])-3) && k<2)
                {
                  if(i==j  && k==(i-(p->n[0])) )                  
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,((p->n[0])-3),((p->n[1])-3),2,f)];                                                  
                }
                else if(i<2 && j<2  && k>((p->n[2])-3) )
                {
                 if(i==j && k==(i-(p->n[0]))  )
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,2,2,(p->n[2])-3,f)];
                } 
                else if(i>((p->n[0])-3) && j<2  && k>((p->n[2])-3)  )
                {
                 if(i==k && j==(i-(p->n[0]))  )
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-3,2,(p->n[2])-3,f)];
                }
                else if(i<2 && j>((p->n[1])-3)  && k>((p->n[2])-3) )
                {
                 if(j==k && i==(j-(p->n[1]))  )
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,2,(p->n[1])-3,(p->n[2])-3,f)];
                } 
                else if(i>((p->n[0])-3) && j>((p->n[1])-3)  && k>((p->n[2])-3) )
                {
                 if(i==j==k  )
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-3,(p->n[1])-3,(p->n[2])-3,f)];
                }                     
                else*/ if((i==0 || i==1)  && dir==0 && ((p->boundtype[f][dir][0])==3))                
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,2,j,k,f)];              
                else if(((i==((p->n[0])-1)) || (i==((p->n[0])-2)))  && dir==0 && ((p->boundtype[f][dir][0])==3))                
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-3,j,k,f)];                            
                else if((j==0 || j==1)  && dir==1 && ((p->boundtype[f][dir][0])==3))                
                   wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,2,k,f)];                    
                else if(((j==((p->n[1])-1)) || (j==((p->n[1])-2)))  && dir==1  && ((p->boundtype[f][dir][0])==3))                
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,(p->n[1])-3,k,f)];
                else if(k==0 || k==1  && dir==2  && ((p->boundtype[f][dir][0])==3))                
                   wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,2,k,f)];                    
                else if((k==((p->n[2])-1)) || (k==((p->n[2])-2))  && dir==2  && ((p->boundtype[f][dir][0])==3))                
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,(p->n[2])-3,f)];

        #else
             /*if(i<2 && j<2)
                {
                  if(i==j)
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,2,j,f)];
                  else                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,2,f)];                  
                }
                else if(i<2 && j>((p->n[1])-3))
                {
                  if(i==(j-(p->n[1])))                  
                     wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,2,j,f)];                     
                  else                  
                     wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,((p->n[1])-3),f)];                   
                }
                else if(i>((p->n[0])-3) && j<2)
                {
                  if((i-(p->n[0]))==j)                  
                    //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(i-3),j,f)];
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,((p->n[0])-3),j,f)];                  
                  else                  
                   // wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,j+2,f)];
                   wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,2,f)];                        
                }
                else if(i>((p->n[0])-3) && j>((p->n[1])-3))
                {
                  if(i==j)                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(i-3),j,f)];                   
                  else                  
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(j-3),f)];                  
                }                       
                else */

                if((i==0 || i==1)  && dir==0 && ((p->boundtype[f][dir][0])==3))                
                  //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i+2,j,f)];   
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,2,j,f)];              
                else if(((i==((p->n[0])-1)) || (i==((p->n[0])-2)))  && dir==0  && ((p->boundtype[f][dir][0])==3))                
                  //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(i-2),j,f)];    
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-3,j,f)];                            
                else if((j==0 || j==1)  && dir==1  && ((p->boundtype[f][dir][0])==3))                
                  //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,j+2,f)]; 
                   wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,2,f)];                    
                else if(((j==((p->n[1])-1)) || (j==((p->n[1])-2)))  && dir==1  && ((p->boundtype[f][dir][0])==3))                
                  //wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(j-2),f)];
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-3,f)];
                
         #endif



}


__device__ __host__
void bc3_cont_cd4_MODID(real *wt, struct params *p,int *ii, int f) {


int i,j,k;
i=ii[0];
j=ii[1];
k=0;
        #ifdef USE_SAC_3D
          k=ii[2];
            if((p->boundtype[f][0][0])==4)
                if(i==0)              
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,4,j,k,f)];
                else if(i==1)                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,3,j,k,f)];
                else if( i==((p->n[0])-1))               
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-5,j,k,f)];
                else if (i==((p->n[0])-2))                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-4,j,k,f)];
               

            if((p->boundtype[f][1][0])==4)
                if(j==0)               
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,4,k,f)];
                else if(j==1)                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,3,k,f)];
                else if (j== ((p->n[1])-1))               
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,(p->n[1])-5,k,f)];
               else if (j== ((p->n[1])-2))                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,(p->n[1])-4,k,f)];



            if((p->boundtype[f][2][0])==4)
                if(k==0)               
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,4,f)];
                else if(k==1)                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,3,f)];
                else if (k== ((p->n[2])-1))               
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,(p->n[2])-5,f)];
               else if (k== ((p->n[2])-2))                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,(p->n[2])-4,f)];
        #else
        if((p->boundtype[f][0][0])==4)   
                if(i==0)              
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,3,j,f)];
                else if(i==1)                
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,2,j,f)];
                else if( i==((p->n[0])-1))               
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-4,j,f)];
                else if (i==((p->n[0])-2))                
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-3,j,f)];
               

            if((p->boundtype[f][1][0])==4)
                if(j==0)               
                  // wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,4,f)];
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,3,f)];
                else if(j==1)                
                  //  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,3,f)];
                   wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,2,f)];
                else if (j== ((p->n[1])-1))               
                  //  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-5,f)];
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-4,f)];
               else if (j== ((p->n[1])-2))                
                  //  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-4,f)];
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-3,f)];
         #endif

}



__device__ __host__
void bc3_cont_cd4_dir_MODID(real *wt, struct params *p,int *ii, int f, int dir) {


int i,j,k;
i=ii[0];
j=ii[1];
k=0;
        #ifdef USE_SAC_3D

          k=ii[2];
                      if((p->boundtype[f][dir][0])==4)
                      {
                if((i==0 || i==1) && dir==0)              
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,4-i,j,k,f)];             
                else if((( i==((p->n[0])-1)   ))  && dir==0)               
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i-4,j,k,f)];
                else if(((  i==((p->n[0])-2) ))  && dir==0)               
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i-2,j,k,f)];
              

                if((j==0 || j==1) && dir==1)              
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,4-j,k,f)];             
                else if((( j==((p->n[1])-1)   ))  && dir==1)               
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j-4,k,f)];
                else if(((  j==((p->n[1])-2) ))  && dir==1)               
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j-2,k,f)];


                if((k==0 || k==1) && dir==2)              
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,4-k,f)];             
                else if((( k==((p->n[2])-1)   ))  && dir==2)               
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,k-4,f)];
                else if(((  k==((p->n[2])-2) ))  && dir==2)               
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,k-2,f)];

                    }
               
        #else

                          if((p->boundtype[f][dir][0])==4)
                          {
                if((i==0 || i==1) && dir==0)              
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,4-i,j,k,f)];             
                else if((( i==((p->n[0])-1)   ))  && dir==0)               
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i-4,j,k,f)];
                else if(((  i==((p->n[0])-2) ))  && dir==0)               
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i-2,j,k,f)];
              

                if((j==0 || j==1) && dir==1)              
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,4-j,k,f)];             
                else if((( j==((p->n[1])-1)   ))  && dir==1)               
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j-4,k,f)];
                else if(((  j==((p->n[1])-2) ))  && dir==1)               
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j-2,k,f)];
                    }
         #endif

}




__device__ __host__
void bc3_setfixed_dir_MODID(real *wt, struct params *p,struct bparams *bp,int *ii, int f,int dir) {


int i,j,k;
i=ii[0];
j=ii[1];
k=0;


        #ifdef USE_SAC_3D
          k=ii[2];
        #endif

          if((p->boundtype[f][dir][0])==5)   
                if(i==0 || i==1  && dir==0)                
                  bp->fixed1[encodefixed13_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,k,f)];                
                else if((i==((p->n[0])-1)) || (i==((p->n[0])-2)) && dir==0 )                
                  bp->fixed1[encodefixed13_MODID(p,1+(p->n[0])-i,j,k,f)]=wt[encode3_MODID(p,i,j,k,f)];                
                else if(j==0 || j==1  && dir==1 )                
                  bp->fixed2[encodefixed23_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,k,f)];                
                else if((j==((p->n[1])-1)) || (j==((p->n[1])-2))  && dir==1)                
                  bp->fixed2[encodefixed23_MODID(p,i,1+(p->n[1])-j,k,f)]=wt[encode3_MODID(p,i,j,k,f)];
           #ifdef USE_SAC_3D
                else if(k==0 || k==1 && dir==2)                
                  bp->fixed3[encodefixed33_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,k,f)];                
                else if((k==((p->n[2])-1)) || (k==((p->n[2])-2))  && dir==2)                
                  bp->fixed3[encodefixed33_MODID(p,i,j,1+(p->n[2])-k,f)]=wt[encode3_MODID(p,i,j,k,f)];
           #endif

}


__device__ __host__
void bc3_fixed_dir_MODID(real *wt, struct params *p, struct bparams *bp,int *ii, int f,int dir) {


int i,j,k;
i=ii[0];
j=ii[1];
k=0;


        #ifdef USE_SAC_3D
          k=ii[2];
        #endif

     if((p->boundtype[f][dir][0])==5)         
                if(i==0 || i==1  && dir==0)                
                  wt[encode3_MODID(p,i,j,k,f)]=bp->fixed1[encodefixed13_MODID(p,i,j,k,f)];                
                else if((i==((p->n[0])-1)) || (i==((p->n[0])-2)) && dir==0 )                
                  wt[encode3_MODID(p,i,j,k,f)]=bp->fixed1[encodefixed13_MODID(p,1+(p->n[0])-i,j,k,f)];                
                else if(j==0 || j==1  && dir==1 )                
                  wt[encode3_MODID(p,i,j,k,f)]=bp->fixed2[encodefixed23_MODID(p,i,j,k,f)];                
                else if((j==((p->n[1])-1)) || (j==((p->n[1])-2))  && dir==1)                
                  wt[encode3_MODID(p,i,j,k,f)]=bp->fixed2[encodefixed23_MODID(p,i,1+(p->n[1])-j,k,f)];
           #ifdef USE_SAC_3D
                else if(k==0 || k==1 && dir==2)                
                  wt[encode3_MODID(p,i,j,k,f)]=bp->fixed3[encodefixed33_MODID(p,i,j,k,f)];                
                else if((k==((p->n[2])-1)) || (k==((p->n[2])-2))  && dir==2)                
                  wt[encode3_MODID(p,i,j,k,f)]=bp->fixed3[encodefixed33_MODID(p,i,j,1+(p->n[2])-k,f)];
           #endif
               




}


__device__ __host__
void bc3_periodic1_dir_MODID(real *wt, struct params *p,int *ii, int f,int dir) {

int i,j,k;
i=ii[0];
j=ii[1];
k=0;
        #ifdef USE_SAC_3D
          k=ii[2];
          if((p->boundtype[f][dir][0])==0)   
                if((i==0 || i==1) && dir==0)                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-4+i,j,k,f)];
             
                else if(((i==((p->n[0])-1)) || (i==((p->n[0])-2))) && dir==0)                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,4-(p->n[0])+i,j,k,f)];

          if((p->boundtype[f][dir][0])==0)   
                if((j==0 || j==1) && dir==1)                
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,(p->n[1])-4+j,k,f)];

                else if(((j==((p->n[1])-1)) || (j==((p->n[1])-2))) && dir==1)                 
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,4-(p->n[1])+j,k,f)];

           if((p->boundtype[f][dir][0])==0)   
                if((k==0 || k==1) && dir==2)                
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,(p->n[2])-4+j,f)];

                else if(((k==((p->n[2])-1)) || (k==((p->n[2])-2))) && dir==2)                 
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,4-(p->n[2])+k,f)];
       #else


  	wt[encode3_MODID(p,i,j,k,f)]=(((i==0 || i==1  || i==((p->n[0])-1) || i==((p->n[0])-2)) && dir==0)?((i==0 || i==1) && dir==0)*wt[encode3_MODID(p,(p->n[0])-4+i,j,k,f)]+(((i==((p->n[0])-1)) || (i==((p->n[0])-2))) && dir==0)*wt[encode3_MODID(p,4-(p->n[0])+i,j,k,f)]:wt[encode3_MODID(p,i,j,k,f)]);
          // if((p->boundtype[f][dir][0])==0)   
            /*    if((i==0 || i==1) && dir==0   )                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-4+i,j,k,f)];

                else if(((i==((p->n[0])-1)) || (i==((p->n[0])-2))) && dir==0   )                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,4-(p->n[0])+i,j,k,f)];*/

          // if((p->boundtype[f][dir][0])==0)   
           /*     if((j==0 || j==1) && dir==1  )                
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,(p->n[1])-4+j,k,f)];

                else if(((j==((p->n[1])-1)) || (j==((p->n[1])-2))) && dir==1  )                 
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,4-(p->n[1])+j,k,f)];*/

       #endif


}


__device__ __host__
void bc3_symm_dir_MODID(real *wt, struct params *p,int *ii, int f,int dir) {

int i,j,k;
i=ii[0];
j=ii[1];
k=0;
        #ifdef USE_SAC_3D
          k=ii[2];
                    if((p->boundtype[f][dir][0])==6)   
                if(i==0  && dir==0)              
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,4,j,k,f)];
                else if(i==1 && dir==0)                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,3,j,k,f)];
                else if( i==((p->n[0])-1)  && dir==0)               
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-5,j,k,f)];
                else if (i==((p->n[0])-2)  && dir==0)                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-4,j,k,f)];
               

          if((p->boundtype[f][dir][0])==6)   
                if(j==0  && dir==1)               
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,4,k,f)];
                else if(j==1  && dir==1)                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,3,k,f)];
                else if (j== ((p->n[1])-1)  && dir==1)               
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,(p->n[1])-5,k,f)];
               else if (j== ((p->n[1])-2)  && dir==1)                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,(p->n[1])-4,k,f)];



          if((p->boundtype[f][dir][0])==6)   
                if(k==0 && dir==2)               
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,4,f)];
                else if(k==1 && dir==2)                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,3,f)];
                else if (k== ((p->n[2])-1) && dir==2)               
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,(p->n[2])-5,f)];
               else if (k== ((p->n[2])-2) && dir==2)                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,(p->n[2])-4,f)];
        #else
        if((p->boundtype[f][dir][0])==6)   
                if(i==0  && dir==0)              
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,3,j,f)];
                else if(i==1  && dir==0)                
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,2,j,f)];
                else if( i==((p->n[0])-1)  && dir==0)               
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-4,j,f)];
                else if (i==((p->n[0])-2)  && dir==0)                
                    wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,(p->n[0])-3,j,f)];
               

          if((p->boundtype[f][dir][0])==6)   
                if(j==0  && dir==1)               
                  // wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,4,f)];
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,3,f)];
                else if(j==1  && dir==1)                
                  //  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,3,f)];
                   wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,2,f)];
                else if (j== ((p->n[1])-1)  && dir==1)               
                  //  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-5,f)];
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-4,f)];
               else if (j== ((p->n[1])-2)  && dir==1)                
                  //  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-4,f)];
                  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-3,f)];
         #endif


}


__device__ __host__
void bc3_asymm_dir_MODID(real *wt, struct params *p,int *ii, int f,int dir) {

int i,j,k;
i=ii[0];
j=ii[1];
k=0;
        #ifdef USE_SAC_3D
          k=ii[2];
          if((p->boundtype[f][dir][0])==7)   
                if(i==0  && dir==0)              
                    wt[encode3_MODID(p,i,j,k,f)]=-wt[encode3_MODID(p,4,j,k,f)];
                else if(i==1 && dir==0)                
                    wt[encode3_MODID(p,i,j,k,f)]=-wt[encode3_MODID(p,3,j,k,f)];
                else if( i==((p->n[0])-1)  && dir==0)               
                    wt[encode3_MODID(p,i,j,k,f)]=-wt[encode3_MODID(p,(p->n[0])-5,j,k,f)];
                else if (i==((p->n[0])-2)  && dir==0)                
                    wt[encode3_MODID(p,i,j,k,f)]=-wt[encode3_MODID(p,(p->n[0])-4,j,k,f)];
               

          if((p->boundtype[f][dir][0])==7)   
                if(j==0  && dir==1)               
                    wt[encode3_MODID(p,i,j,k,f)]=-wt[encode3_MODID(p,i,4,k,f)];
                else if(j==1  && dir==1)                
                    wt[encode3_MODID(p,i,j,k,f)]=-wt[encode3_MODID(p,i,3,k,f)];
                else if (j== ((p->n[1])-1)  && dir==1)               
                    wt[encode3_MODID(p,i,j,k,f)]=-wt[encode3_MODID(p,i,(p->n[1])-5,k,f)];
               else if (j== ((p->n[1])-2)  && dir==1)                
                    wt[encode3_MODID(p,i,j,k,f)]=-wt[encode3_MODID(p,i,(p->n[1])-4,k,f)];



          if((p->boundtype[f][dir][0])==7)   
                if(k==0 && dir==2)               
                    wt[encode3_MODID(p,i,j,k,f)]=-wt[encode3_MODID(p,i,j,4,f)];
                else if(k==1 && dir==2)                
                    wt[encode3_MODID(p,i,j,k,f)]=-wt[encode3_MODID(p,i,j,3,f)];
                else if (k== ((p->n[2])-1) && dir==2)               
                    wt[encode3_MODID(p,i,j,k,f)]=-wt[encode3_MODID(p,i,j,(p->n[2])-5,f)];
               else if (k== ((p->n[2])-2) && dir==2)                
                    wt[encode3_MODID(p,i,j,k,f)]=-wt[encode3_MODID(p,i,j,(p->n[2])-4,f)];
        #else
           if((p->boundtype[f][dir][0])==7)     
                if(i==0  && dir==0)              
                    wt[fencode_MODID(p,i,j,f)]=-wt[fencode_MODID(p,3,j,f)];
                else if(i==1  && dir==0)                
                    wt[fencode_MODID(p,i,j,f)]=-wt[fencode_MODID(p,2,j,f)];
                else if( i==((p->n[0])-1)  && dir==0)               
                    wt[fencode_MODID(p,i,j,f)]=-wt[fencode_MODID(p,(p->n[0])-4,j,f)];
                else if (i==((p->n[0])-2)  && dir==0)                
                    wt[fencode_MODID(p,i,j,f)]=-wt[fencode_MODID(p,(p->n[0])-3,j,f)];
               

          if((p->boundtype[f][dir][0])==7)   
                if(j==0  && dir==1)               
                  // wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,4,f)];
                  wt[fencode_MODID(p,i,j,f)]=-wt[fencode_MODID(p,i,3,f)];
                else if(j==1  && dir==1)                
                  //  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,3,f)];
                   wt[fencode_MODID(p,i,j,f)]=-wt[fencode_MODID(p,i,2,f)];
                else if (j== ((p->n[1])-1)  && dir==1)               
                  //  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-5,f)];
                  wt[fencode_MODID(p,i,j,f)]=-wt[fencode_MODID(p,i,(p->n[1])-4,f)];
               else if (j== ((p->n[1])-2)  && dir==1)                
                  //  wt[fencode_MODID(p,i,j,f)]=wt[fencode_MODID(p,i,(p->n[1])-4,f)];
                  wt[fencode_MODID(p,i,j,f)]=-wt[fencode_MODID(p,i,(p->n[1])-3,f)];
         #endif


}



__device__ __host__
void bc3_periodic1_MODID(real *wt, struct params *p,int *ii, int f) {

int i,j,k;
i=ii[0];
j=ii[1];
k=0;
        #ifdef USE_SAC_3D
          k=ii[2];
          if((p->boundtype[f][0][0])==0)   
                if(i==0 || i==1 )                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-4+i,j,k,f)];
             
                else if((i==((p->n[0])-1)) || (i==((p->n[0])-2)))                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,4-(p->n[0])+i,j,k,f)];

          if((p->boundtype[f][1][0])==0)
                if(j==0 || j==1 )                
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,(p->n[1])-4+j,k,f)];

                else if((j==((p->n[1])-1)) || (j==((p->n[1])-2)) )                 
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,4-(p->n[1])+j,k,f)];
                  
          if((p->boundtype[f][2][0])==0)
                if(k==0 || k==1 )                
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,(p->n[2])-4+j,f)];

                else if((k==((p->n[2])-1)) || (k==((p->n[2])-2)) )                 
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,4-(p->n[2])+k,f)];
       #else
          if((p->boundtype[f][0][0])==0)
                if(i==0 || i==1 )                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-4+i,j,k,f)];

                else if((i==((p->n[0])-1)) || (i==((p->n[0])-2)))                
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,4-(p->n[0])+i,j,k,f)];

           if((p->boundtype[f][0][0])==0)
                if(j==0 || j==1 )                
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,(p->n[1])-4+j,k,f)];

                else if((j==((p->n[1])-1)) || (j==((p->n[1])-2)) )                 
                  wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,4-(p->n[1])+j,k,f)];

       #endif


}


__device__ __host__
void bc3_periodic2_MODID(real *wt, struct params *p,int *ii, int f) {

int i,j,k;
i=ii[0];
j=ii[1];
k=0;
        #ifdef USE_SAC_3D
          k=ii[2];
          
                if(i<2 && j<2  && k<2)
                {
                 if(i==j==k )
                   if((p->boundtype[f][1][0])==0)
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,(p->n[1])-4+j,k,f)];

                }
                else if(i<2 && j>((p->n[1])-3)  &&  k<2)
                {
                  if(i==(j-(p->n[1]))  && k==(j-(p->n[1])))
                  if((p->boundtype[f][1][0])==0)                  
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-4+i,j,k,f)];                                     
             
                }
                else if(i>((p->n[0])-3) && j<2 && k<2)
                {
                     if((p->boundtype[f][0][0])==0)
                  if(j==(i-(p->n[0]))  && k==(i-(p->n[0])))                  
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,4-(p->n[0])+i,j,k,f)];                                    
           
                }
                else if(i>((p->n[0])-3) && j>((p->n[1])-3) && k<2)
                {
                     if((p->boundtype[f][1][0])==0)
                  if(i==j  && k==(i-(p->n[0])))   
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,4-(p->n[1])+j,k,f)];                                    
                                 
                                            
                }
                else if(i<2 && j<2  && k>((p->n[2])-3))
                {
                     if((p->boundtype[f][2][0])==0)
                 if(i==j && i==(k-(p->n[2]))  )                 
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,j,(p->n[2])-4+k,f)];                                     

                } 
                else if(i>((p->n[0])-3) && j<2  && k>((p->n[2])-3))
                {
                     if((p->boundtype[f][0][0])==2)
                 if(i==k && j==(i-(p->n[0]))  )
                     wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,4-(p->n[0])+i,j,4-(p->n[2])+k,f)];   
                }
                else if(i<2 && j>((p->n[1])-3)  && k>((p->n[2])-3) )
                {
                     if((p->boundtype[f][2][0])==0)
                 if(j==k && i==(j-(p->n[1]))  )
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,4-(p->n[1])+j,4-(p->n[2])+k,f)]; 
                } 
                else if(i>((p->n[0])-3) && j>((p->n[1])-3)  && k>((p->n[2])-3) )
                {
                     if((p->boundtype[f][2][0])==0)
                 if(i==j==k  )
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,4-(p->n[0])+i,4-(p->n[1])+j,4-(p->n[2])+k,f)]; 
                }   

        #else

               if(i<2 && j<2)
                {
                      if((p->boundtype[f][0][0])==1)
                  if(i==j)
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,(p->n[1])-4+j,k,f)];
                  else              
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-4+i,j,k,f)];                                    
                }
                else if(i<2 && j>((p->n[1])-3))
                {
                     if((p->boundtype[f][0][0])==0)
                  if(i==(j-(p->n[1])))                  
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,(p->n[0])-4+i,j,k,f)];                                     
                  else                  
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,4-(p->n[1])+j,k,f)];                                     
                }
                else if(i>((p->n[0])-3) && j<2)
                {
                     if((p->boundtype[f][1][0])==0)
                  if((i-(p->n[0]))==j)                  
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,4-(p->n[0])+i,j,k,f)];                                    
                  else                  
                   wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,(p->n[1])-4+j,k,f)];                                    
                }
                else if(i>((p->n[0])-3) && j>((p->n[1])-3))
                {
                     if((p->boundtype[f][1][0])==0)
                  if(i==j)                  
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,i,4-(p->n[1])+j,k,f)];                                    
                  else                  
                    wt[encode3_MODID(p,i,j,k,f)]=wt[encode3_MODID(p,4-(p->n[0])+i,j,k,f)];                                    
                }                       
                 
       #endif         




}


__device__ __host__
real sacdabs_MODID(real val) {
   return(fabs(val));
}
