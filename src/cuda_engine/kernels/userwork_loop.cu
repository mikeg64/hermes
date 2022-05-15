//Only for eclipse parsers
#ifdef __CDT_PARSER__
#define __global__
#define __device__
#define __shared__
#define __const__
#endif

// Started as:
// userwork_loop_dev_dev<<<nnBlocks, BLOCK_SIZE>>>(pG->U, is-1, ie+1, js-1, je+1, sizex, hdtodx1, hdtodx2);
__global__ void userwork_loop_dev( Gas *U, int is, int ie, int js, int je, int sizex, Real dtodx2) {
  
  
  Real x1,x2;
  int i = (blockIdx.x * blockDim.x) + threadIdx.x;
  int j = i / sizex;
  i = i % sizex;

  /* Check bounds */
  if(i < is || i > ie || j < js || j > je) return;

  int ind = j*sizex+i;
   
  
  cc_pos_dev(pG_gpu,  i, j, &x1, &x2);
  
  
  		xp=x1-xxmin;
		yp=x3-yymin;
		zp=x2;

		r2=(zp-xcz)*(zp-xcz);
        r1=(xp-xcx)*(xp-xcx);
		
        exp_y=exp(-r1/(delta_x*delta_x));
		exp_z=exp(-r2/(delta_z*delta_z));
        exp_x=exp(-r1/(delta_y*delta_y));

        exp_xyz=exp_x*exp_z;

		vvz=AA*exp_xyz*tdep;
  
  
  		U[ind].M2 += dt*vvz*(U[ind].d);
		U[ind].E += dt*vvz*vvz*(U[ind].d)/2.0;
  
  
  
  

  
}
