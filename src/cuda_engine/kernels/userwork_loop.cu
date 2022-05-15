//Only for eclipse parsers
#ifdef __CDT_PARSER__
#define __global__
#define __device__
#define __shared__
#define __const__
#endif

// Started as:
// userwork_loop_dev_dev<<<nnBlocks, BLOCK_SIZE>>>(pG->U, is-1, ie+1, js-1, je+1, sizex, hdtodx1, hdtodx2);
__global__ void userwork_loop_dev(Grid_gpu *pG_gpu, Gas *U, int is, int ie, int js, int je, Real tdep, Real xxmax, Real yymax, Real xxmin, Real yymin, Real AA, Real delta_x, Real delta_y, Real xcz, Real xcx, Real delta_x, Real delta_y) {
  
  
  Real x1,x2;
  int i = (blockIdx.x * blockDim.x) + threadIdx.x;
  int j = i / sizex;
  i = i % sizex;

  /* Check bounds */
  if(i < is || i > ie || j < js || j > je) return;

  int ind = j*sizex+i;
   
  
  cc_pos_dev(pG_gpu,  i, j, &x1, &x2);
  
  
  		xp=x1-xxmin;
		yp=x2-yymin;


		r2=(yp-xcz)*(yp-xcz);
        r1=(xp-xcx)*(xp-xcx);
		
 
		exp_z=exp(-r2/(delta_y*delta_y));
        exp_x=exp(-r1/(delta_x*delta_x));

        exp_xyz=exp_x*exp_z;

		vvz=AA*exp_xyz*tdep;
  
  
  		U[ind].M2 += dt*vvz*(U[ind].d);
		U[ind].E += dt*vvz*vvz*(U[ind].d)/2.0;
  
  
  
  

  
}
