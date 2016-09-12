
char **hlines; //header lines for vac config files 
hlines=(char **)calloc(5, sizeof(char*));
// Define time-domain
real dt;


real **d_gw, **d_gwnew;
real **d_gwtemp,**d_gwtemp1,**d_gwtemp2;
real **d_gwmod,  **d_gdwn1,  **d_gwd;


real *d_w;
real *d_wnew;

real *d_wmod,  *d_dwn1,  *d_dwn2,  *d_dwn3,  *d_dwn4,  *d_wd;

real *w,*wnew,*wd,*wdnew, *temp2,*wmod;
real *d_wtemp,*d_wtemp1,*d_wtemp2;

#ifdef USE_MULTIGPU
  real *gmpivisc0,*gmpivisc1,*gmpivisc2, *gmpiw, *gmpiwmod, *gmpiw0, *gmpiwmod0, *gmpiw1, *gmpiwmod1, *gmpiw2, *gmpiwmod2;
  real *d_gmpivisc0,*d_gmpivisc1,*d_gmpivisc2, *d_gmpiw, *d_gmpiwmod, *d_gmpiw0, *d_gmpiwmod0, *d_gmpiw1, *d_gmpiwmod1, *d_gmpiw2,
  *d_gmpiwmod2;
  real *d_gmpiviscr0,*d_gmpiviscr1,*d_gmpiviscr2, *d_gmpiwr0, *d_gmpiwmodr0, *d_gmpiwr1, *d_gmpiwmodr1, *d_gmpiwr2, *d_gmpiwmodr2;

  real **d_ggmpivisc0,**d_ggmpivisc1,**d_ggmpivisc2, **d_ggmpiw, **d_ggmpiwmod, **d_ggmpiw0, **d_ggmpiwmod0, **d_ggmpiw1, **d_ggmpiwmod1, **d_ggmpiw2, **d_ggmpiwmod2;
#endif


#ifdef USE_MULTIGPU
//buffers to use on GPU
  real *d_gmpisendbuffer;
  real *d_gmpirecvbuffer;

   
  real *d_gmpisrcbufferl;
  real *d_gmpisrcbufferr;
  real *d_gmpitgtbufferl;
  real *d_gmpitgtbufferr;
#endif

