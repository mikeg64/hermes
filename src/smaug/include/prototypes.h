
#include "iotypes.h"
#include "iobparams.h"




int cuinit(struct params **p, struct bparams **bp, real **wmod,real **wnew, real **wd, struct state **state, struct params **d_p, struct bparams **d_bp, real **d_wnew, real **d_wmod, real **d_dwn1, real **d_wd, struct state **d_state, real **d_wtemp, real **d_wtemp1, real **d_wtemp2);
int cuupdatemod(struct params **p, struct bparams **bp,real **w, real **wnew, real **wd, struct state **state, struct params **d_p, struct bparams **d_bp,real **d_w, real **d_wnew, real **d_wmod, real **d_dwn1, real **d_wd, struct state **d_state, real **d_wtemp, real **d_wtemp1, real **d_wtemp2);
int initgrid(struct params **p,   struct state **state, real **wd, struct params **d_p, real **d_dwn1, real **d_wd, struct state **d_state, real **d_wtemp, real **d_wtemp1, real **d_wtemp2);

int cufinish(struct params **p, real **w, real **wnew,   struct state **state, struct params **d_p, struct bparams **d_bp,real **d_w, real **d_wnew, real **d_wmod, real **d_dwn1, real **d_wd, struct state **d_state, real **d_wtemp, real **d_wtemp1, real **d_wtemp2);

int cuboundary(struct params **p, struct bparams **bp,struct params **d_p, struct bparams **d_bp,struct state **d_s,  real **d_wmod,  int order,int idir,int field);
int cuupdate(struct params **p, real **w, real **wmod,real **wtemp2,  struct state **state,struct params **d_p, real **d_w,  real **d_wmod,  real **d_wtemp2,  struct state **d_state,int step);

int cuupdatehostwd(struct params **p, real **wd, real **wmod,real **wtemp2,  struct state **state,struct params **d_p, real **d_wd,  real **d_wmod,  real **d_wtemp2,  struct state **d_state,int step);
int cuupdatedevicewd(struct params **p, real **wd, real **wmod,real **wtemp2,  struct state **state,struct params **d_p, real **d_wd,  real **d_wmod,  real **d_wtemp2,  struct state **d_state,int step);



int cucentdiff1(struct params **p, struct params **d_p, struct state **d_s,real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero, real dt, int field, int dir);
int cucentdiff2(struct params **p, struct params **d_p, struct state **d_s, real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd, int order,int ordero, real dt, int field,int dir);

int cugrav(struct params **p, struct params **d_p, struct state **d_s, real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd, int order,int ordero, real dt);
int cusource(struct params **p, struct params **d_p, struct state **d_s, real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd, int order,int ordero, real dt);

int cuadvance(struct params **p, struct params **d_p,    real **d_wmod, real **d_w,int order);
int cucomputedervfields(struct params **p, struct params **d_p,  real **d_wmod, real **d_wd, int order);
int cucomputevels(struct params **p, struct params **d_p,  real **d_wmod, real **d_wd, int order,int dir);
int cucomputemaxc(struct params **p, struct params **d_p,  real **d_wmod, real **d_wd, int order,int dir, real **wd, real **d_wtemp);
int cucomputemaxcourant(struct params **p, struct params **d_p,  real **d_wmod, real **d_wd, int order,int dir, real **wd, real **d_wtemp);
int cucomputec(struct params **p, struct params **d_p,  real **d_wmod, real **d_wd, int order,int dir);
int cucomputept(struct params **p, struct params **d_p,  real **d_wmod, real **d_wd, int order,int dir);
int cucomputepk(struct params **p, struct params **d_p,  real **d_wmod, real **d_wd, int order,int dir);
int cucomputepbg(struct params **p, struct params **d_p,  real **d_wmod, real **d_wd, int order,int dir);
int cudivb(struct params **p,struct params **d_p, real **d_w,  real **d_wmod, real **d_dwn1, real **d_wd, int order,int ordero, real dt);

int cuhyperdifmomsource1(struct params **p,  struct params **d_p,   real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero, real **d_wtemp, int field, int dim, int ii, int ii0, real dt);

int cuhyperdifmomsourcene1(struct params **p,  struct params **d_p,   real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero, real **d_wtemp, int field, int dim, int ii, int ii0, real dt);

int cuhyperdifesource1(struct params **p,  struct params **d_p,   real **d_wmod, real **d_dwn1, real **d_wd, int order,int ordero, real **d_wtemp, int field, int dim, real dt);

int cuhyperdifbsource1(struct params **p,  struct params **d_p,   real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero, real **d_wtemp, int field, int dim, int jj, int ii0,int mm,real sb,real dt);

int cuhyperdifbsourcene1(struct params **p, struct params **d_p,   real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero, real **d_wtemp, int field, int dim, int jj, int ii0,int mm,real sb,real dt);


int cuhyperdifrhosource1(struct params **p,  struct params **d_p,   real **d_wmod, real **d_dwn1, real **d_wd, int order, int ordero, real **d_wtemp, int field, int dim, real dt);



int cunushk1(struct params **p,  struct params **d_p,   real **d_wmod,  real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2);
int cugetdtvisc1(struct params **p,  struct params **d_p,   real **d_wmod, real **wd, real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2);

int cuhyperdifvisc1(struct params **p,  struct params **d_p,   real **d_wmod,  real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2, int field, int dim,int hand);
int cuhyperdifvisc1r(struct params **p,  struct params **d_p,   real **d_wmod,real **wd,  real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2, int field, int dim);
int cuhyperdifvisc1l(struct params **p,  struct params **d_p,   real **d_wmod, real **wd, real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2, int field, int dim);
int cuhyperdifvisc1ir(struct params **p,  struct params **d_p,   real **d_wmod,  real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2, int field, int dim);
int cuhyperdifvisc1il(struct params **p,  struct params **d_p,   real **d_wmod,  real **d_wd, int order, real **d_wtemp, real **d_wtemp1, real **d_wtemp2, int field, int dim);


int cuhyperdifviscmax(struct params **p, real **w, real **wnew, struct params **d_p, real **d_w, real **d_wnew,  real **d_wmod, real **d_dwn1, real **d_wd, int order, real **d_wtemp, int field, int dim);
int cusync(struct params **p);
#ifdef USE_MULTIGPU

int cuinitmgpurbuffers(struct params **p,    
real **d_gmpiviscr0,    
real **d_gmpiviscr1,    
real **d_gmpiviscr2,   
real **d_gmpiwr0, 
real **d_gmpiwmodr0,   
real **d_gmpiwr1, 
real **d_gmpiwmodr1,   
real **d_gmpiwr2, 
real **d_gmpiwmodr2);

 int cuinitmgpubuffers(struct params **p,real **w, real **wmod, real **temp2, real **gmpivisc0, real **gmpivisc1, real **gmpivisc2,   real **gmpiw0, real **gmpiwmod0,   real **gmpiw1, real **gmpiwmod1,   real **gmpiw2, real **gmpiwmod2, struct params **d_p,   real **d_w, real **d_wmod,real **d_wtemp2,    real **d_gmpivisc0,    real **d_gmpivisc1,    real **d_gmpivisc2,   real **d_gmpiw0, real **d_gmpiwmod0,   real **d_gmpiw1, real **d_gmpiwmod1,   real **d_gmpiw2, real **d_gmpiwmod2);
 int cufinishmgpu(struct params **p,real **w, real **wmod, real **temp2, real **gmpivisc0, real **gmpivisc1, real **gmpivisc2,   real **gmpiw0, real **gmpiwmod0,   real **gmpiw1, real **gmpiwmod1,   real **gmpiw2, real **gmpiwmod2, struct params **d_p,   real **d_w, real **d_wmod,real **d_wtemp2,    real **d_gmpivisc0,    real **d_gmpivisc1,    real **d_gmpivisc2,   real **d_gmpiw0, real **d_gmpiwmod0,   real **d_gmpiw1, real **d_gmpiwmod1,   real **d_gmpiw2, real **d_gmpiwmod2);
int cucopywtompiw(struct params **p,real **w, real **wmod,    real **gmpiw0, real **gmpiwmod0,    real **gmpiw1, real **gmpiwmod1,    real **gmpiw2, real **gmpiwmod2, struct params **d_p  ,real **d_w, real **d_wmod,   real **d_gmpiw0, real **d_gmpiwmod0,   real **d_gmpiw1, real **d_gmpiwmod1,   real **d_gmpiw2, real **d_gmpiwmod2, int order, int idir);
int cucopywfrommpiw(struct params **p,real **w, real **wmod,    real **gmpiw0, real **gmpiwmod0,    real **gmpiw1, real **gmpiwmod1,    real **gmpiw2, real **gmpiwmod2, struct params **d_p  ,real **d_w, real **d_wmod,   real **d_gmpiw0, real **d_gmpiwmod0,   real **d_gmpiw1, real **d_gmpiwmod1,   real **d_gmpiw2, real **d_gmpiwmod2, int order, int idir);

int cucopywdtompiwd(struct params **p,real **wd,     real **gmpiw0,    real **gmpiw1,    real **gmpiw2, struct params **d_p  ,real **d_wd,    real **d_gmpiw0,   real **d_gmpiw1,   real **d_gmpiw2,  int order,int idir);

int cucopywdfrommpiwd(struct params **p,real **wd,     real **gmpiw0,    real **gmpiw1,    real **gmpiw2, struct params **d_p  ,real **d_wd,    real **d_gmpiw0,   real **d_gmpiw1,   real **d_gmpiw2,  int order, int idir);


int cucopytompivisc(struct params **p,real **temp2, real **gmpivisc0, real **gmpivisc1, real **gmpivisc2,  struct params **d_p,real **d_wtemp2,    real **d_gmpivisc0,    real **d_gmpivisc1,    real **d_gmpivisc2);
int cucopyfrommpivisc(struct params **p,real **temp2, real **gmpivisc0, real **gmpivisc1, real **gmpivisc2,  struct params **d_p,real **d_wtemp2,    real **d_gmpivisc0,    real **d_gmpivisc1,    real **d_gmpivisc2);
int cusetgpu(struct params **p);



int cucopywtompiwmod(struct params **p,real **w, real **wmod,    real **gmpiw0, real **gmpiwmod0,    real **gmpiw1, real **gmpiwmod1,    real **gmpiw2, real **gmpiwmod2, struct params **d_p  ,real **d_w, real **d_wmod,   real **d_gmpiw0, real **d_gmpiwmod0,   real **d_gmpiw1, real **d_gmpiwmod1,   real **d_gmpiw2, real **d_gmpiwmod2, int order, int idir);

int cucopywmodfrommpiw(struct params **p,real **w, real **wmod,    real **gmpiw0, real **gmpiwmod0,    real **gmpiw1, real **gmpiwmod1,    real **gmpiw2, real **gmpiwmod2, struct params **d_p  ,real **d_w, real **d_wmod,   real **d_gmpiw0, real **d_gmpiwmod0,   real **d_gmpiw1, real **d_gmpiwmod1,   real **d_gmpiw2, real **d_gmpiwmod2, int order, int idir);


#endif
