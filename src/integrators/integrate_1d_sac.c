#include "../copyright.h"


/*============================================================================*/
/*! \file integrate_1d_sac.c
 *  \Compute fluxes using sac . 
 *
 * PURPOSE: Integrate MHD equations using 3D version of the directionally
 *   .  The variables updated are:
 *   -  U.[d,M1,M2,E,B1c,B2c,s] -- where U is of type ConsS
 *   -  B1i, B2i  -- interface magnetic field
 *   Also adds gravitational source terms.
 *   - For adb hydro, requires (9*ConsS +  3*Real) = 48 3D arrays
 *   - For adb mhd, requires   (9*ConsS + 10*Real) = 73 3D arrays
 *   
 *  Source terms added are hyperdiffusion terms
 *
 * REFERENCES:
 * - Magnetohydrodynamic code for gravitationally-stratified media : SAC
 *   Shelyag et al 
 *   http://adsabs.harvard.edu/abs/2008A%26A...486..655S
 *
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - integrate_1d_sac()
 * - integrate_init_1d()
 * - integrate_destruct_1d() */

/*
 * PROGRESS
 * Initial boiler plate based on athena integrate_1d_ctu.c

 */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"
#ifdef PARTICLES
#include "../particles/particle.h"
#endif

#if defined(SAC_INTEGRATOR)
#ifdef SPECIAL_RELATIVITY
#error : The SAC integrator cannot be used for special relativity.
#endif /* SPECIAL_RELATIVITY */

/* The L/R states of conserved variables and fluxes at each cell face */
static Cons1DS *Uc_x1=NULL, *Ur_x1Face=NULL, *x1Flux=NULL;
static ConsS ***Uinit=NULL; /*Uinit used to store initial fields*/
/* 1D scratch vectors used by lr_states and flux functions */
static Real *Bx=NULL, *Bxb=NULL, *Bxc=NULL, *temp=NULL, *grad=NULL;
static Prim1DS *W=NULL, *Wl=NULL, *Wr=NULL;
static Cons1DS *U1d=NULL;

/* Variables at t^{n+1/2} used for source terms */
static Real *dhalf = NULL, *phalf = NULL;

/* Variables needed for cylindrical coordinates */
#ifdef CYLINDRICAL
static Real *geom_src=NULL;
#endif



static void computemaxc(ConsS ***Uint, GridS *pG, int dim);


static void hyperdifviscr(int fieldi,int dim,ConsS ***Uint, GridS *pG);
static void hyperdifviscl(int fieldi,int dim,ConsS ***Uint, GridS *pG);

static void hyperdifrhosource(int dim,Real dt,ConsS ***Uint, GridS *pG);



static void hyperdifesource(int dim,Real dt,ConsS ***Uint, GridS *pG);
//static void hyperdifmomsource(int k,int l, int ii,int ii0,Real dt,ConsS ***Uint, GridS *pG);
//static void hyperdifmomsourcene(int k,int l, int ii,int ii0, Real dt,ConsS ***Uint, GridS *pG);

static void hyperdifmomsource(int kf,int nf,int ii,int ii0,Real dt,ConsS ***Uint, GridS *pG);
static void hyperdifmomsourcene(int kf,int nf, int ii,int ii0, Real dt,ConsS ***Uint, GridS *pG);

#ifdef MHD
static void hyperdifbsource(int kf,int lf,int jj,int ii0, int mm, Real sb,  Real dt, ConsS ***Uint, GridS *pG);
static void hyperdifbsourcene(int kf,int lf,int jj,int ii0, int mm, Real sb,  Real dt, ConsS ***Uint, GridS *pG);
#endif /* MHD */



/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void integrate_1d_sac(DomainS *pD)
 *  \brief 1D version of CTU unsplit integrator for MHD
 *
 *   The numbering of steps follows the numbering in the 3D version.
 *   NOT ALL STEPS ARE NEEDED IN 1D.
 */

void integrate_1d_sac(DomainS *pD)
{
  GridS *pG=(pD->Grid);
  Real dtodx1 = pG->dt/pG->dx1, hdtodx1 = 0.5*pG->dt/pG->dx1;
  int i,il,iu, is = pG->is, ie = pG->ie;
  int js = pG->js;
  int ks = pG->ks;
  Real x1,x2,x3,phicl,phicr,phifc,phil,phir,phic,M1h,M2h,M3h;
#ifndef BAROTROPIC
  Real coolfl,coolfr,coolfc,coolf,Eh=0.0;
#endif
#if defined(MHD) 
  Real B1ch,B2ch,B3ch;
#endif
#if (NSCALARS > 0)
  int n;
#endif
#ifdef SELF_GRAVITY
  Real gxl,gxr,flux_m1l,flux_m1r;
#endif
#ifdef FEEDBACK
  Real dt1 = 1.0/pG->dt;
#endif
#ifdef STATIC_MESH_REFINEMENT
  int ncg,npg,dim;
#endif

#ifdef CYLINDRICAL
#ifndef ISOTHERMAL
  Real Pavgh;
#endif
  Real g,gl,gr,rinv;
  Real hdt = 0.5*pG->dt;
  Real geom_src_d,geom_src_Vx,geom_src_Vy,geom_src_P,geom_src_By,geom_src_Bz;
  const Real *r=pG->r, *ri=pG->ri;
#endif /* CYLINDRICAL */
  Real lsf=1.0, rsf=1.0;

/* With particles, one more ghost cell must be updated in predict step */
#ifdef PARTICLES
  Real d1;
  il = is - 3;
  iu = ie + 3;
#else
  il = is - 1;
  iu = ie + 1;
#endif

/* Compute predictor feedback from particle drag */
#ifdef FEEDBACK
  feedback_predictor(pD);
  exchange_gpcouple(pD,1);
#endif


/*Used for hyperdiffusion computations*/
int ii1, dim, ii, ii0, mm;
int jj1, jj, jj0;
int fieldi; /*integers map to following index rho, mom1, mom2, energy, b1, b2,energyb,rhob,b1b,b2b*/
Real sb;
int size1;

size1=1+ie+2*nghost-is;


/*=== STEP 1: Compute L/R x1-interface states and 1D x1-Fluxes ===============*/
printf("step1 \n");
ks=0;
js=0;

/*Store the initial variables for use in hyperdiffusion computations*/
      for (i=il; i<=iu; i++) {

        Uinit[0][0][i].d  = pG->U[ks][js][i].d;
        Uinit[0][0][i].M1 = pG->U[ks][js][i].M1;
        Uinit[0][0][i].M2 = pG->U[ks][js][i].M2;
        Uinit[0][0][i].M3 = pG->U[ks][js][i].M3;
#ifndef BAROTROPIC
        Uinit[0][0][i].E  = pG->U[ks][js][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
        Uinit[0][0][i].B1c = pG->U[ks][js][i].B1c;
        Uinit[0][0][i].B2c = pG->U[ks][js][i].B2c;
        Uinit[0][0][i].B3c = pG->U[ks][js][i].B3c;
        //Bxc[i] = pG->U[ks][js][i].B1c;
        //Bxb[i] = pG->B1cb[k][j][i];
        //B1_x1Face[k][j][i] = pG->B1i[k][j][i];
#endif /* MHD */


#ifdef SAC_INTEGRATOR
        Uinit[0][0][i].db  = pG->U[ks][js][i].db;

#ifdef MHD
        Uinit[0][0][i].B1cb = pG->U[ks][js][i].B1cb;        
        Uinit[0][0][i].B2cb = pG->U[ks][js][i].B2cb;
        Uinit[0][0][i].B3cb = pG->U[ks][js][i].B3cb;
#endif
#endif


#ifdef SMAUG_INTEGRATOR
        Uinit[0][0][i].db  = pG->U[ks][js][i].db;

#ifdef MHD
        Uinit[0][0][i].B1cb = pG->U[ks][js][i].B1cb;        
        Uinit[0][0][i].B2cb = pG->U[ks][js][i].B2cb;
        Uinit[0][0][i].B3cb = pG->U[ks][js][i].B3cb;
#endif
#endif


#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) Uinit[0][0][i].s[n] = pG->U[ks][js][i].s[n];
#endif


}



printf("step1a\n");

/*--- Step 1a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M1, M2, M3, E, B2c, B3c, s[n])
 */

  for (i=is-nghost; i<=ie+nghost; i++) {
    U1d[i].d  = pG->U[ks][js][i].d;
    U1d[i].Mx = pG->U[ks][js][i].M1;
    U1d[i].My = pG->U[ks][js][i].M2;
    U1d[i].Mz = pG->U[ks][js][i].M3;
#ifndef BAROTROPIC
    U1d[i].E  = pG->U[ks][js][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
    U1d[i].By = pG->U[ks][js][i].B2c;
    U1d[i].Bz = pG->U[ks][js][i].B3c;
    Bxc[i] = pG->U[ks][js][i].B1c;
    /*Bxi[i] = pG->B1i[ks][js][i];*/
#endif /* MHD */
#ifdef SAC_INTEGRATOR
        U1d[i].db  = pG->U[ks][js][i].db;
        
#ifdef MHD
        U1d[i].Byb = pG->U[ks][js][i].B2cb;
        U1d[i].Bzb = pG->U[ks][js][i].B3cb;
        Bxb[i] = pG->U[ks][js][i].B1cb;
#endif /* MHD */

#endif
#ifdef SMAUG_INTEGRATOR
        U1d[i].db  = pG->U[ks][js][i].db;
        
#ifdef MHD
        U1d[i].Byb = pG->U[ks][js][i].B2cb;
        U1d[i].Bzb = pG->U[ks][js][i].B3cb;
        Bxb[i] = pG->U[ks][js][i].B1cb;
#endif /* MHD */


//printf("here\n");
#ifdef MHD
 Bxb[i] =0.0; //temporary debug seg fault
#endif

#endif
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) U1d[i].s[n] = pG->U[ks][js][i].s[n];
#endif
  }

 

/*--- Step 1b ------------------------------------------------------------------
 * Compute L and R states at X1-interfaces, add "MHD source terms" for 0.5*dt
 */


  for (i=is-nghost; i<=ie+nghost; i++) {
    W[i] = Cons1D_to_Prim1D(&U1d[i], &Bxc[i],&Bxb[i]);
  }

 /* lr_states(pG,W,Bxc,pG->dt,pG->dx1,il+1,iu-1,Wl,Wr,1);*/

/* Apply density floor */
 /* for (i=il+1; i<=iu; i++){
    if (Wl[i].d < d_MIN) {
      Wl[i].d = d_MIN;
    }
    if (Wr[i].d < d_MIN) {
      Wr[i].d = d_MIN;
    }
   }  */






 
printf("step1c\n");
/*--- Step 1c ------------------------------------------------------------------
 * Add source terms from static gravitational potential for 0.5*dt to L/R states
 */

  if (StaticGravPot != NULL){
    for (i=il+1; i<=iu; i++) {
      cc_pos(pG,i,js,ks,&x1,&x2,&x3);
// #ifdef CYLINDRICAL
//       gl = (*x1GravAcc)(x1vc(pG,i-1),x2,x3);
//       gr = (*x1GravAcc)(x1vc(pG,i),x2,x3);
//       gl = (*x1GravAcc)(x1-pG->dx1,x2,x3);
//       gr = (*x1GravAcc)(x1,x2,x3);
      /* APPLY GRAV. SOURCE TERMS TO V1 USING ACCELERATION FOR (dt/2) */
      // Wl[i].Vx -= hdt*gl;
      // Wr[i].Vx -= hdt*gr;
//      W[i].Vx -= 0.5*hdt*(gr+gl);
// #else
      phicr = (*StaticGravPot)( x1             ,x2,x3);
      phicl = (*StaticGravPot)((x1-    pG->dx1),x2,x3);
//      phifc = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

      W[i].Vx -= dtodx1*(phicr - phicl);
      //Wr[i].Vx -= dtodx1*(phicr - phifc);
// #endif /* CYLINDRICAL */
    }
  }




/*--- Step 1c (cont) -----------------------------------------------------------
 * Add source terms for self-gravity for 0.5*dt to L/R states
 */

#ifdef SELF_GRAVITY
  for (i=il+1; i<=iu; i++) {
    //Wl[i].Vx -= hdtodx1*(pG->Phi[ks][js][i] - pG->Phi[ks][js][i-1]);
    //Wr[i].Vx -= hdtodx1*(pG->Phi[ks][js][i] - pG->Phi[ks][js][i-1]);
    W[i].Vx -= hdtodx1*(pG->Phi[ks][js][i] - pG->Phi[ks][js][i-1]);
  }
#endif






/*--- Step 1c (cont) -----------------------------------------------------------
 * Add source terms from optically-thin cooling for 0.5*dt to L/R states
 */

#ifndef BAROTROPIC
  if (CoolingFunc != NULL){
    for (i=il+1; i<=iu; i++) {
     /* coolfl = (*CoolingFunc)(Wl[i].d,Wl[i].P,(0.5*pG->dt));
      coolfr = (*CoolingFunc)(Wr[i].d,Wr[i].P,(0.5*pG->dt));*/

      /*check cooling function*/
      coolfc = (*CoolingFunc)(W[i].d+W[i].db,W[i].P,(pG->dt));
      W[i].P -= pG->dt*Gamma_1*coolfc;
      /*Wl[i].P -= 0.5*pG->dt*Gamma_1*coolfl;
      Wr[i].P -= 0.5*pG->dt*Gamma_1*coolfr;*/
    }
  }
#endif /* BAROTROPIC */



/*--- Step 1c (cont) -----------------------------------------------------------
 * Add source terms for particle feedback for 0.5*dt to L/R states
 */

#ifdef FEEDBACK
    for (i=il+1; i<=iu; i++) {

      d1 = 1.0/(W[i].d+W[i].db);
      W[i].Vx -= pG->Coup[ks][js][i].fb1*d1;
      W[i].Vy -= pG->Coup[ks][js][i].fb2*d1;
      W[i].Vz -= pG->Coup[ks][js][i].fb3*d1;

#ifndef BAROTROPIC
      //Wl[i].P += pG->Coup[ks][js][i-1].Eloss*Gamma_1;
      //Wr[i].P += pG->Coup[ks][js][i].Eloss*Gamma_1;
      W[i].P += pG->Coup[ks][js][i].Eloss*Gamma_1;
#endif
    }
#endif /* FEEDBACK */




/*--- Step 1c (cont) -----------------------------------------------------------
 * Add the geometric source-terms now using cell-centered primitive
 * variables at time t^n
 */

#ifdef CYLINDRICAL
      for (i=il+1; i<=iu; i++) {
        // left state geometric source term (uses W[i-1])
//         rinv = 1.0/x1vc(pG,i-1);
        rinv = 1.0/r[i-1];
        geom_src_d  = -(W[i-1].d+W[i-1].db)*W[i-1].Vx*rinv;
        geom_src_Vx =  SQR(W[i-1].Vy);
        geom_src_Vy = -W[i-1].Vx*W[i-1].Vy;
#ifdef MHD
        geom_src_Vx -= SQR(W[i-1].By)/(W[i-1].d+W[i-1].db);
        geom_src_Vy += Bxc[i-1]*W[i-1].By/(W[i-1].d+W[i-1].db);
        geom_src_By =  -W[i-1].Vy*Bxc[i-1]*rinv;
        geom_src_Bz =  -W[i-1].Vx*W[i-1].Bz*rinv;
#endif /* MHD */
        geom_src_Vx *= rinv;
        geom_src_Vy *= rinv;
#ifndef ISOTHERMAL
        geom_src_P  = -Gamma*W[i-1].P*W[i-1].Vx*rinv;
#endif /* ISOTHERMAL */

        // add source term to left state
        //Wl[i].d  += hdt*geom_src_d;
        //Wl[i].Vx += hdt*geom_src_Vx;
        //Wl[i].Vy += hdt*geom_src_Vy;
#ifdef MHD
        //Wl[i].By += hdt*geom_src_By;
        //Wl[i].Bz += hdt*geom_src_Bz;
#endif /* MHD */
#ifndef ISOTHERMAL
        //Wl[i].P  += hdt*geom_src_P;
#endif /* ISOTHERMAL */

        // right state geometric source term (uses W[i])
//         rinv = 1.0/x1vc(pG,i);
        rinv = 1.0/r[i];
        geom_src_d  += -(W[i].d+W[i].db)*W[i].Vx*rinv;
        geom_src_Vx +=  SQR(W[i].Vy);
        geom_src_Vy += -W[i].Vx*W[i].Vy;
#ifdef MHD
        geom_src_Vx -= SQR(W[i].By)/(W[i].d+W[i].db);
        geom_src_Vy += Bxc[i]*W[i].By/(W[i].d+W[i].db);
        geom_src_By +=  -W[i].Vy*Bxc[i]*rinv;
        geom_src_Bz +=  -W[i].Vx*W[i].Bz*rinv;
#endif /* MHD */
        geom_src_Vx *= rinv;
        geom_src_Vy *= rinv;
#ifndef ISOTHERMAL
        geom_src_P  += -Gamma*W[i].P*W[i].Vx*rinv;
#endif /* ISOTHERMAL */

        // add source term to right state
        W[i].d  += hdt*geom_src_d;
        W[i].Vx += hdt*geom_src_Vx;
        W[i].Vy += hdt*geom_src_Vy;
#ifdef MHD
        W[i].By += hdt*geom_src_By;
        W[i].Bz += hdt*geom_src_Bz;
#endif /* MHD */
#ifndef ISOTHERMAL
        W[i].P  += hdt*geom_src_P;
#endif /* ISOTHERMAL */
      }
#endif /* CYLINDRICAL */




#ifdef MHD
 Bxb[i] =0.0; //temporary debug seg fault
#endif


printf("step 1d \n");

/*--- Step 1d ------------------------------------------------------------------
 * Compute 1D fluxes in x1-direction, storing into 3D array
 */
    for (i=il+1; i<=iu; i++) {

#ifdef MHD
      Uc_x1[i] = Prim1D_to_Cons1D(&W[i],&Bxc[i],&Bxb[i]);
      fluxes(Uc_x1[i],Uc_x1[i],W[i],W[i],Bxc[i],Bxb[i],&x1Flux[i]);
#else
      Uc_x1[i] = Prim1D_to_Cons1D(&W[i],NULL,NULL);
      fluxes(Uc_x1[i],Uc_x1[i],W[i],W[i],0,0,&x1Flux[i]);
#endif

    }
  

printf("step 8\n");
/*checked for sac to here*/
/*=== STEP 8: Compute cell-centered values at n+1/2 ==========================*/

/*--- Step 8a ------------------------------------------------------------------
 * Calculate d^{n+1/2} (needed with static potential, or cooling)
 */

#ifndef PARTICLES
  if ((StaticGravPot != NULL) || (CoolingFunc != NULL))
#endif
  {
    for (i=il+1; i<=iu-1; i++) {
      dhalf[i] = pG->U[ks][js][i].d - hdtodx1*(x1Flux[i+1].d - x1Flux[i].d );
      if ((dhalf[i] < d_MIN) || (dhalf[i] != dhalf[i])) {
        dhalf[i] = d_MIN;
      }
#ifdef PARTICLES
      pG->Coup[ks][js][i].grid_d = dhalf[i];
#endif
    }
  }

/*--- Step 8b ------------------------------------------------------------------
 * Calculate P^{n+1/2} (needed with cooling)
 */

#ifndef PARTICLES
  if (CoolingFunc != NULL)
#endif /* PARTICLES */
  {
    for (i=il+1; i<=iu-1; i++) {
      M1h = pG->U[ks][js][i].M1 - hdtodx1*(x1Flux[i+1].Mx - x1Flux[i].Mx);
      M2h = pG->U[ks][js][i].M2 - hdtodx1*(x1Flux[i+1].My - x1Flux[i].My);
      M3h = pG->U[ks][js][i].M3 - hdtodx1*(x1Flux[i+1].Mz - x1Flux[i].Mz);
#ifndef BAROTROPIC
      Eh  = pG->U[ks][js][i].E  - hdtodx1*(x1Flux[i+1].E  - x1Flux[i].E );
#endif

/* Add source terms for fixed gravitational potential */
      if (StaticGravPot != NULL){
        cc_pos(pG,i,js,ks,&x1,&x2,&x3);
        phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);
        M1h -= hdtodx1*(phir-phil)*pG->U[ks][js][i].d;
      }

/* Add source terms due to self-gravity  */
#ifdef SELF_GRAVITY
      phir = 0.5*(pG->Phi[ks][js][i] + pG->Phi[ks][js][i+1]);
      phil = 0.5*(pG->Phi[ks][js][i] + pG->Phi[ks][js][i-1]);
      M1h -= hdtodx1*(phir-phil)*(pG->U[ks][js][i].d+pG->U[ks][js][i].db);
#endif /* SELF_GRAVITY */

/* Add the particle feedback terms */
#ifdef FEEDBACK
      M1h -= pG->Coup[ks][js][i].fb1;
      M2h -= pG->Coup[ks][js][i].fb2;
      M3h -= pG->Coup[ks][js][i].fb3;
#endif /* FEEDBACK */

#ifndef BAROTROPIC
      phalf[i] = Eh - 0.5*(M1h*M1h + M2h*M2h + M3h*M3h)/dhalf[i];

#ifdef MHD
      B1ch = pG->U[ks][js][i].B1c;
      B2ch = pG->U[ks][js][i].B2c - hdtodx1*(x1Flux[i+1].By - x1Flux[i].By);
      B3ch = pG->U[ks][js][i].B3c - hdtodx1*(x1Flux[i+1].Bz - x1Flux[i].Bz);
      phalf[i] -= 0.5*(B1ch*B1ch + B2ch*B2ch + B3ch*B3ch);
#endif /* MHD */

      phalf[i] *= Gamma_1;
#endif /* BAROTROPIC */

#ifdef PARTICLES
      d1 = 1.0/dhalf[i];
      pG->Coup[ks][js][i].grid_v1 = M1h*d1;
      pG->Coup[ks][js][i].grid_v2 = M2h*d1;
      pG->Coup[ks][js][i].grid_v3 = M3h*d1;
#ifndef BAROTROPIC
      pG->Coup[ks][js][i].grid_cs = sqrt(Gamma*phalf[i]*d1);
#endif  /* BAROTROPIC */
#endif /* PARTICLES */

    }
  }

/*=== STEP 8.5: Integrate the particles, compute the feedback ================*/

#ifdef PARTICLES
  Integrate_Particles(pD);
#ifdef FEEDBACK
  exchange_gpcouple(pD,2);
#endif
#endif

/*=== STEPS 9-10: Not needed in 1D ===*/

/*=== STEP 11: Add source terms for a full timestep using n+1/2 states =======*/

/*--- Step 11a -----------------------------------------------------------------
 * Add geometric source terms
 */

#ifdef CYLINDRICAL
  for (i=il+1; i<=iu-1; i++) {
    rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];

    /* calculate density at time n+1/2 */
    dhalf[i] = pG->U[ks][js][i].d
             - hdtodx1*(rsf*x1Flux[i+1].d - lsf*x1Flux[i].d);

    /* calculate x2-momentum at time n+1/2 */
    M2h = pG->U[ks][js][i].M2
        - hdtodx1*(SQR(rsf)*x1Flux[i+1].My - SQR(lsf)*x1Flux[i].My);

    /* compute geometric source term at time n+1/2 */
    geom_src[i] = SQR(M2h)/dhalf[i];
#ifdef MHD
    B2ch = pG->U[ks][js][i].B2c - hdtodx1*(x1Flux[i+1].By - x1Flux[i].By);
    geom_src[i] -= SQR(B2ch);
#endif
#ifdef ISOTHERMAL
    geom_src[i] += Iso_csound2*dhalf[i];
#ifdef MHD
    B1ch = pG->U[ks][js][i].B1c;
    B3ch = pG->U[ks][js][i].B3c - hdtodx1*(rsf*x1Flux[i+1].Bz - lsf*x1Flux[i].Bz);
    geom_src[i] += 0.5*(SQR(B1ch)+SQR(B2ch)+SQR(B3ch));
#endif
#else /* ISOTHERMAL */
    Pavgh = 0.5*(lsf*x1Flux[i].Pflux + rsf*x1Flux[i+1].Pflux);
    geom_src[i] += Pavgh;
#endif
//     geom_src[i] /= x1vc(pG,i);
    geom_src[i] /= r[i];

    /* add time-centered geometric source term for full dt */
    pG->U[ks][js][i].M1 += pG->dt*geom_src[i];
  }
#endif /* CYLINDRICAL */

/*--- Step 11a -----------------------------------------------------------------
 * Add gravitational source terms as a Static Potential.
 *   The energy source terms computed at cell faces are averaged to improve
 * conservation of total energy.
 *    S_{M} = -(\rho)^{n+1/2} Grad(Phi);   S_{E} = -(\rho v)^{n+1/2} Grad{Phi}
 */

  if (StaticGravPot != NULL){
    for (i=is; i<=ie; i++) {
      cc_pos(pG,i,js,ks,&x1,&x2,&x3);
      phic = (*StaticGravPot)((x1            ),x2,x3);
      phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
      phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

#ifdef CYLINDRICAL
//       g = (*x1GravAcc)(x1vc(pG,i),x2,x3);
      rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
//       pG->U[ks][js][i].M1 -= pG->dt*dhalf[i]*g;
      pG->U[ks][js][i].M1 -= dtodx1*dhalf[i]*(phir-phil);
#else
      pG->U[ks][js][i].M1 -= dtodx1*dhalf[i]*(phir-phil);
#endif

#ifndef BAROTROPIC
      pG->U[ks][js][i].E -= dtodx1*(lsf*x1Flux[i  ].d*(phic - phil) +
                                    rsf*x1Flux[i+1].d*(phir - phic));
#endif
    }
  }

/*--- Step 11b -----------------------------------------------------------------
 * Add source terms for self-gravity.
 * A flux correction using Phi^{n+1} in the main loop is required to make
 * the source terms 2nd order: see selfg_flux_correction().
 */

#ifdef SELF_GRAVITY
  for (i=is; i<=ie; i++) {
      phic = pG->Phi[ks][js][i];
      phil = 0.5*(pG->Phi[ks][js][i-1] + pG->Phi[ks][js][i  ]);
      phir = 0.5*(pG->Phi[ks][js][i  ] + pG->Phi[ks][js][i+1]);

      gxl = (pG->Phi[ks][js][i-1] - pG->Phi[ks][js][i  ])/(pG->dx1);
      gxr = (pG->Phi[ks][js][i  ] - pG->Phi[ks][js][i+1])/(pG->dx1);

/* 1-momentum flux.  2nd term is needed only if Jean's swindle used */
      flux_m1l = 0.5*(gxl*gxl)/four_pi_G + grav_mean_rho*phil;
      flux_m1r = 0.5*(gxr*gxr)/four_pi_G + grav_mean_rho*phir;

      pG->U[ks][js][i].M1 -= dtodx1*(flux_m1r - flux_m1l);
#ifndef BAROTROPIC
      pG->U[ks][js][i].E -= dtodx1*(x1Flux[i  ].d*(phic - phil) +
                                    x1Flux[i+1].d*(phir - phic));
#endif
  }

/* Save mass fluxes in Grid structure for source term correction in main loop */
  for (i=is; i<=ie+1; i++) {
    pG->x1MassFlux[ks][js][i] = x1Flux[i].d;
  }
#endif /* SELF_GRAVITY */

/*--- Step 11c -----------------------------------------------------------------
 * Add source terms for optically thin cooling
 */

#ifndef BAROTROPIC
  if (CoolingFunc != NULL){
    for (i=is; i<=ie; i++) {
      coolf = (*CoolingFunc)(dhalf[i],phalf[i],pG->dt);
      pG->U[ks][js][i].E -= pG->dt*coolf;
    }
  }
#endif /* BAROTROPIC */

printf("step11\n");
/*--- Step 11d -----------------------------------------------------------------
 * Add source terms for particle feedback
 */

#ifdef FEEDBACK
  for (i=is; i<=ie; i++) {
    pG->U[ks][js][i].M1 -= pG->Coup[ks][js][i].fb1;
    pG->U[ks][js][i].M2 -= pG->Coup[ks][js][i].fb2;
    pG->U[ks][js][i].M3 -= pG->Coup[ks][js][i].fb3;
#ifndef BAROTROPIC
    pG->U[ks][js][i].E += pG->Coup[ks][js][i].Eloss;
    pG->Coup[ks][js][i].Eloss *= dt1;
#endif
  }
#endif



printf("start hyperdiffusion\n");
/* compute hyperdiffusion source terms  */

//hyperdifvisc1r

//hyperdifvisc1l

//computec
//computemaxc(Uinit,pG);

//density contribution
for(dim=0; dim<2; dim++) //each direction
{

printf("step12c maxc\n");
;//computemaxc(Uinit,pG,dim);
printf("step12c viscr\n");
;//hyperdifviscr(rho,dim,Uinit, pG);
printf("step12c viscl\n");
;//hyperdifviscl(rho,dim,Uinit, pG);
//hyperdifvisc1ir
//hyperdifvisc1il
//int dim,Real dt,ConsS ***Uint, GridS *pG
printf("step12c rhosource\n");
;//hyperdifrhosource(dim,pG->dt,Uinit, pG) ;
}

//energy hyperdiffusion term
for(dim=0; dim<2; dim++) //each direction
{
//hyperdifvisc1ir
//hyperdifvisc1il
//hyperdifesource1


;//computemaxc(Uinit,pG,dim);

;//hyperdifviscr(energy,dim,Uinit, pG);

;//hyperdifviscl(energy,dim,Uinit, pG);

;//hyperdifesource(dim,pG->dt,Uinit, pG) ;

}



       //momentum hyperdiffusion term
for(dim=0; dim<2; dim++) //each direction   //k
for( fieldi=0; fieldi<2; fieldi++)          //l
{
//hyperdifvisc1ir
//hyperdifvisc1il
//hyperdifesource1
;//hyperdifviscr(mom1+fieldi,dim,Uinit, pG);
;//hyperdifviscl(mom1+fieldi,dim,Uinit, pG);

		         for(ii1=0;ii1<=1;ii1++)
		         {
		                  if (ii1 == 0)
		                  {
				           ii=dim;
				           ii0=fieldi;  //f is field
		                  }
		                  else
		                  {
				           ii=fieldi;
				           ii0=dim;
		                   }

				  if(ii==dim)
				    ;//hyperdifmomsource(dim,fieldi,ii,ii0,pG->dt,Uinit, pG);
				  else
				    ;//hyperdifmomsourcene(dim,fieldi,ii,ii0,pG->dt,Uinit, pG);  //off diagonal
		        }


}
#ifdef MHD

  //b field hyperdiffusion term

for(dim=0; dim<1; dim++) //each direction //k
for( fieldi=0; fieldi<2; fieldi++)          //l
{
//hyperdifvisc1ir
//hyperdifvisc1il

;//hyperdifviscr(b1+fieldi,dim,Uinit, pG);
;//hyperdifviscl(b1+fieldi,dim,Uinit, pG);


if(fieldi != dim)
{



		         for(ii1=0;ii1<=1;ii1++)
		         {

                        if (ii1 == 0)
				          {
						   jj=dim;
						   mm=fieldi;
						   sb=-1.0;
						   ii0=dim;
				          }
				          else
				          {
						   ii0=fieldi;
						   mm=dim;
						   sb=1.0;
						   jj=fieldi;
				          }




				  if(mm==dim)
                    ;//hyperdifbsource(fieldi,dim,jj,ii0,mm,sb,pG->dt,Uinit, pG);
				  else
                    ;//hyperdifbsourcene(fieldi,dim,jj,ii0,mm,sb,pG->dt,Uinit, pG);  //off diagonal

		        }
}


}

#endif  /*hyperdiffusion source term for bfield*/






/*--- Step 12c -----------------------------------------------------------------
 * Update cell-centered variables in pG using 3D x3-Fluxes
 */

/*--- Step 12a -----------------------------------------------------------------
 * Update cell-centered variables in pG using 1D x1-fluxes
 */

printf("flux divergence \n");

/*for (i=is; i<=ie; i++) {*/
/*
Use conserved Uc_x1[i] ConsS
and primitive W[i] Prim1DS

First add divergence of flux
*/

/*compute gradient of velocity*/
for (i=is; i<=ie; i++)
{
   temp[i]=W[i].Vx;
   grad[i]=0;
}
gradient4(temp, size1 ,pG->dx1,grad);

for (i=is+2; i<=ie-2; i++) {

#ifdef CYLINDRICAL
    rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
#endif
    pG->U[ks][js][i].d  -= dtodx1*(lsf*x1Flux[i-2].d+8*rsf*x1Flux[i+1].d  - 8*lsf*x1Flux[i-1].d -rsf*x1Flux[i+2].d)/12;
    pG->U[ks][js][i].M1  -= dtodx1*(lsf*x1Flux[i-2].Mx+8*rsf*x1Flux[i+1].Mx  - 8*lsf*x1Flux[i-1].Mx -rsf*x1Flux[i+2].Mx)/12;
    pG->U[ks][js][i].M2  -= dtodx1*(SQR(lsf)*x1Flux[i-2].My+8*SQR(rsf)*x1Flux[i+1].My  - 8*SQR(lsf)*x1Flux[i-1].My -SQR(rsf)*x1Flux[i+2].My)/12;
    pG->U[ks][js][i].M3  -= dtodx1*(lsf*x1Flux[i-2].Mz+8*rsf*x1Flux[i+1].Mz  - 8*lsf*x1Flux[i-1].Mz -rsf*x1Flux[i+2].Mz)/12;

#ifndef BAROTROPIC
    pG->U[ks][js][i].E  -= dtodx1*(lsf*x1Flux[i-2].E+8*rsf*x1Flux[i+1].E  - 8*lsf*x1Flux[i-1].E -rsf*x1Flux[i+2].E)/12;
    pG->U[ks][js][i].E  -= dtodx1*grad[i]*W[i].Pb;//(pG->dt)*grad[i]*W[i].Pb;

#ifdef MHD
    pG->U[ks][js][i].E  += dtodx1*grad[i]*Bxb[i]*Bxb[i]; //(pG->dt)*grad[i]*Bxb[i]*Bxb[i];
  /*remember energy contribution from background b-field*/
#endif



#endif /* BAROTROPIC */
#ifdef MHD
    pG->U[ks][js][i].B2c  -= dtodx1*(lsf*x1Flux[i-2].By+8*rsf*x1Flux[i+1].By  - 8*lsf*x1Flux[i-1].By -rsf*x1Flux[i+2].By)/12;
    pG->U[ks][js][i].B3c  -= dtodx1*(lsf*x1Flux[i-2].Bz+8*rsf*x1Flux[i+1].Bz  - 8*lsf*x1Flux[i-1].Bz -rsf*x1Flux[i+2].Bz)/12;

/* For consistency, set B2i and B3i to cell-centered values.  */
//    pG->B2i[ks][js][i] = pG->U[ks][js][i].B2c;
//    pG->B3i[ks][js][i] = pG->U[ks][js][i].B3c;
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++)
       pG->U[ks][js][i].s[n]  -= dtodx1*(lsf*x1Flux[i-2].s[n]+8*rsf*x1Flux[i+1].s[n]  - 8*lsf*x1Flux[i-1].s[n] -rsf*x1Flux[i+2].s[n])/12;
#endif



}





 
/*--- Step 12b: Not needed in 1D ---*/
/*--- Step 12c: Not needed in 1D ---*/
/*--- Step 12d: Not needed in 1D ---*/












#ifdef STATIC_MESH_REFINEMENT
/*--- Step 12e -----------------------------------------------------------------
 * With SMR, store fluxes at boundaries of child and parent grids.
 */

/* x1-boundaries of child Grids (interior to THIS Grid) */
  for (ncg=0; ncg<pG->NCGrid; ncg++) {
    for (dim=0; dim<2; dim++){
      if (pG->CGrid[ncg].myFlx[dim] != NULL) {

        if (dim==0) i = pG->CGrid[ncg].ijks[0];
        if (dim==1) i = pG->CGrid[ncg].ijke[0] + 1;

        pG->CGrid[ncg].myFlx[dim][ks][js].d  = x1Flux[i].d; 
        pG->CGrid[ncg].myFlx[dim][ks][js].M1 = x1Flux[i].Mx; 
        pG->CGrid[ncg].myFlx[dim][ks][js].M2 = x1Flux[i].My;
        pG->CGrid[ncg].myFlx[dim][ks][js].M3 = x1Flux[i].Mz; 
#ifndef BAROTROPIC
        pG->CGrid[ncg].myFlx[dim][ks][js].E  = x1Flux[i].E; 
#endif /* BAROTROPIC */
#ifdef MHD
        pG->CGrid[ncg].myFlx[dim][ks][js].B1c = 0.0;
        pG->CGrid[ncg].myFlx[dim][ks][js].B2c = x1Flux[i].By; 
        pG->CGrid[ncg].myFlx[dim][ks][js].B3c = x1Flux[i].Bz; 
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++)
          pG->CGrid[ncg].myFlx[dim][ks][js].s[n]  = x1Flux[i].s[n]; 
#endif
      }
    }
  }

/* x1-boundaries of parent Grids (at boundaries of THIS Grid)  */
  for (npg=0; npg<pG->NPGrid; npg++) {
    for (dim=0; dim<2; dim++){
      if (pG->PGrid[npg].myFlx[dim] != NULL) {

        if (dim==0) i = pG->PGrid[npg].ijks[0];
        if (dim==1) i = pG->PGrid[npg].ijke[0] + 1;

        pG->PGrid[npg].myFlx[dim][ks][js].d  = x1Flux[i].d; 
        pG->PGrid[npg].myFlx[dim][ks][js].M1 = x1Flux[i].Mx; 
        pG->PGrid[npg].myFlx[dim][ks][js].M2 = x1Flux[i].My;
        pG->PGrid[npg].myFlx[dim][ks][js].M3 = x1Flux[i].Mz; 
#ifndef BAROTROPIC
        pG->PGrid[npg].myFlx[dim][ks][js].E  = x1Flux[i].E; 
#endif /* BAROTROPIC */
#ifdef MHD
        pG->PGrid[npg].myFlx[dim][ks][js].B1c = 0.0;
        pG->PGrid[npg].myFlx[dim][ks][js].B2c = x1Flux[i].By; 
        pG->PGrid[npg].myFlx[dim][ks][js].B3c = x1Flux[i].Bz; 
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++)
          pG->PGrid[npg].myFlx[dim][ks][js].s[n]  = x1Flux[i].s[n]; 
#endif
      }
    }
  }

#endif /* STATIC_MESH_REFINEMENT */







  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_init_1d(MeshS *pM)
 *  \brief Allocate temporary integration arrays 
*/
void integrate_init_1d(MeshS *pM)
{

  int size1=0,nl,nd;



/* Cycle over all Grids on this processor to find maximum Nx1 */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        if ((pM->Domain[nl][nd].Grid->Nx[0]) > size1){
          size1 = pM->Domain[nl][nd].Grid->Nx[0];
        }
      }
    }
  }


  size1 = size1 + 2*nghost;

#ifdef MHD
  /*if ((Bxc = (Real*)malloc(size1*sizeof(Real))) == NULL) goto on_error;*/
  if ((Bxc = (Real*)malloc(size1*sizeof(Real))) == NULL) goto on_error;
  if ((Bxb = (Real*)malloc(size1*sizeof(Real))) == NULL) goto on_error;
#endif
  if ((temp = (Real*)malloc(size1*sizeof(Real))) == NULL) goto on_error;
  if ((grad = (Real*)malloc(size1*sizeof(Real))) == NULL) goto on_error;

  //if ((Bxi = (Real*)malloc(size1*sizeof(Real))) == NULL) goto on_error;
  if ((Uinit   =(ConsS***)calloc_3d_array(1,1,size1,sizeof(ConsS)))
    == NULL) goto on_error;
  if ((U1d       =(Cons1DS*)malloc(size1*sizeof(Cons1DS)))==NULL) goto on_error;
  if ((Uc_x1 =(Cons1DS*)malloc(size1*sizeof(Cons1DS)))==NULL) goto on_error;
/*if ((Ur_x1Face =(ConsS*)malloc(size1*sizeof(ConsS)))==NULL) goto on_error;*/
  if ((x1Flux    =(Cons1DS*)malloc(size1*sizeof(Cons1DS)))==NULL) goto on_error;

  if ((W  = (Prim1DS*)malloc(size1*sizeof(Prim1DS))) == NULL) goto on_error;
  /*if ((Wl = (Prim1DS*)malloc(size1*sizeof(Prim1DS))) == NULL) goto on_error;
  if ((Wr = (Prim1DS*)malloc(size1*sizeof(Prim1DS))) == NULL) goto on_error;*/



#ifdef CYLINDRICAL
  if((StaticGravPot != NULL) || (CoolingFunc != NULL))
#endif
  {
    if ((dhalf  = (Real*)malloc(size1*sizeof(Real))) == NULL) goto on_error;
  }
  if(CoolingFunc != NULL){
    if ((phalf  = (Real*)malloc(size1*sizeof(Real))) == NULL) goto on_error;
  }

#ifdef CYLINDRICAL
  if ((geom_src = (Real*)calloc_1d_array(size1, sizeof(Real))) == NULL)
    goto on_error;
#endif


  return;

  on_error:

    integrate_destruct();
    ath_error("[integrate_init_1d]: malloc returned a NULL pointer\n");




}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_destruct_2d(void)
 *  \brief Free temporary integration arrays 
 */
void integrate_destruct_1d(void)
{

  
  /*if (Bxc != NULL) free(Bxc);*/

#ifdef MHD
  if (Bxc != NULL) free(Bxc);
  if (Bxb != NULL) free(Bxb);
#endif
  if (temp != NULL) free(temp);
  if (grad != NULL) free(grad);
  //if (Bxi != NULL) free(Bxi);
  if (Uinit    != NULL) free_3d_array(Uinit);
  if (U1d != NULL) free(U1d);
  if (Uc_x1 != NULL) free(Uc_x1);
  /*if (Ur_x1Face != NULL) free(Ur_x1Face);*/
  if (x1Flux != NULL) free(x1Flux);

  if (W  != NULL) free(W);
  /*if (Wl != NULL) free(Wl);
  if (Wr != NULL) free(Wr);*/

  if (dhalf != NULL) free(dhalf);
  if (phalf != NULL) free(phalf);

#ifdef CYLINDRICAL
  if (geom_src != NULL) free_1d_array((void*)geom_src);
#endif

  return;
  
  
  
  
}

/*=========================== PRIVATE FUNCTIONS ==============================*/



static void computemaxc(ConsS ***Uint, GridS *pG, int dim)
{

	/* Calculate cmax_idim=cfast_i+abs(v_idim) within ix^L
	! where cfast_i=sqrt(0.5*(cf**2+sqrt(cf**4-4*cs**2*b_i**2/rho)))
	! and cf**2=b**2/rho+cs**2/rho is the square of the speed of the fast wave
	! perpendicular to the magnetic field, and cs is the sound speed.*/

	int il,iu, is = pG->is, ie = pG->ie;
	int jl,ju, js = pG->js, je = pG->je;
	int kl,ku, ks = pG->ks, ke = pG->ke;

        int i1,i2,i3,n1z,n2z,n3z;
        int iss,jss,kss;

	register Real cmax=0;
	Real rhotot,rhototsq,pthermal,cs2,cfast2;
	Real cfasttemp,momfield,bfield;

        /*rho, mom1, mom2, mom3, energy, b1, b2, b3*/
	/* With particles, one more ghost cell must be updated in predict step */
	#ifdef PARTICLES
	  Real d1;
	  il = is - 3;
	  iu = ie + 3;
	  jl = js - 3;
	  ju = je + 3;
	  kl = ks - 3;
	  ku = ke + 3;
	#else
	  il = is - 2;
	  iu = ie + 2;
	  jl = js - 2;
	  ju = je + 2;
	  kl = ks - 2;
	  ku = ke + 2;
	#endif
        kl=0;
        ku=0;
        jl=0;
        ju=0;


	if (pG->Nx[0] > 1)
		n1z = pG->Nx[0] + 2*nghost;
	else
		n1z = 1;

	if (pG->Nx[1] > 1)
		n2z = pG->Nx[1] + 2*nghost;
	else
		n2z = 1;

	if (pG->Nx[2] > 1)
		n3z = pG->Nx[2] + 2*nghost;
	else
		n3z = 1;


        // TODO
        //getcmax();



	//for(i1=0;i1<n1z;i1++)
	//for(i2=0;i2<n2z;i2++)
	//for(i3=0;i1<n3z;i3++)
	//i3=ks;
        //i2=js;
	//for (i1=il; i1<=iu; i1++)
  for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {



		rhotot=(Uinit[i3][i2][i1].d+Uinit[i3][i2][i1].db);
		rhototsq=rhotot*rhotot;
		pthermal=Uinit[i3][i2][i1].E -((Uinit[i3][i2][i1].M1*Uinit[i3][i2][i1].M1+Uinit[i3][i2][i1].M2*Uinit[i3][i2][i1].M2+Uinit[i3][i2][i1].M3*Uinit[i3][i2][i1].M3)/rhotot)  ;


#ifdef MHD
		pthermal-=0.5*((Uinit[i3][i2][i1].B1c*Uinit[i3][i2][i1].B1c+Uinit[i3][i2][i1].B2c*Uinit[i3][i2][i1].B2c+Uinit[i3][i2][i1].B3c*Uinit[i3][i2][i1].B3c));
		pthermal-=0.5*((Uinit[i3][i2][i1].B1c*Uinit[i3][i2][i1].B1cb+Uinit[i3][i2][i1].B2c*Uinit[i3][i2][i1].B2cb+Uinit[i3][i2][i1].B3c*Uinit[i3][i2][i1].B3cb));
#endif
		pthermal*=(Gamma_1);
                //printf("cdens, cmax=%f %f %f %f\n",rhototsq, pthermal,cs2,Gamma_1);


#ifdef MHD
		cs2=Gamma*(pthermal+(Gamma_1)*(Uinit[i3][i2][i1].Eb-0.5*(    ((Uinit[i3][i2][i1].B1cb*Uinit[i3][i2][i1].B1cb+Uinit[i3][i2][i1].B2cb*Uinit[i3][i2][i1].B2cb+Uinit[i3][i2][i1].B3cb*Uinit[i3][i2][i1].B3cb)  ))));
#else
		cs2=Gamma*(pthermal+(Gamma_1)*(Uinit[i3][i2][i1].Eb));
#endif
		cs2/=rhototsq;




		pG->Hv[i3][i2][i1].csound=sqrt(cs2);
                cmax=MAX(cmax,pG->Hv[i3][i2][i1].csound);

#ifdef MHD
		cfast2=cs2+((Uinit[i3][i2][i1].B1c+Uinit[i3][i2][i1].B1cb)*(Uinit[i3][i2][i1].B1c+Uinit[i3][i2][i1].B1cb)+
                        (Uinit[i3][i2][i1].B2c+Uinit[i3][i2][i1].B2cb)*(Uinit[i3][i2][i1].B2c+Uinit[i3][i2][i1].B2cb)+
			(Uinit[i3][i2][i1].B3c+Uinit[i3][i2][i1].B3cb)*(Uinit[i3][i2][i1].B3c+Uinit[i3][i2][i1].B3cb))/(rhotot);




		pG->Hv[i3][i2][i1].cfast=sqrt(cfast2);
		cmax=MAX(cmax,pG->Hv[i3][i2][i1].cfast);



		switch(dim)
		{
			case 0:
				bfield=Uinit[i3][i2][i1].B1c+Uinit[i3][i2][i1].B1cb;
				momfield=Uinit[i3][i2][i1].M1;
			break;

			case 1:
				bfield=Uinit[i3][i2][i1].B2c+Uinit[i3][i2][i1].B2cb;
				momfield=Uinit[i3][i2][i1].M2;
			break;

			case 2:
				bfield=Uinit[i3][i2][i1].B3c+Uinit[i3][i2][i1].B3cb;
				momfield=Uinit[i3][i2][i1].M3;
			break;
		}
		cfasttemp=cfast2+sqrt(cfast2*cfast2-4*cs2*bfield*bfield/rhotot);
                pG->Hv[i3][i2][i1].cmaxd=sqrt(cfasttemp/2)+(momfield/rhotot);
		cmax=MAX(cmax,cfasttemp);
#endif

				}


        		}
		}

	pG->cmax=cmax;

}






static void hyperdifviscr(int fieldi,int dim,ConsS ***Uinit, GridS *pG)
{
	Real ***wtemp1=NULL, ***wtemp2=NULL, ***wtemp3=NULL, ***tmpnu=NULL, ***d1=NULL, ***d3=NULL, ***fieldd=NULL;
        Real maxt1,maxt2;

	int n1z,n2z,n3z;
        int i,j,k;

	int il,iu, is = pG->is, ie = pG->ie;
	int jl,ju, js = pG->js, je = pG->je;
	int kl,ku, ks = pG->ks, ke = pG->ke;

        int i1,i2,i3;
        int iss,jss,kss;

	Real dtodx1 = pG->dt/pG->dx1, dtodx2 = pG->dt/pG->dx2;

        /*rho, mom1, mom2, mom3, energy, b1, b2, b3*/
	/* With particles, one more ghost cell must be updated in predict step */
	#ifdef PARTICLES
	  Real d1;
	  il = is - 3;
	  iu = ie + 3;
	  jl = js - 3;
	  ju = je + 3;
	  kl = ks - 3;
	  ku = ke + 3;
	#else
	  il = is - 2;
	  iu = ie + 2;
	  jl = js - 2;
	  ju = je + 2;
	  kl = ks - 2;
	  ku = ke + 2;
	#endif

        kl=0;
        ku=0;

	if (pG->Nx[0] > 1)
		n1z = pG->Nx[0] + 2*nghost;
	else
		n1z = 1;

	if (pG->Nx[1] > 1)
		n2z = pG->Nx[1] + 2*nghost;
	else
		n2z = 1;

	if (pG->Nx[2] > 1)
		n3z = pG->Nx[2] + 2*nghost;
	else
		n3z = 1;

	wtemp1 = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	wtemp2 = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	wtemp3 = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	fieldd = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmpnu = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	d3 = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	d1 = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));


        printf("viscr maxc %d %d %d\n",n3z,n2z,n1z);
        // TODO
        //getcmax();
        computemaxc(Uinit, pG, dim);

	//for(i1=0;i1<n1z;i1++)
	//for(i2=0;i2<n2z;i2++)
	//for(i3=0;i1<n3z;i3++)
	//i3=ks;
        //i2=js;
	//for (i1=il; i1<=iu; i1++)
        printf("viscr after maxc\n");
  for (i3=kl; i3<ku; i3++) {
    for (i2=jl; i2<ju; i2++) {
    	for (i1=il; i1<iu; i1++) {

		switch(fieldi)
		{
		case rho:
			fieldd[i3][i2][i1]=Uinit[i3][i2][i1].d;
		break;
		case mom1:
			fieldd[i3][i2][i1]=Uinit[i3][i2][i1].M1;
		break;
		case mom2:
			fieldd[i3][i2][i1]=Uinit[i3][i2][i1].M2;
		break;
		case mom3:
			fieldd[i3][i2][i1]=Uinit[i3][i2][i1].M3;
		break;
		case energy:
			fieldd[i3][i2][i1]=Uinit[i3][i2][i1].E;
		break;
#ifdef MHD
		case b1:
			fieldd[i3][i2][i1]=Uinit[i3][i2][i1].B1c;
		break;
		case b2:
			fieldd[i3][i2][i1]=Uinit[i3][i2][i1].B2c;
		break;
		case b3:
			fieldd[i3][i2][i1]=Uinit[i3][i2][i1].B3c;
		break;
#endif
		}
        }
}
}

printf("fields define\n");


  for (i3=kl; i3<ku; i3++) {
    for (i2=jl; i2<ju; i2++) {
    	for (i1=il; i1<iu; i1++) {
		wtemp1[i3][i2][i1]=0.0;
		wtemp2[i3][i2][i1]=0.0;
		d3[i3][i2][i1]=0.0;

		pG->Hv[i3][i2][i1].hdnur[dim][fieldi]=0.0;

	       if(fieldi==energy)

#ifdef MHD
		wtemp1[i3][i2][i1]=fieldd[i3][i2][i1]-0.5*((Uinit[i3][i2][i1].B1c*Uinit[i3][i2][i1].B1c+Uinit[i3][i2][i1].B2c*Uinit[i3][i2][i1].B2c+Uinit[i3][i2][i1].B3c*Uinit[i3][i2][i1].B3c)
	+(Uinit[i3][i2][i1].M1*Uinit[i3][i2][i1].M1+Uinit[i3][i2][i1].M2*Uinit[i3][i2][i1].M2+Uinit[i3][i2][i1].M3*Uinit[i3][i2][i1].M3)/(Uinit[i3][i2][i1].d+Uinit[i3][i2][i1].db ));
#else
		wtemp1[i3][i2][i1]=fieldd[i3][i2][i1]-0.5*((Uinit[i3][i2][i1].M1*Uinit[i3][i2][i1].M1+Uinit[i3][i2][i1].M2*Uinit[i3][i2][i1].M2+Uinit[i3][i2][i1].M3*Uinit[i3][i2][i1].M3)/(Uinit[i3][i2][i1].d+Uinit[i3][i2][i1].db ));

#endif
	       else
	       {
		wtemp1[i3][i2][i1]=fieldd[i3][i2][i1];
		if((fieldi ==mom1 || fieldi == mom2 || fieldi == mom3))
			wtemp1[i3][i2][i1]=fieldd[i3][i2][i1]/(Uinit[i3][i2][i1].d+Uinit[i3][i2][i1].db);

		}

		//comment removed below to test mpi 29/10/2013
		tmpnu[i3][i2][i1]=wtemp1[i3][i2][i1];


	}
}
}

printf("temp1 tmpnu fields define\n");

        // TODO boundary terms

  for (i3=kl; i3<ku; i3++) {
    for (i2=jl; i2<ju; i2++) {
    	for (i1=il; i1<iu; i1++) {

		   d3[i3][i2][i1]=fabs(3.0*(tmpnu[i3+(dim==2)][i2+(dim==1)][i1+(dim==0)] - tmpnu[i3][i2][i1] ) - (tmpnu[i3+2*(dim==2)][i2+2*(dim==1)][i1+2*(dim==0)] - tmpnu[i3-(dim==2)][i2-(dim==1)][i1-(dim==0)]   ));
}
}
}


  for (i3=kl; i3<ku; i3++) {
    for (i2=jl; i2<ju; i2++) {
    	for (i1=il; i1<iu; i1++) {
		   d1[i3][i2][i1]=fabs((tmpnu[i3+(dim==2)][i2+(dim==1)][i1+(dim==0)] - tmpnu[i3][i2][i1] ));
}
}
}

   /*to here*/
        maxt2=0.0;
        maxt1=0.0;
  for (i3=kl; i3<ku; i3++) {
    for (i2=jl; i2<ju; i2++) {
    	for (i1=il; i1<iu; i1++) {

           for(kss=-(dim==2); kss<=(dim==2); kss++)
           for(jss=-(dim==1); jss<=(dim==1); jss++)
           for(iss=-(dim==0); iss<=(dim==0); iss++)
           {
{
{
                   if(d3[i3+1+kss][i2+1+jss][i1+1+iss]>maxt1)
                         maxt1=d3[i3+1+kss][i2+1+jss][i1+1+iss];
}
}
}

           wtemp2[i3][i2][i1]=maxt1;
           for(kss=-2*(dim==2); kss<=2*(dim==2); kss++)
           for(jss=-2*(dim==1); jss<=2*(dim==1); jss++)
           for(iss=-2*(dim==0); iss<=2*(dim==0); iss++)
{
{
{
                                     if(d1[i3+1+kss][i2+1+jss][i1+1+iss]>maxt2)
                                            maxt2=d1[i3+1+kss][i2+1+jss][i1+1+iss];

}
}
}


            wtemp3[i3][i2][i1]=maxt2;
        }
        }
}


  for (i3=kl; i3<ku; i3++) {
    for (i2=jl; i2<ju; i2++) {
    	for (i1=il; i1<iu; i1++) {
           if(wtemp3[i3][i2][i1]>0)
           {

	pG->Hv[i3][i2][i1].hdnur[dim][fieldi]=((dim==0)*(pG->dx1)+(dim==1)*(pG->dx2)+(dim==2)*(pG->dx3))*(pG->cmax)*(pG->chyp[fieldi])*wtemp2[i3][i2][i1]/wtemp3[i3][i2][i1];


         }
         else
           	pG->Hv[i3][i2][i1].hdnur[dim][fieldi]=0;


            if(pG->Hv[i3][i2][i1].hdnur[dim][fieldi]>(pG->maxviscoef))
                        pG->maxviscoef=pG->Hv[i3][i2][i1].hdnur[dim][fieldi];
         //}
         }
}
}

printf("free memory \n");


	if (wtemp1 != NULL) free(wtemp1);
	if (wtemp2 != NULL) free(wtemp2);
	if (wtemp3 != NULL) free(wtemp3);
	if (fieldd != NULL) free(fieldd);
	if (tmpnu != NULL) free(tmpnu);
	if (d3 != NULL) free(d3);
	if (d1 != NULL) free(d1);

printf("memory freed\n");
	return;

}

static void hyperdifviscl(int fieldi,int dim,ConsS ***Uinit, GridS *pG)
{
	Real ***wtemp1=NULL, ***wtemp2=NULL, ***wtemp3=NULL, ***tmpnu=NULL, ***d1=NULL, ***d3=NULL, ***fieldd=NULL;
        Real maxt1,maxt2;

	int n1z,n2z,n3z;
        int i,j,k;

	int il,iu, is = pG->is, ie = pG->ie;
	int jl,ju, js = pG->js, je = pG->je;
	int kl,ku, ks = pG->ks, ke = pG->ke;


        int i1,i2,i3;
        int iss,jss,kss;

        /*rho, mom1, mom2, mom3, energy, b1, b2, b3*/

	/* With particles, one more ghost cell must be updated in predict step */
	#ifdef PARTICLES
	  Real d1;
	  il = is - 3;
	  iu = ie + 3;
	  jl = js - 3;
	  ju = je + 3;
	  kl = ks - 3;
	  ku = ke + 3;
	#else
	  il = is - 2;
	  iu = ie + 2;
	  jl = js - 2;
	  ju = je + 2;
	  kl = ks - 2;
	  ku = ke + 2;
	#endif

        kl=0;
        ku=0;

	if (pG->Nx[0] > 1)
		n1z = pG->Nx[0] + 2*nghost;
	else
		n1z = 1;

	if (pG->Nx[1] > 1)
		n2z = pG->Nx[1] + 2*nghost;
	else
		n2z = 1;

	if (pG->Nx[2] > 1)
		n3z = pG->Nx[2] + 2*nghost;
	else
		n3z = 1;

	wtemp1 = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	wtemp2 = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	wtemp3 = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	fieldd = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmpnu = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	d3 = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	d1 = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));

        // TODO
        //getcmax();
	computemaxc(Uinit, pG, dim);

	//for(i1=0;i1<n1z;i1++)
	//for(i2=0;i2<n2z;i2++)
	//for(i3=0;i1<n3z;i3++)
	//i3=ks;
        //i2=js;
	//for (i1=il; i1<=iu; i1++)
  for (i3=kl; i3<ku; i3++) {
    for (i2=jl; i2<ju; i2++) {
    	for (i1=il; i1<iu; i1++) {

		switch(fieldi)
		{
		case rho:
			fieldd[i3][i2][i1]=Uinit[i3][i2][i1].d;
		break;
		case mom1:
			fieldd[i3][i2][i1]=Uinit[i3][i2][i1].M1;
		break;
		case mom2:
			fieldd[i3][i2][i1]=Uinit[i3][i2][i1].M2;
		break;
		case mom3:
			fieldd[i3][i2][i1]=Uinit[i3][i2][i1].M3;
		break;
		case energy:
			fieldd[i3][i2][i1]=Uinit[i3][i2][i1].E;
		break;
#ifdef MHD
		case b1:
			fieldd[i3][i2][i1]=Uinit[i3][i2][i1].B1c;
		break;
		case b2:
			fieldd[i3][i2][i1]=Uinit[i3][i2][i1].B2c;
		break;
		case b3:
			fieldd[i3][i2][i1]=Uinit[i3][i2][i1].B3c;
		break;
#endif
		}
        }
}
}




  for (i3=kl; i3<ku; i3++) {
    for (i2=jl; i2<ju; i2++) {
    	for (i1=il; i1<iu; i1++) {
		wtemp1[i3][i2][i1]=0.0;
		wtemp2[i3][i2][i1]=0.0;
		d3[i3][i2][i1]=0.0;

		pG->Hv[i3][i2][i1].hdnul[dim][fieldi]=0.0;

	       if(fieldi==energy)
#ifdef MHD
		wtemp1[i3][i2][i1]=fieldd[i3][i2][i1]-0.5*((Uinit[i3][i2][i1].B1c*Uinit[i3][i2][i1].B1c+Uinit[i3][i2][i1].B2c*Uinit[i3][i2][i1].B2c+Uinit[i3][i2][i1].B3c*Uinit[i3][i2][i1].B3c)
	+(Uinit[i3][i2][i1].M1*Uinit[i3][i2][i1].M1+Uinit[i3][i2][i1].M2*Uinit[i3][i2][i1].M2+Uinit[i3][i2][i1].M3*Uinit[i3][i2][i1].M3)/(Uinit[i3][i2][i1].d+Uinit[i3][i2][i1].db ));
#else
		wtemp1[i3][i2][i1]=fieldd[i3][i2][i1]-0.5*((Uinit[i3][i2][i1].M1*Uinit[i3][i2][i1].M1+Uinit[i3][i2][i1].M2*Uinit[i3][i2][i1].M2+Uinit[i3][i2][i1].M3*Uinit[i3][i2][i1].M3)/(Uinit[i3][i2][i1].d+Uinit[i3][i2][i1].db ));


#endif
	       else
	       {
		wtemp1[i3][i2][i1]=fieldd[i3][i2][i1];
		if((fieldi ==mom1 || fieldi == mom2 || fieldi == mom3))
			wtemp1[i3][i2][i1]=fieldd[i3][i2][i1]/(Uinit[i3][i2][i1].d+Uinit[i3][i2][i1].db);

		}

		//comment removed below to test mpi 29/10/2013
		tmpnu[i3][i2][i1]=wtemp1[i3][i2][i1];


	}
}
}




        // TODO boundary terms

  for (i3=kl; i3<ku; i3++) {
    for (i2=jl; i2<ju; i2++) {
    	for (i1=il; i1<iu; i1++) {

		   d3[i3][i2][i1]=fabs(3.0*( tmpnu[i3][i2][i1] -  tmpnu[i3+(dim==2)][i2+(dim==1)][i1+(dim==0)]  ) - ( tmpnu[i3-(dim==2)][i2-(dim==1)][i1-(dim==0)] - tmpnu[i3+2*(dim==2)][i2+2*(dim==1)][i1+2*(dim==0)]    ));
}
}
}





  for (i3=kl; i3<ku; i3++) {
    for (i2=jl; i2<ju; i2++) {
    	for (i1=il; i1<iu; i1++) {
		   d1[i3][i2][i1]=fabs(( tmpnu[i3][i2][i1] -   tmpnu[i3+(dim==2)][i2+(dim==1)][i1+(dim==0)]  ));
}
}
}

   /*to here*/
        maxt2=0.0;
        maxt1=0.0;
  for (i3=kl; i3<ku; i3++) {
    for (i2=jl; i2<ju; i2++) {
    	for (i1=il; i1<iu; i1++) {

           for(kss=-(dim==2); kss<=(dim==2); kss++)
           for(jss=-(dim==1); jss<=(dim==1); jss++)
           for(iss=-(dim==0); iss<=(dim==0); iss++)
           {
{
{
                   if(d3[i3+1+kss][i2+1+jss][i1+1+iss]>maxt1)
                         maxt1=d3[i3+1+kss][i2+1+jss][i1+1+iss];
}
}
}

           wtemp2[i3][i2][i1]=maxt1;
           for(kss=-2*(dim==2); kss<=2*(dim==2); kss++)
           for(jss=-2*(dim==1); jss<=2*(dim==1); jss++)
           for(iss=-2*(dim==0); iss<=2*(dim==0); iss++)
{
{
{
                                     if(d1[i3+1+kss][i2+1+jss][i1+1+iss]>maxt2)
                                            maxt2=d1[i3+1+kss][i2+1+jss][i1+1+iss];

}
}
}


            wtemp3[i3][i2][i1]=maxt2;
        }
        }
}


  for (i3=kl; i3<ku; i3++) {
    for (i2=jl; i2<ju; i2++) {
    	for (i1=il; i1<iu; i1++) {
           if(wtemp3[i3][i2][i1]>0)
           {

	pG->Hv[i3][i2][i1].hdnul[dim][fieldi]=((dim==0)*(pG->dx1)+(dim==1)*(pG->dx2)+(dim==2)*(pG->dx3))*(pG->cmax)*(pG->chyp[fieldi])*wtemp2[i3][i2][i1]/wtemp3[i3][i2][i1];


         }
         else
           	pG->Hv[i3][i2][i1].hdnul[dim][fieldi]=0;


            if(pG->Hv[i3][i2][i1].hdnul[dim][fieldi]>(pG->maxviscoef))
                        pG->maxviscoef=pG->Hv[i3][i2][i1].hdnul[dim][fieldi];
         //}
         }
}
}




	if (wtemp1 != NULL) free(wtemp1);
	if (wtemp2 != NULL) free(wtemp2);
	if (wtemp3 != NULL) free(wtemp3);
	if (fieldd != NULL) free(fieldd);

	if (tmpnu != NULL) free(tmpnu);
	if (d3 != NULL) free(d3);
	if (d1 != NULL) free(d1);


	return;


}





static void hyperdifrhosource(int dim,Real dt,ConsS ***Uinit, GridS *pG)
{
	Real ***wtempr=NULL, ***wtempl=NULL, ***wtemp3=NULL, ***tmp=NULL, ***tmp2=NULL, ***fieldd=NULL;
    Real maxt1,maxt2;
	Real nur,nul;

	int n1z,n2z,n3z;
        int i,j,k;

	int il,iu, is = pG->is, ie = pG->ie;
	int jl,ju, js = pG->js, je = pG->je;
	int kl,ku, ks = pG->ks, ke = pG->ke;


        int i1,i2,i3;
        int iss,jss,kss;

	int fieldi=rho;
        Real dtodx1 = pG->dt/pG->dx1, dtodx2 = pG->dt/pG->dx2, dtodx3 = pG->dt/pG->dx3;
        /*rho, mom1, mom2, mom3, energy, b1, b2, b3*/

	/* With particles, one more ghost cell must be updated in predict step */
	#ifdef PARTICLES
	  Real d1;
	  il = is - 3;
	  iu = ie + 3;
	  jl = js - 3;
	  ju = je + 3;
	  kl = ks - 3;
	  ku = ke + 3;
	#else
	  il = is - 2;
	  iu = ie + 2;
	  jl = js - 2;
	  ju = je + 2;
	  kl = ks - 2;
	  ku = ke + 2;
	#endif

        kl=0;
        ku=0;

	if (pG->Nx[0] > 1)
		n1z = pG->Nx[0] + 2*nghost;
	else
		n1z = 1;

	if (pG->Nx[1] > 1)
		n2z = pG->Nx[1] + 2*nghost;
	else
		n2z = 1;

	if (pG->Nx[2] > 1)
		n3z = pG->Nx[2] + 2*nghost;
	else
		n3z = 1;
switch(dim)
{
case 1:
	fieldd = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	wtempr = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	wtempl = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	wtemp3 = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmp = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmp2 = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
break;
case 2:
	fieldd = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	wtempr = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	wtempl = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	wtemp3 = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmp = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmp2 = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
break;
case 3:
	fieldd = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	wtempr = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	wtempl = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	wtemp3 = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmp = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmp2 = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
break;

}



     //CALL setnu(w,rho_,idim,ixOmin1,ixOmin2,ixOmax1,ixOmax2,nuR,nuL)

     //    tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)=w(ixImin1:ixImax1,&
     //   ixImin2:ixImax2,ixImin3:ixImax3,rho_)




 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			fieldd[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=Uinit[i3][i2][i1].d;
					}
				}
			}

//need gradient for different dimensions
     //CALL gradient1L(tmp,ixmin1,ixmin2,ixmax1,ixmax2,idim,tmp2)

switch(dim)
{

case 1:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1l(fieldd[i3][i2], n1z,pG->dx1,tmp2[i3][i2]);
				}
			}
break;

case 2:
 for (i3=kl; i3<=ku; i3++) {
    for (i1=il; i1<=iu; i1++) {
			gradient1l(fieldd[i3][i1], n2z,pG->dx2,tmp2[i3][i1]);
				}
			}
break;

case 3:
 for (i1=il; i1<=iu; i1++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1l(fieldd[i1][i2], n3z,pG->dx3,tmp2[i1][i2]);
				}
			}
break;

}



/*nur=pG->Hv[i3][i2][i1].hdnur[dim][fieldi];
nul=pG->Hv[i3][i2][i1].hdnur[dim][fieldi];*/

     /*tmpL(ixImin1:ixImax1,ixImin2:ixImax2)=(nuL(ixImin1:ixImax1,&
        ixImin2:ixImax2)+nushk(ixImin1:ixImax1,ixImin2:ixImax2,idim))&
        *tmp2(ixImin1:ixImax1,ixImin2:ixImax2) */

for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			wtempl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=(pG->Hv[i3][i2][i1].hdnul[dim][fieldi])*tmp2[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)];
					}
				}
			}

        //need gradient for different dimensions
     //CALL gradient1R(tmp,ixmin1,ixmin2,ixmax1,ixmax2,idim,tmp2)

switch(dim)
{

case 1:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1r(fieldd[i3][i2], n1z,pG->dx1,tmp2[i3][i2]);
				}
			}
break;

case 2:
 for (i3=kl; i3<=ku; i3++) {
    for (i1=il; i1<=iu; i1++) {
			gradient1r(fieldd[i3][i1], n2z,pG->dx2,tmp2[i3][i1]);
				}
			}
break;

case 3:
 for (i1=il; i1<=iu; i1++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1r(fieldd[i1][i2], n3z,pG->dx3,tmp2[i1][i2]);
				}
			}
break;

}







     /*tmpR(ixImin1:ixImax1,ixImin2:ixImax2)=(nuR(ixImin1:ixImax1,&
        ixImin2:ixImax2)+nushk(ixImin1:ixImax1,ixImin2:ixImax2,idim))&
        *tmp2(ixImin1:ixImax1,ixImin2:ixImax2)*/
for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			wtempr[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]  =(pG->Hv[i3][i2][i1].hdnur[dim][fieldi])*tmp2[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)];
					}
				}
			}



     /*wnew(ixImin1:ixImax1,ixImin2:ixImax2,rho_)=wnew(ixImin1:ixImax1,&
        ixImin2:ixImax2,rho_)+(tmpR(ixImin1:ixImax1,ixImin2:ixImax2)&
        -tmpL(ixImin1:ixImax1,ixImin2:ixImax2))/dx(ixImin1:ixImax1,&
        ixImin2:ixImax2,idim)*qdt*/
for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
     pG->U[i3][i2][i1].d  +=  dt*(wtempr[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]-wtempl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]);
                   }
	}
	}


	if (wtempr != NULL) free(wtempr);
	if (wtempl != NULL) free(wtempl);
	if (wtemp3 != NULL) free(wtemp3);
	if (tmp != NULL) free(tmp);
	if (tmp2 != NULL) free(tmp2);
	if (fieldd != NULL) free(fieldd);



	return;
}


static void hyperdifesource(int dim,Real dt,ConsS ***Uint, GridS *pG)
{
	Real ***wtempr=NULL, ***wtempl=NULL, ***wtemp3=NULL, ***tmp=NULL, ***tmp2=NULL, ***fieldd=NULL;
        Real maxt1,maxt2;
	Real nur,nul;

	int n1z,n2z,n3z;
        int i,j,k;

	int il,iu;
	int jl,ju;
	int kl,ku;

	int is,ie,js,je,ks,ke;

        int i1,i2,i3;
        int iss,jss,kss;

	int fieldi=energy;
        Real dtodx1 = pG->dt/pG->dx1, dtodx2 = pG->dt/pG->dx2, dtodx3 = pG->dt/pG->dx3;

	is = pG->is,
	ie = pG->ie;
	js = pG->js,
	je = pG->je;
	ks = pG->ks,
	ke = pG->ke;


        /*rho, mom1, mom2, mom3, energy, b1, b2, b3*/

	/* With particles, one more ghost cell must be updated in predict step */
	#ifdef PARTICLES
	  Real d1;
	  il = is - 3;
	  iu = ie + 3;
	  jl = js - 3;
	  ju = je + 3;
	  kl = ks - 3;
	  ku = ke + 3;
	#else
	  il = is - 2;
	  iu = ie + 2;
	  jl = js - 2;
	  ju = je + 2;
	  kl = ks - 2;
	  ku = ke + 2;
	#endif

        kl=0;
        ku=0;

	if (pG->Nx[0] > 1)
		n1z = pG->Nx[0] + 2*nghost;
	else
		n1z = 1;

	if (pG->Nx[1] > 1)
		n2z = pG->Nx[1] + 2*nghost;
	else
		n2z = 1;

	if (pG->Nx[2] > 1)
		n3z = pG->Nx[2] + 2*nghost;
	else
		n3z = 1;


switch(dim)
{
case 1:
	fieldd = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	wtempr = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	wtempl = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	wtemp3 = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmp = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmp2 = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
break;
case 2:
	fieldd = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	wtempr = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	wtempl = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	wtemp3 = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmp = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmp2 = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
break;
case 3:
	fieldd = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	wtempr = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	wtempl = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	wtemp3 = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmp = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmp2 = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
break;

}





     //CALL setnu(w,rho_,idim,ixOmin1,ixOmin2,ixOmax1,ixOmax2,nuR,nuL)


 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			fieldd[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=Uinit[i3][i2][i1].E;
					}
				}
			}



    // tmp(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
   //     e_)-half*((w(ixImin1:ixImax1,ixImin2:ixImax2,b1_)**2&
   //     +w(ixImin1:ixImax1,ixImin2:ixImax2,b2_)**2)+(w(ixImin1:ixImax1,&
   //     ixImin2:ixImax2,m1_)**2+w(ixImin1:ixImax1,ixImin2:ixImax2,m2_)**2)&
   //     /(w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)+w(ixImin1:ixImax1,&
   //     ixImin2:ixImax2,rhob_)))
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {


#ifdef MHD
			tmp[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]= fieldd[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)] -((Uinit[i3][i2][i1].B1c*Uinit[i3][i2][i1].B1c+Uinit[i3][i2][i1].B2c*Uinit[i3][i2][i1].B2c+Uinit[i3][i2][i1].B3c*Uinit[i3][i2][i1].B3c)/2)+((Uinit[i3][i2][i1].M1*Uinit[i3][i2][i1].M1+  Uinit[i3][i2][i1].M2*Uinit[i3][i2][i1].M2 +  Uinit[i3][i2][i1].M3*Uinit[i3][i2][i1].M3 ))/(Uinit[i3][i2][i1].d+Uinit[i3][i2][i1].db) ;
#else
			tmp[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]= fieldd[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)] +((Uinit[i3][i2][i1].M1*Uinit[i3][i2][i1].M1+  Uinit[i3][i2][i1].M2*Uinit[i3][i2][i1].M2 +  Uinit[i3][i2][i1].M3*Uinit[i3][i2][i1].M3 ))/(Uinit[i3][i2][i1].d+Uinit[i3][i2][i1].db) ;

#endif




					}
				}
			}


     //   CALL gradient1L(tmp,ixmin1,ixmin2,ixmax1,ixmax2,idim,tmp2)
//need gradient for different dimensions
switch(dim)
{

case 1:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1l(tmp[i3][i2], n1z,pG->dx1,tmp2[i3][i2]);
				}
			}
break;

case 2:
 for (i3=kl; i3<=ku; i3++) {
    for (i1=il; i1<=iu; i1++) {
			gradient1l(tmp[i3][i1], n2z,pG->dx2,tmp2[i3][i1]);
				}
			}
break;

case 3:
 for (i1=il; i1<=iu; i1++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1l(tmp[i1][i2], n3z,pG->dx3,tmp2[i1][i2]);
				}
			}
break;

}






  //   tmpL(ixImin1:ixImax1,ixImin2:ixImax2)=(nuL(ixImin1:ixImax1,&
  //      ixImin2:ixImax2)+nushk(ixImin1:ixImax1,ixImin2:ixImax2,idim))&
  //      *tmp2(ixImin1:ixImax1,ixImin2:ixImax2)
for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			wtempl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=(pG->Hv[i3][i2][i1].hdnul[dim][fieldi])*tmp2[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)];
					}
				}
			}



    //CALL gradient1R(tmp,ixmin1,ixmin2,ixmax1,ixmax2,idim,tmp2)
//need gradient for different dimensions

switch(dim)
{

case 1:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1r(tmp[i3][i2], n1z,pG->dx1,tmp2[i3][i2]);
				}
			}
break;

case 2:
 for (i3=kl; i3<=ku; i3++) {
    for (i1=il; i1<=iu; i1++) {
			gradient1r(tmp[i3][i1], n2z,pG->dx2,tmp2[i3][i1]);
				}
			}
break;

case 3:
 for (i1=il; i1<=iu; i1++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1r(tmp[i1][i2], n3z,pG->dx3,tmp2[i1][i2]);
				}
			}
break;

}






    // tmpR(ixImin1:ixImax1,ixImin2:ixImax2)=(nuR(ixImin1:ixImax1,&
    //    ixImin2:ixImax2)+nushk(ixImin1:ixImax1,ixImin2:ixImax2,idim))&
    //    *tmp2(ixImin1:ixImax1,ixImin2:ixImax2)
for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			wtempr[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=(pG->Hv[i3][i2][i1].hdnur[dim][fieldi])*tmp2[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)];
					}
				}
			}








/*nur=pG->Hv[i3][i2][i1].hdnur[dim][fieldi];
nul=pG->Hv[i3][i2][i1].hdnur[dim][fieldi];*/

    // wnew(ixImin1:ixImax1,ixImin2:ixImax2,e_)=wnew(ixImin1:ixImax1,&
    //    ixImin2:ixImax2,e_)+(tmpR(ixImin1:ixImax1,ixImin2:ixImax2)&
    //    -tmpL(ixImin1:ixImax1,ixImin2:ixImax2))/dx(ixImin1:ixImax1,&
    //    ixImin2:ixImax2,idim)*qdt
for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
     pG->U[k][j][i].E  += (dtodx1*(dim==1)+dtodx2*(dim==2)+dtodx3*(dim==3))*(wtempr[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]-wtempl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]);
	}
	}
	}


	if (wtempr != NULL) free(wtempr);
	if (wtempl != NULL) free(wtempl);
	if (wtemp3 != NULL) free(wtemp3);
	if (tmp != NULL) free(tmp);
	if (tmp2 != NULL) free(tmp2);
	if (fieldd != NULL) free(fieldd);



	return;

}


static void hyperdifmomsource(int kf,int nf,int ii,int ii0,Real dt,ConsS ***Uint, GridS *pG)
{

//ii maps to the dimension - (k below in sac code)

// ii0 maps to the field id

int dim=kf;
int fieldi=nf;
int lf=nf;


	//Real ***tmprhor=NULL, ***tmprhol=NULL, ***tmpvr=NULL, ***tmpvl=NULL, ***tmpr=NULL, ***tmpl=NULL, ***tmp=NULL, ***tmp2=NULL, ***fieldd=NULL;
        Real maxt1,maxt2;
	Real nur,nul;

	int n1z,n2z,n3z;
        int i,j,k;

	int il,iu;
	int jl,ju;
	int kl,ku;

	int is,ie,js,je,ks,ke;

        int i1,i2,i3;
        int iss,jss,kss;

	//fieldi=energy;
        Real dtodx1 = pG->dt/pG->dx1, dtodx2 = pG->dt/pG->dx2, dtodx3 = pG->dt/pG->dx3;

	is = pG->is,
	ie = pG->ie;
	js = pG->js,
	je = pG->je;
	ks = pG->ks,
	ke = pG->ke;
        Real dx;


        /*rho, mom1, mom2, mom3, energy, b1, b2, b3*/

	/* With particles, one more ghost cell must be updated in predict step */
	#ifdef PARTICLES
	  Real d1;
	  il = is - 3;
	  iu = ie + 3;
	  jl = js - 3;
	  ju = je + 3;
	  kl = ks - 3;
	  ku = ke + 3;
	#else
	  il = is - 2;
	  iu = ie + 2;
	  jl = js - 2;
	  ju = je + 2;
	  kl = ks - 2;
	  ku = ke + 2;
	#endif

        kl=0;
        ku=0;

        dx= (pG->dx1)*(dim==1)+(pG->dx2)*(dim==2)+(pG->dx3)*(dim==3);


	if (pG->Nx[0] > 1)
		n1z = pG->Nx[0] + 2*nghost;
	else
		n1z = 1;

	if (pG->Nx[1] > 1)
		n2z = pG->Nx[1] + 2*nghost;
	else
		n2z = 1;

	if (pG->Nx[2] > 1)
		n3z = pG->Nx[2] + 2*nghost;
	else
		n3z = 1;


	Real ***tmprhor=NULL, ***tmprhol=NULL, ***tmpvr=NULL, ***tmpvl=NULL, ***tmpr=NULL, ***tmpl=NULL, ***tmp=NULL, ***tmp2=NULL, ***fieldf=NULL,***fieldl=NULL;


switch(dim)
{
case 1:
	fieldl = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	fieldf = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmpr = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmpl = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmprhor = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmprhol = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmpvr = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmpvl = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmp = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmp2 = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
break;
case 2:
	fieldf = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	fieldl = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmpr = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmpl = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmprhor = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmprhol = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmpvr = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmpvl = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmp = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmp2 = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
break;
case 3:
	fieldf = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	fieldl = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmpr = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmpl = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmprhor = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmprhol = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmpvr = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmpvl = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmp = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmp2 = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
break;

}





switch(fieldi)
{

case 1:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			fieldf[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=Uinit[i3][i2][i1].M1;
					}
				}
			}
break;

case 2:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			fieldf[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=Uinit[i3][i2][i1].M2;
					}
				}
			}
break;

case 3:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			fieldf[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=Uinit[i3][i2][i1].M3;
					}
				}
			}
break;


}



switch(ii0)
{

case 1:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			fieldl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=Uinit[i3][i2][i1].M1;
					}
				}
			}
break;

case 2:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			fieldl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=Uinit[i3][i2][i1].M2;
					}
				}
			}
break;

case 3:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			fieldl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=Uinit[i3][i2][i1].M3;
					}
				}
			}
break;


}

//for(dim=0; dim<2; dim++) //each direction   //k
//for( fieldi=0; fieldi<2; fieldi++)          //l

    /*j is + h is -
    /* tmprhoL(ixmin1:ixmax1,ixmin2:ixmax2)=((w(ixmin1:ixmax1,ixmin2:ixmax2,&
        rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_))+(w(hxmin1:hxmax1,&
        hxmin2:hxmax2,rho_)+w(hxmin1:hxmax1,hxmin2:hxmax2,rhob_)))/two  */
 for (i3=(kl+(dim==3)); i3<=ku; i3++) {
    for (i2=(jl+(dim==2)); i2<=ju; i2++) {
    	for (i1=(il+(dim==1)); i1<=iu; i1++) {

			tmprhol[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=(Uinit[i3][i2][i1].d+Uinit[i3][i2][i1].db+Uinit[i3-(dim==3)][i2-(dim==2)][i1-(dim==1)].d+Uinit[i3-(dim==3)][i2-(dim==2)][i1-(dim==1)].db);
					}
				}
			}

   /*  tmprhoR(ixmin1:ixmax1,ixmin2:ixmax2)=((w(jxmin1:jxmax1,jxmin2:jxmax2,&
        rho_)+w(jxmin1:jxmax1,jxmin2:jxmax2,rhob_))+(w(ixmin1:ixmax1,&
        ixmin2:ixmax2,rho_)+w(ixmin1:ixmax1,ixmin2:ixmax2,rhob_)))/two  */
 for (i3=(kl); i3<=(ku-(dim==3)); i3++) {
    for (i2=(jl); i2<=(ju-(dim==2)); i2++) {
    	for (i1=(il); i1<=(iu-(dim==1)); i1++) {

			tmprhor[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=(Uinit[i3][i2][i1].d+Uinit[i3][i2][i1].db+Uinit[i3+(dim==3)][i2+(dim==2)][i1+(dim==1)].d+Uinit[i3+(dim==3)][i2+(dim==2)][i1+(dim==1)].db);
					}
				}
			}


       /* tmp(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,m0_+l)/(w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)&
           +w(ixImin1:ixImax1,ixImin2:ixImax2,rhob_))*/
 for (i3=(kl); i3<=(ku-(dim==3)); i3++) {
    for (i2=(jl); i2<=(ju-(dim==2)); i2++) {
    	for (i1=(il); i1<=(iu-(dim==1)); i1++) {

			tmp[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=fieldf[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]/(Uinit[i3][i2][i1].d+Uinit[i3][i2][i1].db);
					}
				}
			}


            /*  tmpVL(ixmin1:ixmax1,ixmin2:ixmax2)=(w(ixmin1:ixmax1,&
                 ixmin2:ixmax2,m0_+ii0)+w(hxmin1:hxmax1,hxmin2:hxmax2,m0_&
                 +ii0))/two   */
 for (i3=(kl+(dim==3)); i3<=(ku); i3++) {
    for (i2=(jl+(dim==2)); i2<=(ju); i2++) {
    	for (i1=(il+(dim==1)); i1<=(iu); i1++) {

			tmpvl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=fieldl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]+fieldl[AIN3(i1-(dim==1),i2-(dim==2),i3-(dim==3),dim)][AIN2((dim==1),i2-(dim==2),i3-(dim==3),dim)][AIN1((dim==1),i2-(dim==2),i3-(dim==3),dim)]/2;
					}
				}
			}

             /* tmpVR(ixmin1:ixmax1,ixmin2:ixmax2)=(w(jxmin1:jxmax1,&
                 jxmin2:jxmax2,m0_+ii0)+w(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
                 +ii0))/two */
 for (i3=(kl); i3<=(ku-(dim==3)); i3++) {
    for (i2=(jl); i2<=(ju-(dim==2)); i2++) {
    	for (i1=(il); i1<=(iu-(dim==1)); i1++) {

			tmpvr[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=fieldl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]+fieldl[AIN3(i1+(dim==1),i2+(dim==2),i3+(dim==3),dim)][AIN2((dim==1),i2+(dim==2),i3+(dim==3),dim)][AIN1((dim==1),i2+(dim==2),i3+(dim==3),dim)]/2;
					}
				}
			}


           /*   CALL gradient1L(tmp,ixmin1,ixmin2,ixmax1,ixmax2,k,tmp2)*/
switch(dim)
{

case 1:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1l(tmp[i3][i2], n1z,pG->dx1,tmp2[i3][i2]);
				}
			}
break;

case 2:
 for (i3=kl; i3<=ku; i3++) {
    for (i1=il; i1<=iu; i1++) {
			gradient1l(tmp[i3][i1], n2z,pG->dx2,tmp2[i3][i1]);
				}
			}
break;

case 3:
 for (i1=il; i1<=iu; i1++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1l(tmp[i1][i2], n3z,pG->dx3,tmp2[i1][i2]);
				}
			}
break;

}




           /*   tmpL(ixImin1:ixImax1,ixImin2:ixImax2)=(nuL(ixImin1:ixImax1,&
                 ixImin2:ixImax2)+nushk(ixImin1:ixImax1,ixImin2:ixImax2,k))&
                 *tmp2(ixImin1:ixImax1,ixImin2:ixImax2)*/
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			tmpl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=tmp2[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]*(pG->Hv[i3][i2][i1].hdnul[dim][fieldi]);

					}
				}
			}




           /*   CALL gradient1R(tmp,ixmin1,ixmin2,ixmax1,ixmax2,k,tmp2)*/

switch(dim)
{

case 1:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1r(tmp[i3][i2], n1z,pG->dx1,tmp2[i3][i2]);
				}
			}
break;

case 2:
 for (i3=kl; i3<=ku; i3++) {
    for (i1=il; i1<=iu; i1++) {
			gradient1r(tmp[i3][i1], n2z,pG->dx2,tmp2[i3][i1]);
				}
			}
break;

case 3:
 for (i1=il; i1<=iu; i1++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1r(tmp[i1][i2], n3z,pG->dx3,tmp2[i1][i2]);
				}
			}
break;

}


           /*   tmpR(ixImin1:ixImax1,ixImin2:ixImax2)=(nuR(ixImin1:ixImax1,&
                 ixImin2:ixImax2)+nushk(ixImin1:ixImax1,ixImin2:ixImax2,k))&
                 *tmp2(ixImin1:ixImax1,ixImin2:ixImax2) */
for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			tmpr[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=tmp2[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]*(pG->Hv[i3][i2][i1].hdnur[dim][fieldi]);

					}
				}
			}


            /*  tmp2(ixImin1:ixImax1,ixImin2:ixImax2)=(tmprhoR(ixImin1:ixImax1,&
                 ixImin2:ixImax2)*tmpR(ixImin1:ixImax1,ixImin2:ixImax2)&
                 -tmprhoL(ixImin1:ixImax1,ixImin2:ixImax2)*tmpL&
                 (ixImin1:ixImax1,ixImin2:ixImax2))/dx(ixImin1:ixImax1,&
                 ixImin2:ixImax2,k)/two */
for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			tmp2[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=(tmprhor[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]*tmpr[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]-tmprhol[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]*tmpl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)])/2;

					}
				}
			}

            /*  wnew(ixImin1:ixImax1,ixImin2:ixImax2,m0_+ii0)&
                 =wnew(ixImin1:ixImax1,ixImin2:ixImax2,m0_+ii0)&
                 +tmp2(ixImin1:ixImax1,ixImin2:ixImax2)*qdt */



for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {

switch(ii0)
{
case 1:
     pG->U[i3][i2][i1].M1  += dt*(tmp2[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]);
break;

case 2:
     pG->U[i3][i2][i1].M2  += dt*(tmp2[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]);
break;

case 3:
     pG->U[i3][i2][i1].M3  += dt*(tmp2[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]);
break;

}


	}
	}
	}






           /*   tmp2(ixImin1:ixImax1,ixImin2:ixImax2)=(tmpVR(ixImin1:ixImax1,&
                 ixImin2:ixImax2)*tmpR(ixImin1:ixImax1,ixImin2:ixImax2)&
                 -tmpVL(ixImin1:ixImax1,ixImin2:ixImax2)*tmpL(ixImin1:ixImax1,&
                 ixImin2:ixImax2))/dx(ixImin1:ixImax1,ixImin2:ixImax2,k)/two   */
for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			tmp2[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=(tmpvr[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]*tmpr[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]-tmpvl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]*tmpl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)])/2;

					}
				}
			}



          /*    wnew(ixImin1:ixImax1,ixImin2:ixImax2,e_)=wnew(ixImin1:ixImax1,&
                 ixImin2:ixImax2,e_)+tmp2(ixImin1:ixImax1,ixImin2:ixImax2)*qdt  */

for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {

 pG->U[i3][i2][i1].E  += dt*(tmp2[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]);

	}
	}
	}





	if (tmprhor != NULL) free(tmprhor);
	if (tmprhol != NULL) free(tmprhol);
	if (tmpvr != NULL) free(tmpvr);
	if (tmpvl != NULL) free(tmpvl);
	if (tmpr != NULL) free(tmpr);
	if (tmpl != NULL) free(tmpl);
	if (tmp != NULL) free(tmp);
	if (tmp2 != NULL) free(tmp2);
	if (fieldf != NULL) free(fieldf);
	if (fieldl != NULL) free(fieldl);



	return;
}



static void hyperdifmomsourcene(int kf,int nf, int ii,int ii0, Real dt,ConsS ***Uint, GridS *pG)
{




//ii maps to the dimension - (k below in sac code)

// ii0 maps to the field id



int dim=kf;
int fieldi=nf;


	Real ***tmprhor=NULL, ***tmprhol=NULL, ***tmprhoc=NULL, ***tmpvr=NULL, ***tmpvl=NULL, ***tmpr=NULL, ***tmpl=NULL, ***tmp=NULL, ***tmp2=NULL, ***tmpc=NULL, ***fieldd=NULL, ***fieldl=NULL, ***fieldf=NULL;
        Real maxt1,maxt2;
	Real nur,nul;

	int n1z,n2z,n3z;
        int i,j,k;

	int il,iu;
	int jl,ju;
	int kl,ku;

	int is,ie,js,je,ks,ke;

        int i1,i2,i3;
        int iss,jss,kss;

	//int fieldi=energy;
        Real dtodx1 = pG->dt/pG->dx1, dtodx2 = pG->dt/pG->dx2, dtodx3 = pG->dt/pG->dx3;

	is = pG->is,
	ie = pG->ie;
	js = pG->js,
	je = pG->je;
	ks = pG->ks,
	ke = pG->ke;
        Real dx;


        /*rho, mom1, mom2, mom3, energy, b1, b2, b3*/

	/* With particles, one more ghost cell must be updated in predict step */
	#ifdef PARTICLES
	  Real d1;
	  il = is - 3;
	  iu = ie + 3;
	  jl = js - 3;
	  ju = je + 3;
	  kl = ks - 3;
	  ku = ke + 3;
	#else
	  il = is - 2;
	  iu = ie + 2;
	  jl = js - 2;
	  ju = je + 2;
	  kl = ks - 2;
	  ku = ke + 2;
	#endif

        kl=0;
        ku=0;

        dx= (pG->dx1)*(dim==1)+(pG->dx2)*(dim==2)+(pG->dx3)*(dim==3);


	if (pG->Nx[0] > 1)
		n1z = pG->Nx[0] + 2*nghost;
	else
		n1z = 1;

	if (pG->Nx[1] > 1)
		n2z = pG->Nx[1] + 2*nghost;
	else
		n2z = 1;

	if (pG->Nx[2] > 1)
		n3z = pG->Nx[2] + 2*nghost;
	else
		n3z = 1;


	//Real ***tmprhor=NULL, ***tmprhol=NULL, ***tmpvr=NULL, ***tmpvl=NULL, ***tmpr=NULL, ***tmpl=NULL, ***tmpc=NULL, ***tmp=NULL, ***tmp2=NULL, ***fieldf=NULL, ***fieldl=NULL;


switch(dim)
{
case 1:
	fieldf = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	fieldl = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmpr = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmpl = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmpc = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmprhor = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmprhol = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmprhoc = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmpvr = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmpvl = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmp = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmp2 = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
break;
case 2:
	fieldf = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	fieldl = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmpr = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmpl = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmpc = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmprhoc = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmprhor = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmprhol = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmpvr = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmpvl = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmp = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmp2 = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
break;
case 3:
	fieldf = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	fieldl = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmpr = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmpl = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmpc = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmprhoc = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmprhor = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmprhol = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmpvr = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmpvl = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmp = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmp2 = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
break;

}

switch(fieldi)
{

case 1:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			fieldf[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=Uinit[i3][i2][i1].M1;
					}
				}
			}
break;

case 2:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			fieldf[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=Uinit[i3][i2][i1].M2;
					}
				}
			}
break;

case 3:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			fieldf[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=Uinit[i3][i2][i1].M3;
					}
				}
			}
break;


}


switch(ii0)
{

case 1:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			fieldl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=Uinit[i3][i2][i1].M1;
					}
				}
			}
break;

case 2:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			fieldl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=Uinit[i3][i2][i1].M2;
					}
				}
			}
break;

case 3:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			fieldl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=Uinit[i3][i2][i1].M3;
					}
				}
			}
break;


}













 /* tmprhoC(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,ixImin2:ixImax2,&
     rho_)+w(ixImin1:ixImax1,ixImin2:ixImax2,rhob_)*/
for (i3=(kl); i3<=(ku-(dim==3)); i3++) {
    for (i2=(jl); i2<=(ju-(dim==2)); i2++) {
    	for (i1=(il); i1<=(iu-(dim==1)); i1++) {

			tmprhoc[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=(Uinit[i3][i2][i1].d+Uinit[i3][i2][i1].db);
					}
				}
			}

       /* tmp(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,&
           ixImin2:ixImax2,m0_+l)/(w(ixImin1:ixImax1,ixImin2:ixImax2,rho_)&
           +w(ixImin1:ixImax1,ixImin2:ixImax2,rhob_))*/
for (i3=(kl); i3<=(ku-(dim==3)); i3++) {
    for (i2=(jl); i2<=(ju-(dim==2)); i2++) {
    	for (i1=(il); i1<=(iu-(dim==1)); i1++) {

			tmp[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=fieldl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]/(Uinit[i3][i2][i1].d+Uinit[i3][i2][i1].db);
					}
				}
			}



            /*  tmpVL(ixmin1:ixmax1,ixmin2:ixmax2)=(w(ixmin1:ixmax1,&
                 ixmin2:ixmax2,m0_+ii0)+w(hxmin1:hxmax1,hxmin2:hxmax2,m0_&
                 +ii0))/two   */
 for (i3=(kl+(dim==3)); i3<=(ku); i3++) {
    for (i2=(jl+(dim==2)); i2<=(ju); i2++) {
    	for (i1=(il+(dim==1)); i1<=(iu); i1++) {

			tmpvl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=fieldf[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]+fieldf[AIN3(i1-(dim==1),i2-(dim==2),i3-(dim==3),dim)][AIN2((dim==1),i2-(dim==2),i3-(dim==3),dim)][AIN1((dim==1),i2-(dim==2),i3-(dim==3),dim)]/2;
					}
				}
			}

             /* tmpVR(ixmin1:ixmax1,ixmin2:ixmax2)=(w(jxmin1:jxmax1,&
                 jxmin2:jxmax2,m0_+ii0)+w(ixmin1:ixmax1,ixmin2:ixmax2,m0_&
                 +ii0))/two */
 for (i3=(kl); i3<=(ku-(dim==3)); i3++) {
    for (i2=(jl); i2<=(ju-(dim==2)); i2++) {
    	for (i1=(il); i1<=(iu-(dim==1)); i1++) {

			tmpvr[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=fieldf[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]+fieldf[AIN3(i1+(dim==1),i2+(dim==2),i3+(dim==3),dim)][AIN2((dim==1),i2+(dim==2),i3+(dim==3),dim)][AIN1((dim==1),i2+(dim==2),i3+(dim==3),dim)]/2;
					}
				}
			}














              /*CALL gradient1(tmp,ixmin1,ixmin2,ixmax1,ixmax2,k,tmp2)*/

switch(dim)
{

case 1:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1(tmp[i3][i2], n1z,pG->dx1,tmp2[i3][i2]);
				}
			}
break;

case 2:
 for (i3=kl; i3<=ku; i3++) {
    for (i1=il; i1<=iu; i1++) {
			gradient1(tmp[i3][i1], n2z,pG->dx2,tmp2[i3][i1]);
				}
			}
break;

case 3:
 for (i1=il; i1<=iu; i1++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1(tmp[i1][i2], n3z,pG->dx3,tmp2[i1][i2]);
				}
			}
break;

}


		/*tmp2(ixImin1:ixImax1,ixImin2:ixImax2)=tmp2(ixImin1:ixImax1,&
                 ixImin2:ixImax2)*(nuL(ixImin1:ixImax1,ixImin2:ixImax2)&
                 +nuR(ixImin1:ixImax1,ixImin2:ixImax2)+two*nushk&
                 (ixImin1:ixImax1,ixImin2:ixImax2,k))/two/two  */
 for (i3=(kl); i3<=(ku-(dim==3)); i3++) {
    for (i2=(jl); i2<=(ju-(dim==2)); i2++) {
    	for (i1=(il); i1<=(iu-(dim==1)); i1++) {

			tmp2[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=(tmp2[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]*(pG->Hv[i3][i2][i1].hdnul[dim][fieldi]))/2/2;
					}
				}
			}


             /* tmp(ixImin1:ixImax1,ixImin2:ixImax2)=tmprhoC(ixImin1:ixImax1,&
                 ixImin2:ixImax2)*tmp2(ixImin1:ixImax1,ixImin2:ixImax2)*/
for (i3=(kl); i3<=(ku-(dim==3)); i3++) {
    for (i2=(jl); i2<=(ju-(dim==2)); i2++) {
    	for (i1=(il); i1<=(iu-(dim==1)); i1++) {

			tmp[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=(tmprhoc[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)])*(tmp2[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]);
					}
				}
			}



             /* CALL gradient1(tmp,ixmin1,ixmin2,ixmax1,ixmax2,i,tmpC)*/
switch(dim)
{

case 1:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1(tmp[i3][i2], n1z,pG->dx1,tmpc[i3][i2]);
				}
			}
break;

case 2:
 for (i3=kl; i3<=ku; i3++) {
    for (i1=il; i1<=iu; i1++) {
			gradient1(tmp[i3][i1], n2z,pG->dx2,tmpc[i3][i1]);
				}
			}
break;

case 3:
 for (i1=il; i1<=iu; i1++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1(tmp[i1][i2], n3z,pG->dx3,tmpc[i1][i2]);
				}
			}
break;

}


		/*  wnew(ixImin1:ixImax1,ixImin2:ixImax2,m0_+ii0)&
                 =wnew(ixImin1:ixImax1,ixImin2:ixImax2,m0_+ii0)&
                 +tmpC(ixImin1:ixImax1,ixImin2:ixImax2)*qdt*/
for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {

switch(fieldi)
{
case 1:
     pG->U[i3][i2][i1].M1  += dt*(tmpc[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]);
break;

case 2:
     pG->U[i3][i2][i1].M2  += dt*(tmpc[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]);
break;

case 3:
     pG->U[i3][i2][i1].M3  += dt*(tmpc[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]);
break;

}


	}
	}
	}







             /* tmp(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,&
                 ixImin2:ixImax2,m0_+ii0)*tmp2(ixImin1:ixImax1,&
                 ixImin2:ixImax2)*/


for (i3=(kl); i3<=(ku-(dim==3)); i3++) {
    for (i2=(jl); i2<=(ju-(dim==2)); i2++) {
    	for (i1=(il); i1<=(iu-(dim==1)); i1++) {

			tmp[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=(fieldf[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]*(tmp2[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]));
					}
				}
			}



             /* CALL gradient1(tmp,ixmin1,ixmin2,ixmax1,ixmax2,i,tmpC)*/
switch(dim)
{

case 1:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1(tmp[i3][i2], n1z,pG->dx1,tmpc[i3][i2]);
				}
			}
break;

case 2:
 for (i3=kl; i3<=ku; i3++) {
    for (i1=il; i1<=iu; i1++) {
			gradient1(tmp[i3][i1], n2z,pG->dx2,tmpc[i3][i1]);
				}
			}
break;

case 3:
 for (i1=il; i1<=iu; i1++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1(tmp[i1][i2], n3z,pG->dx3,tmpc[i1][i2]);
				}
			}
break;

}



           /*   wnew(ixImin1:ixImax1,ixImin2:ixImax2,e_)=wnew(ixImin1:ixImax1,&
                 ixImin2:ixImax2,e_)+tmpC(ixImin1:ixImax1,ixImin2:ixImax2)*qdt*/
for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {

 pG->U[i3][i2][i1].E  += dt*(tmpc[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]);

	}
	}
	}

	if (tmprhoc != NULL) free(tmprhoc);
	if (tmprhor != NULL) free(tmprhor);
	if (tmprhol != NULL) free(tmprhol);
	if (tmpvr != NULL) free(tmpvr);
	if (tmpvl != NULL) free(tmpvl);
	if (tmpr != NULL) free(tmpr);
	if (tmpl != NULL) free(tmpl);
	if (tmpc != NULL) free(tmpc);
	if (tmp != NULL) free(tmp);
	if (tmp2 != NULL) free(tmp2);
	if (fieldf != NULL) free(fieldf);
	if (fieldl != NULL) free(fieldl);




	return;
}



#ifdef MHD
static void hyperdifbsource(int kf,int lf,int jj,int ii0, int mm, Real sb, Real dt, ConsS ***Uint, GridS *pG)
{




//ii maps to the dimension - (k below in sac code)

// ii0 maps to the field id

int dim=kf;
int fieldi=jj;


	//Real ***tmprhor=NULL, ***tmprhol=NULL, ***tmpvr=NULL, ***tmpvl=NULL, ***tmpr=NULL, ***tmpl=NULL, ***tmp=NULL, ***tmp2=NULL, ***fieldd=NULL;
        Real maxt1,maxt2;
	Real nur,nul;

	int n1z,n2z,n3z;
        int i,j,k;

	int il,iu;
	int jl,ju;
	int kl,ku;

	int is,ie,js,je,ks,ke;

        int i1,i2,i3;
        int iss,jss,kss;

	//fieldi=energy;
        Real dtodx1 = pG->dt/pG->dx1, dtodx2 = pG->dt/pG->dx2, dtodx3 = pG->dt/pG->dx3;

	is = pG->is,
	ie = pG->ie;
	js = pG->js,
	je = pG->je;
	ks = pG->ks,
	ke = pG->ke;
        Real dx;


        /*rho, mom1, mom2, mom3, energy, b1, b2, b3*/

	/* With particles, one more ghost cell must be updated in predict step */
	#ifdef PARTICLES
	  Real d1;
	  il = is - 3;
	  iu = ie + 3;
	  jl = js - 3;
	  ju = je + 3;
	  kl = ks - 3;
	  ku = ke + 3;
	#else
	  il = is - 2;
	  iu = ie + 2;
	  jl = js - 2;
	  ju = je + 2;
	  kl = ks - 2;
	  ku = ke + 2;
	#endif

        kl=0;
        ku=0;

        dx= (pG->dx1)*(dim==1)+(pG->dx2)*(dim==2)+(pG->dx3)*(dim==3);


	if (pG->Nx[0] > 1)
		n1z = pG->Nx[0] + 2*nghost;
	else
		n1z = 1;

	if (pG->Nx[1] > 1)
		n2z = pG->Nx[1] + 2*nghost;
	else
		n2z = 1;

	if (pG->Nx[2] > 1)
		n3z = pG->Nx[2] + 2*nghost;
	else
		n3z = 1;


	Real  ***tmpr=NULL, ***tmpl=NULL, ***tmp=NULL, ***tmp2=NULL, ***fieldr=NULL,***fieldl=NULL;


switch(dim)
{
case 1:
	fieldl = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	fieldr = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmpr = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmpl = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmp = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmp2 = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
break;
case 2:
	fieldr = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	fieldl = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmpr = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmpl = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmp = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmp2 = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
break;
case 3:
	fieldr = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	fieldl = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmpr = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmpl = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmp = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmp2 = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
break;

}



             //    tmp(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,&
             //       ixImin2:ixImax2,b0_+l)





          /*      tmpBL(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)&
                    =(w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b0_+j)&
                    +w(hxmin1:hxmax1,hxmin2:hxmax2,hxmin3:hxmax3,b0_+j))/two
                 tmpBR(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3)&
                    =(w(jxmin1:jxmax1,jxmin2:jxmax2,jxmin3:jxmax3,b0_+j)&
                    +w(ixmin1:ixmax1,ixmin2:ixmax2,ixmin3:ixmax3,b0_+j))/two    */

    /*j is + h is -*/
switch(jj)
{

case 1:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			fieldr[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=(Uinit[i3][i2][i1].B1c+Uinit[i3+(dim==3)][i2+(dim==2)][i1+(dim==1)].B1c)/2.0;
					}
				}
			}
break;

case 2:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			fieldr[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=(Uinit[i3][i2][i1].B2c+Uinit[i3+(dim==3)][i2+(dim==2)][i1+(dim==1)].B2c)/2.0;
					}
				}
			}
break;

case 3:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			fieldr[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=(Uinit[i3][i2][i1].B3c+Uinit[i3+(dim==3)][i2+(dim==2)][i1+(dim==1)].B3c)/2.0;
					}
				}
			}
break;


}




    /*j is + h is -*/
switch(jj)
{

case 1:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			fieldl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=(Uinit[i3][i2][i1].B1c+Uinit[i3-(dim==3)][i2-(dim==2)][i1-(dim==1)].B1c)/2.0;
					}
				}
			}
break;

case 2:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			fieldl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=(Uinit[i3][i2][i1].B2c+Uinit[i3-(dim==3)][i2-(dim==2)][i1-(dim==1)].B2c)/2.0;
					}
				}
			}
break;

case 3:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			fieldl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=(Uinit[i3][i2][i1].B3c+Uinit[i3-(dim==3)][i2-(dim==2)][i1-(dim==1)].B3c)/2.0;
					}
				}
			}
break;


}



          /*       tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
                    =w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,b0_+l)   */


 switch(lf)
{

case 1:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			tmp[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=Uinit[i3][i2][i1].B1c;
					}
				}
			}
break;

case 2:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			tmp[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=Uinit[i3][i2][i1].B2c;
					}
				}
			}
break;

case 3:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			tmp[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=Uinit[i3][i2][i1].B3c;
					}
				}
			}
break;


}



           /*      CALL gradient1L(tmp,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,&
                    ixmax3,k,tmp2)*/
 switch(dim)
{

case 1:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1l(tmp[i3][i2], n1z,pG->dx1,tmp2[i3][i2]);
				}
			}
break;

case 2:
 for (i3=kl; i3<=ku; i3++) {
    for (i1=il; i1<=iu; i1++) {
			gradient1l(tmp[i3][i1], n2z,pG->dx2,tmp2[i3][i1]);
				}
			}
break;

case 3:
 for (i1=il; i1<=iu; i1++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1l(tmp[i1][i2], n3z,pG->dx3,tmp2[i1][i2]);
				}
			}
break;

}



            /*     tmpL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
                    =(nuL(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))&
                    *tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)*/
    for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			tmpl[AIN3(i1,i2,i3,kf)][AIN2(i1,i2,i3,kf)][AIN1(i1,i2,i3,kf)]=tmp2[AIN3(i1,i2,i3,kf)][AIN2(i1,i2,i3,kf)][AIN1(i1,i2,i3,kf)]*(pG->Hv[i3][i2][i1].hdnul[kf][lf]);

					}
				}
			}




           /*      CALL gradient1R(tmp,ixmin1,ixmin2,ixmin3,ixmax1,ixmax2,&
                    ixmax3,k,tmp2)*/
  switch(dim)
{

case 1:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1r(tmp[i3][i2], n1z,pG->dx1,tmp2[i3][i2]);
				}
			}
break;

case 2:
 for (i3=kl; i3<=ku; i3++) {
    for (i1=il; i1<=iu; i1++) {
			gradient1r(tmp[i3][i1], n2z,pG->dx2,tmp2[i3][i1]);
				}
			}
break;

case 3:
 for (i1=il; i1<=iu; i1++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1r(tmp[i1][i2], n3z,pG->dx3,tmp2[i1][i2]);
				}
			}
break;

}



               /*  tmpR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
                    =(nuR(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3))&
                    *tmp2(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)*/

    for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			tmpr[AIN3(i1,i2,i3,kf)][AIN2(i1,i2,i3,kf)][AIN1(i1,i2,i3,kf)]=tmp2[AIN3(i1,i2,i3,kf)][AIN2(i1,i2,i3,kf)][AIN1(i1,i2,i3,kf)]*(pG->Hv[i3][i2][i1].hdnur[kf][lf]);

					}
				}
			}




              /*   wnew(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,b0_&
                    +ii0)=wnew(ixImin1:ixImax1,ixImin2:ixImax2,&
                    ixImin3:ixImax3,b0_+ii0)+sB*(tmpR(ixImin1:ixImax1,&
                    ixImin2:ixImax2,ixImin3:ixImax3)-tmpL(ixImin1:ixImax1,&
                    ixImin2:ixImax2,ixImin3:ixImax3))/dx(ixImin1:ixImax1,&
                    ixImin2:ixImax2,ixImin3:ixImax3,k)*qdt*/

for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {

switch(ii0)
{
case 1:
     pG->U[i3][i2][i1].B1c  += sb*(dtodx1*(dim==1)+dtodx2*(dim==2)+dtodx3*(dim==3))*(tmpr[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]-tmpl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]);
break;

case 2:
     pG->U[i3][i2][i1].B2c  += sb*(dtodx1*(dim==1)+dtodx2*(dim==2)+dtodx3*(dim==3))*(tmpr[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]-tmpl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]);
break;

case 3:
     pG->U[i3][i2][i1].B3c  += sb*(dtodx1*(dim==1)+dtodx2*(dim==2)+dtodx3*(dim==3))*(tmpr[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]-tmpl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]);
break;

}

    	}}}


               /*  wnew(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,e_)&
                    =wnew(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,e_)&
                    +sB*(tmpR(ixImin1:ixImax1,ixImin2:ixImax2,&
                    ixImin3:ixImax3)*tmpBR(ixImin1:ixImax1,ixImin2:ixImax2,&
                    ixImin3:ixImax3)-tmpL(ixImin1:ixImax1,ixImin2:ixImax2,&
                    ixImin3:ixImax3)*tmpBL(ixImin1:ixImax1,ixImin2:ixImax2,&
                    ixImin3:ixImax3))/dx(ixImin1:ixImax1,ixImin2:ixImax2,&
                    ixImin3:ixImax3,k)*qdt   */


for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {

 pG->U[i3][i2][i1].E  += sb*(dtodx1*(dim==1)+dtodx2*(dim==2)+dtodx3*(dim==3))*((tmpr[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]*fieldl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)])-(tmpl[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]*fieldr[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]));

	}
	}
	}

               //upto ear 15/2/2019






	if (tmpr != NULL) free(tmpr);
	if (tmpl != NULL) free(tmpl);
	if (tmp != NULL) free(tmp);
	if (tmp2 != NULL) free(tmp2);
	if (fieldl != NULL) free(fieldl);
	if (fieldr != NULL) free(fieldr);



	return;





    	}

//static void hyperdifbsourcene(int kf,int lf,int ii,int ii0,  int mm, Real sb,  Real dt, ConsS ***Uint, GridS *pG);
#endif /* MHD */

#ifdef MHD



static void hyperdifbsourcene(int kf,int lf,int ii,int ii0, int mm, Real sb, Real dt,ConsS ***Uint, GridS *pG)
{


//ii maps to the dimension - (k below in sac code)

// ii0 maps to the field id

int dim=kf;
int fieldi=lf;
int jj=ii;


	//Real ***tmprhor=NULL, ***tmprhol=NULL, ***tmpvr=NULL, ***tmpvl=NULL, ***tmpr=NULL, ***tmpl=NULL, ***tmp=NULL, ***tmp2=NULL, ***fieldd=NULL;
        Real maxt1,maxt2;
	Real nur,nul;

	int n1z,n2z,n3z;
        int i,j,k;

	int il,iu;
	int jl,ju;
	int kl,ku;

	int is,ie,js,je,ks,ke;

        int i1,i2,i3;
        int iss,jss,kss;

	//fieldi=energy;
        Real dtodx1 = pG->dt/pG->dx1, dtodx2 = pG->dt/pG->dx2, dtodx3 = pG->dt/pG->dx3;

	is = pG->is,
	ie = pG->ie;
	js = pG->js,
	je = pG->je;
	ks = pG->ks,
	ke = pG->ke;
        Real dx;


        /*rho, mom1, mom2, mom3, energy, b1, b2, b3*/

	/* With particles, one more ghost cell must be updated in predict step */
	#ifdef PARTICLES
	  Real d1;
	  il = is - 3;
	  iu = ie + 3;
	  jl = js - 3;
	  ju = je + 3;
	  kl = ks - 3;
	  ku = ke + 3;
	#else
	  il = is - 2;
	  iu = ie + 2;
	  jl = js - 2;
	  ju = je + 2;
	  kl = ks - 2;
	  ku = ke + 2;
	#endif

        kl=0;
        ku=0;

        dx= (pG->dx1)*(dim==1)+(pG->dx2)*(dim==2)+(pG->dx3)*(dim==3);


	if (pG->Nx[0] > 1)
		n1z = pG->Nx[0] + 2*nghost;
	else
		n1z = 1;

	if (pG->Nx[1] > 1)
		n2z = pG->Nx[1] + 2*nghost;
	else
		n2z = 1;

	if (pG->Nx[2] > 1)
		n3z = pG->Nx[2] + 2*nghost;
	else
		n3z = 1;


	Real  ***tmpc=NULL,  ***tmp=NULL, ***tmp2=NULL, ***fieldr=NULL,***fieldl=NULL;


switch(dim)
{
case 1:
	tmpc = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmp = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmp2 = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
break;
case 2:
	tmpc = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmp = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
	tmp2 = (Real***)calloc_3d_array(n3z, n1z, n2z, sizeof(Real));
break;
case 3:
	tmpc = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmp = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
	tmp2 = (Real***)calloc_3d_array(n1z, n2z, n3z, sizeof(Real));
break;

}




    /*j is + h is -*/

 /*tmp(ixImin1:ixImax1,ixImin2:ixImax2)=w(ixImin1:ixImax1,&
                    ixImin2:ixImax2,b0_+l)*/





          /*       tmp(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3)&
                    =w(ixImin1:ixImax1,ixImin2:ixImax2,ixImin3:ixImax3,b0_+l)   */


 switch(lf)
{

case 1:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			tmp[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=Uinit[i3][i2][i1].B1c;
					}
				}
			}
break;

case 2:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			tmp[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=Uinit[i3][i2][i1].B2c;
					}
				}
			}
break;

case 3:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			tmp[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=Uinit[i3][i2][i1].B3c;
					}
				}
			}
break;


}






            /*    CALL gradient1(tmp,ixmin1,ixmin2,ixmax1,ixmax2,k,tmp2)  */

 switch(kf)
{

case 1:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1l(tmp[i3][i2], n1z,pG->dx1,tmp2[i3][i2]);
				}
			}
break;

case 2:
 for (i3=kl; i3<=ku; i3++) {
    for (i1=il; i1<=iu; i1++) {
			gradient1l(tmp[i3][i1], n2z,pG->dx2,tmp2[i3][i1]);
				}
			}
break;

case 3:
 for (i1=il; i1<=iu; i1++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1l(tmp[i1][i2], n3z,pG->dx3,tmp2[i1][i2]);
				}
			}
break;

}




             /*    tmp2(ixImin1:ixImax1,ixImin2:ixImax2)=tmp2(ixImin1:ixImax1,&
                    ixImin2:ixImax2)*(nuL(ixImin1:ixImax1,ixImin2:ixImax2)&
                    +nuR(ixImin1:ixImax1,ixImin2:ixImax2))/two */


 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			tmp2[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=tmp2[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]*(pG->Hv[i3][i2][i1].hdnul[dim][fieldi]+pG->Hv[i3][i2][i1].hdnur[dim][fieldi])/2.0;

					}
				}
			}





                /* CALL gradient1(tmp2,ixmin1,ixmin2,ixmax1,ixmax2,m,tmpC)*/
switch(dim)
{

case 1:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1(tmp2[i3][i2], n1z,pG->dx1,tmpc[i3][i2]);
				}
			}
break;

case 2:
 for (i3=kl; i3<=ku; i3++) {
    for (i1=il; i1<=iu; i1++) {
			gradient1(tmp2[i3][i1], n2z,pG->dx2,tmpc[i3][i1]);
				}
			}
break;

case 3:
 for (i1=il; i1<=iu; i1++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1(tmp2[i1][i2], n3z,pG->dx3,tmpc[i1][i2]);
				}
			}
break;

}





               /*  wnew(ixImin1:ixImax1,ixImin2:ixImax2,b0_+ii0)&
                    =wnew(ixImin1:ixImax1,ixImin2:ixImax2,b0_+ii0)&
                    +sB*tmpC(ixImin1:ixImax1,ixImin2:ixImax2)*qdt*/


for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {

switch(ii0)
{
case 1:
     pG->U[i3][i2][i1].B1c  += sb*(dtodx1*(dim==1)+dtodx2*(dim==2)+dtodx3*(dim==3))*(tmpc[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]);
break;

case 2:
    pG->U[i3][i2][i1].B2c  += sb*(dtodx1*(dim==1)+dtodx2*(dim==2)+dtodx3*(dim==3))*(tmpc[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]);
break;

case 3:
    pG->U[i3][i2][i1].B3c  += sb*(dtodx1*(dim==1)+dtodx2*(dim==2)+dtodx3*(dim==3))*(tmpc[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]);
break;

}

    	}}}









               /*  tmp2(ixImin1:ixImax1,ixImin2:ixImax2)=tmp2(ixImin1:ixImax1,&
                    ixImin2:ixImax2)*w(ixImin1:ixImax1,ixImin2:ixImax2,b0_+j)*/
switch(jj)
{

case 1:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			tmp2[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=tmp2[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]*Uinit[i3][i2][i1].B1c;
					}
				}
			}
break;

case 2:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			tmp2[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=tmp2[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]*Uinit[i3][i2][i1].B2c;
					}
				}
			}
break;

case 3:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			tmp2[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]=tmp2[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]*Uinit[i3][i2][i1].B3c;
					}
				}
			}
break;


}









              /*   CALL gradient1(tmp2,ixmin1,ixmin2,ixmax1,ixmax2,m,tmpC)*/
switch(mm)
{

case 1:
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1(tmp2[i3][i2], n1z,pG->dx1,tmpc[i3][i2]);
				}
			}
break;

case 2:
 for (i3=kl; i3<=ku; i3++) {
    for (i1=il; i1<=iu; i1++) {
			gradient1(tmp2[i3][i1], n2z,pG->dx2,tmpc[i3][i1]);
				}
			}
break;

case 3:
 for (i1=il; i1<=iu; i1++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1(tmp2[i1][i2], n3z,pG->dx3,tmpc[i1][i2]);
				}
			}
break;

}







              /*   wnew(ixImin1:ixImax1,ixImin2:ixImax2,e_)=wnew&
                    (ixImin1:ixImax1,ixImin2:ixImax2,e_)+sB&
                    *tmpC(ixImin1:ixImax1,ixImin2:ixImax2)*qdt  */

for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {

 pG->U[i3][i2][i1].E  += sb*(dtodx1*(dim==1)+dtodx2*(dim==2)+dtodx3*(dim==3))*((tmpc[AIN3(i1,i2,i3,dim)][AIN2(i1,i2,i3,dim)][AIN1(i1,i2,i3,dim)]));

	}
	}
	}




	if (tmpc != NULL) free(tmpc);
	if (tmp != NULL) free(tmp);
	if (tmp2 != NULL) free(tmp2);



	return;

}


#endif /* MHD */


#endif /* SAC_INTEGRATOR */
