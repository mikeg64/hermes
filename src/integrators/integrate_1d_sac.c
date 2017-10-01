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
 *   - For adb hydro, requires (9*Cons1DS +  3*Real) = 48 3D arrays
 *   - For adb mhd, requires   (9*Cons1DS + 10*Real) = 73 3D arrays
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

/* 1D scratch vectors used by lr_states and flux functions */
static Real *Bx=NULL, *Bxb=NULL, *Bxi=NULL;
static Prim1DS *W=NULL, *Wl=NULL, *Wr=NULL;
static Cons1DS *U1d=NULL;

/* Variables at t^{n+1/2} used for source terms */
static Real *dhalf = NULL, *phalf = NULL;

/* Variables needed for cylindrical coordinates */
#ifdef CYLINDRICAL
static Real *geom_src=NULL;
#endif



/*static void hyperdifviscr(int field,int dim,const GridS *pG);
static void hyperdifviscl(int field,int dim,const GridS *pG);
static void hyperdifrhosource(int field,int dim,Real dt,const GridS *pG);
static void hyperdifesource(int dim,Real dt,const GridS *pG); 
static void hyperdifmomsource(int field,int dim,int ii,int ii0,Real dt,const GridS *pG);
static void hyperdifmomsourcene(int field,int dim,int ii,int ii0, Real dt,const GridS *pG);*/

#ifdef MHD
//static void hyperdifbsource(int ii,int ii0, Real dt, const GridS *pG);
//static void hyperdifbsourcene(int ii,int ii0, Real dt, const GridS *pG);
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
  Real coolfl,coolfr,coolf,Eh=0.0;
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
int ii1, dim, ii, ii0;
int field; /*integers map to following index rho, mom1, mom2, energy, b1, b2,energyb,rhob,b1b,b2b*/


/*=== STEP 1: Compute L/R x1-interface states and 1D x1-Fluxes ===============*/

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
    Bx[i] = pG->U[ks][js][i].B1c;
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

#endif
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++) U1d[i].s[n] = pG->U[ks][js][i].s[n];
#endif
  }

 

/*--- Step 1b ------------------------------------------------------------------
 * Compute L and R states at X1-interfaces, add "MHD source terms" for 0.5*dt
 */


  for (i=is-nghost; i<=ie+nghost; i++) {
    W[i] = Cons1D_to_Prim1D(&U1d[i], &Bx[i],&Bxb[i]);
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
//       Wl[i].Vx -= hdt*gl;
//       Wr[i].Vx -= hdt*gr;
// #else
    /*  phicr = (*StaticGravPot)( x1             ,x2,x3);
      phicl = (*StaticGravPot)((x1-    pG->dx1),x2,x3);*/
      phifc = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

      W[i].Vx -= dtodx1*(phifc );

     /* Wl[i].Vx -= dtodx1*(phifc - phicl);
      Wr[i].Vx -= dtodx1*(phicr - phifc);*/
// #endif /* CYLINDRICAL */
    }
  }


/*--- Step 1c (cont) -----------------------------------------------------------
 * Add the geometric source-terms now using cell-centered primitive
 * variables at time t^n
 */

#ifdef SELF_GRAVITY
  for (i=il+1; i<=iu-1; i++) {
   /* Wl[i].Vx -= hdtodx1*(pG->Phi[ks][js][i] - pG->Phi[ks][js][i-1]);
    Wr[i].Vx -= hdtodx1*(pG->Phi[ks][js][i] - pG->Phi[ks][js][i-1]);*/
    W[i].Vx -= hdtodx1*(pG->Phi[ks][js][i+1] - pG->Phi[ks][js][i-1]);
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
      coolf = (*CoolingFunc)(W[i].d,W[i].P,(pG->dt));
      W[i].P -= pG->dt*Gamma_1*coolf;
      /*Wl[i].P -= 0.5*pG->dt*Gamma_1*coolfl;
      Wr[i].P -= 0.5*pG->dt*Gamma_1*coolfr;*/
    }
  }
#endif /* BAROTROPIC */


/*--- Step 1d ------------------------------------------------------------------
 * Compute 1D fluxes in x1-direction, storing into 3D array
 */
    for (i=il+1; i<=iu; i++) {
      Uc_x1[i] = Prim1D_to_Cons1D(&W[i],&Bx[i],&Bxb[i]);
      


      fluxes(Uc_x1[i],Uc_x1[i],W[i],W[i],*Bx,*Bxb,&x1Flux[i]);
    }
  

/*=== STEP 8: Compute cell-centered values at n+1/2 ==========================*/

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
      M1h -= hdtodx1*(phir-phil)*pG->U[ks][js][i].d;
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






/*--- Step 12c -----------------------------------------------------------------
 * Update cell-centered variables in pG using 3D x3-Fluxes
 */

/*--- Step 12a -----------------------------------------------------------------
 * Update cell-centered variables in pG using 1D x1-fluxes
 */

  for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
    rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
#endif
    pG->U[ks][js][i].d  -= dtodx1*(rsf*x1Flux[i+1].d  - lsf*x1Flux[i].d );
    pG->U[ks][js][i].M1 -= dtodx1*(rsf*x1Flux[i+1].Mx - lsf*x1Flux[i].Mx);
    pG->U[ks][js][i].M2 -= dtodx1*(SQR(rsf)*x1Flux[i+1].My - SQR(lsf)*x1Flux[i].My);
    pG->U[ks][js][i].M3 -= dtodx1*(rsf*x1Flux[i+1].Mz - lsf*x1Flux[i].Mz);
#ifndef BAROTROPIC
    pG->U[ks][js][i].E  -= dtodx1*(rsf*x1Flux[i+1].E  - lsf*x1Flux[i].E );
#endif /* BAROTROPIC */
#ifdef MHD
    pG->U[ks][js][i].B2c -= dtodx1*(x1Flux[i+1].By - x1Flux[i].By);
    pG->U[ks][js][i].B3c -= dtodx1*(rsf*x1Flux[i+1].Bz - lsf*x1Flux[i].Bz);
/* For consistency, set B2i and B3i to cell-centered values.  */
    pG->B2i[ks][js][i] = pG->U[ks][js][i].B2c;
    pG->B3i[ks][js][i] = pG->U[ks][js][i].B3c;
#endif /* MHD */
#if (NSCALARS > 0)
    for (n=0; n<NSCALARS; n++)
      pG->U[ks][js][i].s[n] -= dtodx1*(rsf*x1Flux[i+1].s[n] - lsf*x1Flux[i].s[n]);
#endif
  }

/*--- Step 12b: Not needed in 1D ---*/
/*--- Step 12c: Not needed in 1D ---*/
/*--- Step 12d: Not needed in 1D ---*/








//hyperdifvisc1r

//hyperdifvisc1l

//computec
//computemaxc

//density contribution
for(dim=0; dim<2; dim++) //each direction
{



//hyperdifvisc1ir
//hyperdifvisc1il
//hyperdifrhosource1 
;
}

//energy hyperdiffusion term
for(dim=0; dim<2; dim++) //each direction
{
//hyperdifvisc1ir
//hyperdifvisc1il
//hyperdifesource1 
;

}

       //momentum hyperdiffusion term
for(dim=0; dim<2; dim++) //each direction
{
//hyperdifvisc1ir
//hyperdifvisc1il
//hyperdifesource1 

		         for(ii1=0;ii1<=1;ii1++)
		         {
		                  if (ii1 == 0)
		                  {
				           ii=dim;
				           ii0=field;  //f is field
		                  }
		                  else
		                  {
				           ii=field;
				           ii0=dim;
		                   }

				  if(ii==dim)
				  ;//  hyperdifmomsource1(ii,ii0,pG->dt);
				  else
				   ;// hyperdifmomsourcene1(ii,ii0,pG->dt);  //off diagonal
		        }


}
#ifdef MHD

  //b field hyperdiffusion term
for(dim=0; dim<2; dim++) //each direction
{
//hyperdifvisc1ir
//hyperdifvisc1il
 

		         for(ii1=0;ii1<=1;ii1++)
		         {
		                  if (ii1 == 0)
		                  {
				           ii=dim;
				           ii0=field;  //f is field
		                  }
		                  else
		                  {
				           ii=field;
				           ii0=dim;
		                   }

				  if(ii==dim)
				   ;// hyperdifbsource(ii,ii0,pG->dt,pG);
				  else
				   ;// hyperdifbsourcene(ii,ii0,pG->dt,pG);  //off diagonal
		        }


}

#endif  /*hyperdiffusion source term for bfield*/


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

  /*if ((Bxc = (Real*)malloc(size1*sizeof(Real))) == NULL) goto on_error;*/
  if ((Bx = (Real*)malloc(size1*sizeof(Real))) == NULL) goto on_error;
  if ((Bxb = (Real*)malloc(size1*sizeof(Real))) == NULL) goto on_error;
  if ((Bxi = (Real*)malloc(size1*sizeof(Real))) == NULL) goto on_error;

  if ((U1d       =(Cons1DS*)malloc(size1*sizeof(Cons1DS)))==NULL) goto on_error;
  if ((Uc_x1 =(Cons1DS*)malloc(size1*sizeof(Cons1DS)))==NULL) goto on_error;
/*if ((Ur_x1Face =(Cons1DS*)malloc(size1*sizeof(Cons1DS)))==NULL) goto on_error;*/
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
  if (Bx != NULL) free(Bx);
  if (Bxb != NULL) free(Bxb);
  if (Bxi != NULL) free(Bxi);

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




/*

static void hyperdifviscr(int field,int dim,const GridS *pG)
{

	return;
}

static void hyperdifviscl(int field,int dim,const GridS *pG)
{

	return;
}

static void hyperdifrhosource(int field,int dim,Real dt,const GridS *pG)
{

	return;
}


static void hyperdifesource(int dim,Real dt,const GridS *pG)
{

	return;
}

 
static void hyperdifmomsource(int field,int dim,int ii,int ii0,Real dt,const GridS *pG)
{

	return;
}


static void hyperdifmomsourcene(int field,int dim,int ii,int ii0, Real dt,const GridS *pG)
{

	return;
}

*/

#ifdef MHD

/*static void hyperdifbsource(int ii,int ii0, Real dt, const GridS *pG)
{

	return;
}


static void hyperdifbsourcene(int ii,int ii0, Real dt, const GridS *pG)
{

	return;
}*/


#endif /* MHD */


#endif /* SAC_INTEGRATOR */
