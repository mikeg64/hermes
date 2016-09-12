#include "../copyright.h"
/*============================================================================*/
/*! \file integrate_1d_smaug.c
 *  \Compute fluxes using smaug . 
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
 * - Magnetohydrodynamic code for gravitationally-stratified media : SMAUG
 *   Shelyag et al 
 *   http://adsabs.harvard.edu/abs/2008A%26A...486..655S
 *
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - integrate_1d_smaug()
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

#if defined(SMAUG_INTEGRATOR)
#ifdef SPECIAL_RELATIVITY
#error : The SMAUG integrator cannot be used for special relativity.
#endif /* SPECIAL_RELATIVITY */



/* 1D scratch vectors used by lr_states and flux functions */
static Real *Bxc=NULL, *Bxi=NULL;
static Prim1DS *W=NULL;
static Cons1DS *U1d=NULL;



/* The L/R states of conserved variables and fluxes at each cell face */
static Cons1DS **Uc_x1=NULL;
static Cons1DS **Uc_x2=NULL;
static Cons1DS **x1Flux=NULL, **x2Flux=NULL;

static Cons1DS **Uc=NULL;


/* The interface magnetic fields and emfs */
#ifdef MHD
static Real **B1_x1=NULL, **B2_x2=NULL;

#endif /* MHD */

/* density and Pressure at t^{n+1/2} needed by MHD, cooling, and gravity */
static Real ***dhalf = NULL, ***phalf=NULL;

/* variables needed for H-correction of Sanders et al (1998) */
extern Real etah;
#ifdef H_CORRECTION
static Real ***eta1=NULL, ***eta2=NULL, ***eta3=NULL;
#endif

/* variables needed to conserve net Bz in shearing box */
#ifdef SHEARING_BOX
static ConsS **Flxiib=NULL, **Flxoib=NULL;
static ConsS **rFlxiib=NULL, **rFlxoib=NULL;
#endif



/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES: 
 *   integrate_emf1_corner() - the upwind CT method in GS05, for emf1
 *   integrate_emf2_corner() - the upwind CT method in GS05, for emf2
 *   integrate_emf3_corner() - the upwind CT method in GS05, for emf3
 *============================================================================*/


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
/*! \fn void integrate_3d_ctu(DomainS *pD)
 *  \brief 3D CTU integrator for MHD using 6-solve method */

void integrate_1d_smaug(DomainS *pD)
{
  GridS *pG=(pD->Grid);
  Real dtodx1=pG->dt/pG->dx1, dtodx2=pG->dt/pG->dx2, dtodx3=pG->dt/pG->dx3;
  Real hdt = 0.5*pG->dt, dx2=pG->dx2;
  Real q1 = 0.5*dtodx1, q2 = 0.5*dtodx2, q3 = 0.5*dtodx3;
  int dir;
  int i,il,iu, is = pG->is, ie = pG->ie;
  int j,jl,ju, js = pG->js, je = pG->je;
  int k,kl,ku, ks = pG->ks, ke = pG->ke;
  Real x1,x2,x3,phicl,phicr,phifc,phil,phir,phic,M1h,M2h,M3h,Bx=0.0,Bxb=0.0;


/*Used for hyperdiffusion computations*/
int ii1, dim, ii, ii0;
int field; /*integers map to following index rho, mom1, mom2, energy, b1, b2,energyb,rhob,b1b,b2b*/

#ifdef MHD
  Real MHD_src_By,MHD_src_Bz,mdb1,mdb2,mdb3;
  Real db1,db2,db3,l1,l2,l3,B1,B2,B3,V1,V2,V3;
  Real B1ch,B2ch,B3ch;
#endif
// #if defined(MHD) || defined(SELF_GRAVITY)
  Real dx1i=1.0/pG->dx1, dx2i=1.0/pG->dx2, dx3i=1.0/pG->dx3;
// #endif

#if (NSCALARS > 0)
  int n;
#endif
#ifdef SELF_GRAVITY
  Real gxl,gxr,gyl,gyr,gzl,gzr,flx_m1l,flx_m1r,flx_m2l,flx_m2r,flx_m3l,flx_m3r;
#endif
#ifdef FEEDBACK
  Real dt1 = 1.0/pG->dt;
#endif


#ifdef STATIC_MESH_REFINEMENT
  int ncg,npg,dim;
  int ii,ics,ice,jj,jcs,jce,kk,kcs,kce,ips,ipe,jps,jpe,kps,kpe;
#endif


  Real g,gl,gr;
  Real lsf=1.0, rsf=1.0;



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

/* Set etah=0 so first calls to flux functions do not use H-correction */
  etah = 0.0;

/* Compute predictor feedback from particle drag */
#ifdef FEEDBACK
  feedback_predictor(pD);
  exchange_gpcouple(pD,1);
#endif

/*=== STEP 1: Compute L/R x1-interface states and 1D x1-Fluxes ===============*/

/*--- Step 1a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M1, M2, M3, E, B2c, B3c, s[n])
 */

 
    for (j=jl; j<=ju; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        U1d[i].d  = pG->U[ks][j][i].d;
        U1d[i].Mx = pG->U[ks][j][i].M1;
        U1d[i].My = pG->U[ks][j][i].M2;
        U1d[i].Mz = pG->U[ks][j][i].M3;
#ifndef BAROTROPIC
        U1d[i].E  = pG->U[ks][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
        U1d[i].By = pG->U[ks][j][i].B2c;
        U1d[i].Bz = pG->U[ks][j][i].B3c;
        Bxc[i] = pG->U[ks][j][i].B1c;
        Bxi[i] = pG->B1i[ks][j][i];
        B1_x1[j][i] = pG->B1i[ks][j][i];       
#endif /* MHD */
#ifdef BKG
        U1d[i].db  = pG->U[ks][j][i].db;
        U1d[i].Byb = pG->U[ks][j][i].B2cb;
        U1d[i].Bzb = pG->U[ks][j][i].B3cb;

#endif
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) U1d[i].s[n] = pG->U[ks][j][i].s[n];
#endif
      }

/*--- Step 1b ------------------------------------------------------------------
 * Compute L and R states at X1-interfaces, add "MHD source terms" for 0.5*dt
 */

   for (i=is-nghost; i<=ie+nghost; i++) {
      W[i] = Cons1D_to_Prim1D(&U1d[i],&Bxc[i],&Bxb);

/*--- Step 1c ------------------------------------------------------------------
 * Add source terms from static gravitational potential for 0.5*dt to L/R states
 */

      /*if (StaticGravPot != NULL){
        for (i=il+1; i<=iu; i++) {
          cc_pos(pG,i,j,ks,&x1,&x2,&x3);

          phicr = (*StaticGravPot)( x1             ,x2,x3);
          phicl = (*StaticGravPot)((x1-    pG->dx1),x2,x3);
          phifc = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

          gl = 2.0*(phifc - phicl)*dx1i;
          gr = 2.0*(phicr - phifc)*dx1i;
#if defined(CYLINDRICAL) && defined(FARGO)
          gl -= r[i-1]*SQR((*OrbitalProfile)(r[i-1]));
          gr -= r[i  ]*SQR((*OrbitalProfile)(r[i  ]));
#endif

          W[i].Vx -= hdt*gl;
         
        }
      }*/



/*--- Step 1c (cont) -----------------------------------------------------------
 * Add the geometric source-terms now using cell-centered primitive
 * variables at time t^n
 */



/*--- Step 1d ------------------------------------------------------------------
 * Compute 1D fluxes in x1-direction, storing into 3D array
 */
    for (i=il+1; i<=iu; i++) {
      Uc_x1[j][i] = Prim1D_to_Cons1D(&W[i],&Bxi[i],&Bxb);
      

#ifdef MHD
      Bx = B1_x1[j][i];
      Bxb=0.0;//?????????????????????????
#endif
      fluxes(Uc_x1[j][i],Uc_x1[j][i],W[i],W[i],Bx,Bxb,&x1Flux[j][i]);
    }
  }

/*=== STEP 2: Compute L/R x2-interface states and 1D x2-Fluxes ===============*/

/*--- Step 2a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M2, M3, M1, E, B3c, B1c, s[n])
 */

  for (i=il; i<=iu; i++) {
#ifdef CYLINDRICAL
    dx2 = r[i]*pG->dx2;
    dx2i = 1.0/dx2;
    dtodx2 = pG->dt*dx2i;
    hdtodx2 = 0.5*dtodx2;
#endif
    for (j=js-nghost; j<=je+nghost; j++) {
      U1d[j].d  = pG->U[ks][j][i].d;
      U1d[j].Mx = pG->U[ks][j][i].M2;
      U1d[j].My = pG->U[ks][j][i].M3;
      U1d[j].Mz = pG->U[ks][j][i].M1;
#ifndef BAROTROPIC
      U1d[j].E  = pG->U[ks][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
      U1d[j].By = pG->U[ks][j][i].B3c;
      U1d[j].Bz = pG->U[ks][j][i].B1c;
      Bxc[j] = pG->U[ks][j][i].B2c;
      Bxi[j] = pG->B2i[ks][j][i];
      B2_x2[j][i] = pG->B2i[ks][j][i];
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) U1d[j].s[n] = pG->U[ks][j][i].s[n];
#endif
    }

/*--- Step 2b ------------------------------------------------------------------
 * Compute L and R states at X2-interfaces, add "MHD source terms" for 0.5*dt
 */



/*--- Step 2c ------------------------------------------------------------------
 * Add source terms from static gravitational potential for 0.5*dt to L/R states
 */

    /*  if (StaticGravPot != NULL){
        for (j=jl+1; j<=ju; j++) {
          cc_pos(pG,i,j,ks,&x1,&x2,&x3);
          phicr = (*StaticGravPot)(x1, x2             ,x3);
          phicl = (*StaticGravPot)(x1,(x2-    pG->dx2),x3);
          phifc = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

          W[j].Vx -= dtodx2*(phifc - phicl);
          
        }
      }*/



/*--- Step 2d ------------------------------------------------------------------
 * Compute 1D fluxes in x2-direction, storing into 3D array
 */

    for (j=jl+1; j<=ju; j++) {
      Uc_x2[j][i] = Prim1D_to_Cons1D(&W[j],&Bxi[j],&Bxb);
      
#ifdef MHD
      Bx = B2_x2[j][i];
      Bxb=0.0;//???????????????????????????????????
#endif
      fluxes(Uc_x2[j][i],Uc_x2[j][i],W[j],W[j],Bx,Bxb,&x2Flux[j][i]);
    }
  }


/*--- Step 3c ------------------------------------------------------------------
 * Add source terms from static gravitational potential for 0.5*dt to L/R states
 */

    /*  if (StaticGravPot != NULL){
        for (k=kl+1; k<=ku; k++) {
          cc_pos(pG,i,j,ks,&x1,&x2,&x3);
          phicr = (*StaticGravPot)(x1,x2, x3             );
          phicl = (*StaticGravPot)(x1,x2,(x3-    pG->dx3));
          phifc = (*StaticGravPot)(x1,x2,(x3-0.5*pG->dx3));

          W[k].Vx -= dtodx3*(phifc - phicl);
          
        }
      }*/



/*--- Step 3d ------------------------------------------------------------------
 * Compute 1D fluxes in x3-direction, storing into 3D array
 */

 /*Not needed here for 2d problem*/

/*=== STEP 4:  Update face-centered B for 0.5*dt =============================*/

/*--- Step 4a ------------------------------------------------------------------
 * Calculate the cell centered value of emf1,2,3 at t^{n} and integrate
 * to corner.
 */



/*--- Step 4b ------------------------------------------------------------------
 * Update the interface magnetic fields using CT for a half time step.
 */



/*=== STEP 5: Correct x1-interface states with transverse flux gradients =====*/

/*--- Step 5a ------------------------------------------------------------------
 * Correct x1-interface states using x2-fluxes computed in Step 2d.
 * Since the fluxes come from an x2-sweep, (x,y,z) on RHS -> (z,x,y) on LHS 
 */



/*--- Step 5b ------------------------------------------------------------------
 * Correct x1-interface states using x3-fluxes computed in Step 3d.
 * Since the fluxes come from an x3-sweep, (x,y,z) on RHS -> (y,z,x) on LHS
 */



/*--- Step 5c ------------------------------------------------------------------
 * Add the "MHD source terms" from the x2- and x3-flux-gradients to the
 * conservative variables on the x1Face.  Limiting is used as in GS (2007)
 */



/*--- Step 5d ------------------------------------------------------------------
 * Add source terms for a static gravitational potential arising from x2-Flux
 * and x3-Flux gradients.  To improve conservation of total energy, average
 * the energy source term computed at cell faces.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */

 /* if (StaticGravPot != NULL){
  for (k=kl+1; k<=ku-1; k++) {
    for (j=jl+1; j<=ju-1; j++) {
      for (i=il+1; i<=iu; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        phic = (*StaticGravPot)(x1, x2             ,x3);
        phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);    ---------*/

        /* correct right states; x2 and x3 gradients */

/*-----------------------------
#ifdef CYLINDRICAL
        q2 = hdt/(r[i]*pG->dx2);
#endif
        Ur_x1Face[k][j][i].My -= q2*(phir-phil)*pG->U[k][j][i].d;
#ifndef BAROTROPIC
        Ur_x1Face[k][j][i].E -= q2*(x2Flux[k][j  ][i  ].d*(phic - phil)
                                  + x2Flux[k][j+1][i  ].d*(phir - phic));
        #ifdef ROTATING_FRAME
                Ur_x1Face[k][j][i].E += hdt * 0.5*(x2Flux[k][j  ][i  ].d*sin(x2-0.5*pG->dx2) +
						   x2Flux[k][j+1][i  ].d*sin(x2+0.5*pG->dx2)) *SQR(Omega_0)*Rc;
        #endif
#endif

        phir = (*StaticGravPot)(x1,x2,(x3+0.5*pG->dx3));
        phil = (*StaticGravPot)(x1,x2,(x3-0.5*pG->dx3));

        Ur_x1Face[k][j][i].Mz -= q3*(phir-phil)*pG->U[k][j][i].d;
#ifndef BAROTROPIC
        Ur_x1Face[k][j][i].E -= q3*(x3Flux[k  ][j][i  ].d*(phic - phil)
                                  + x3Flux[k+1][j][i  ].d*(phir - phic));
#endif

-----------------------*/



        /* correct left states; x2 and x3 gradients */
/*--------------------------------------------
        phic = (*StaticGravPot)((x1-pG->dx1), x2             ,x3);
        phir = (*StaticGravPot)((x1-pG->dx1),(x2+0.5*pG->dx2),x3);
        phil = (*StaticGravPot)((x1-pG->dx1),(x2-0.5*pG->dx2),x3);

#ifdef CYLINDRICAL
        q2 = hdt/(r[i-1]*pG->dx2);
#endif
        Ul_x1Face[k][j][i].My -= q2*(phir-phil)*pG->U[k][j][i-1].d;
#ifndef BAROTROPIC
        Ul_x1Face[k][j][i].E -= q2*(x2Flux[k][j  ][i-1].d*(phic - phil)
                                  + x2Flux[k][j+1][i-1].d*(phir - phic));
        #ifdef ROTATING_FRAME
                Ul_x1Face[k][j][i].E += hdt * 0.5*(x2Flux[k][j  ][i-1].d*sin(x2-0.5*pG->dx2) +
                                                x2Flux[k][j+1][i-1].d*sin(x2+0.5*pG->dx2)) *SQR(Omega_0)*Rc;
        #endif
#endif

        phir = (*StaticGravPot)((x1-pG->dx1),x2,(x3+0.5*pG->dx3));
        phil = (*StaticGravPot)((x1-pG->dx1),x2,(x3-0.5*pG->dx3));

        Ul_x1Face[k][j][i].Mz -= q3*(phir-phil)*pG->U[k][j][i-1].d;
#ifndef BAROTROPIC
        Ul_x1Face[k][j][i].E -= q3*(x3Flux[k  ][j][i-1].d*(phic - phil)
                                  + x3Flux[k+1][j][i-1].d*(phir - phic));
#endif
      }
    }
  }}*/




/*=== STEP 6: Correct x2-interface states with transverse flux gradients =====*/

/*--- Step 6a ------------------------------------------------------------------
 * Correct x2-interface states using x1-fluxes computed in Step 1d.
 * Since the fluxes come from an x1-sweep, (x,y,z) on RHS -> (y,z,x) on LHS
 */



/*--- Step 6b ------------------------------------------------------------------
 * Correct x2-interface states using x3-fluxes computed in Step 3d.
 * Since the fluxes come from an x3-sweep, (x,y,z) on RHS -> (z,x,y) on LHS 
 */



/*--- Step 6c ------------------------------------------------------------------
 * Add the "MHD source terms" from the x1- and x3-flux-gradients to the
 * conservative variables on the x2Face.  Limiting is used as in GS (2007)
 */



/*--- Step 6d ------------------------------------------------------------------
 * Add source terms for a static gravitational potential arising from x1-Flux
 * and x3-Flux gradients. To improve conservation of total energy,
 * average the energy source term computed at cell faces.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */



/*=== STEP 7: Correct x3-interface states with transverse flux gradients =====*/

/*--- Step 7a ------------------------------------------------------------------
 * Correct x3-interface states using x1-fluxes computed in Step 1d.
 * Since the fluxes come from an x1-sweep, (x,y,z) on RHS -> (z,x,y) on LHS 
 */


/*--- Step 7b ------------------------------------------------------------------
 * Correct x3-interface states using x2-fluxes computed in Step 2d.
 * Since the fluxes come from an x2-sweep, (x,y,z) on RHS -> (y,z,x) on LHS 
 */



/*--- Step 7c ------------------------------------------------------------------
 * Add the "MHD source terms" from the x1- and x2-flux-gradients to the
 * conservative variables on the x3Face.  Limiting is used as in GS07.
 */



/*--- Step 7d ------------------------------------------------------------------
 * Add source terms for a static gravitational potential arising from x1-Flux
 * and x2-Flux gradients. To improve conservation of total energy,
 * average the energy source term computed at cell faces.
 *    S_{M} = -(\rho) Grad(Phi);   S_{E} = -(\rho v) Grad{Phi}
 */



/*--- Step 7e ------------------------------------------------------------------
 * Apply density floor
 */


/*=== STEP 8: Compute cell-centered values at n+1/2 ==========================*/


/*=== STEP 9: Compute 3D x1-Flux, x2-Flux, x3-Flux ===========================*/

/*--- Step 9a ------------------------------------------------------------------
 * Compute maximum wavespeeds in multidimensions (eta in eq. 10 from Sanders et
 *  al. (1998)) for H-correction
 */


/*--- Step 9b ------------------------------------------------------------------
 * Compute 3D x1-fluxes from corrected L/R states.
 */


/*--- Step 9c ------------------------------------------------------------------
 * Compute 3D x2-fluxes from corrected L/R states.
 */


/*--- Step 9d ------------------------------------------------------------------
 * Compute 3D x3-fluxes from corrected L/R states.
 */

 

/*=== STEP 10: Update face-centered B for a full timestep ====================*/

/*--- Step 10a -----------------------------------------------------------------
 * Integrate emf*^{n+1/2} to the grid cell corners
 */




/*--- Step 10b -----------------------------------------------------------------
 * Update the interface magnetic fields using CT for a full time step.
 */


/*=== STEP 11: Add source terms for a full timestep using n+1/2 states =======*/





/*=== STEP 12: Update cell-centered values for a full timestep ===============*/

/*--- Step 12a -----------------------------------------------------------------
 * Update cell-centered variables in pG using 3D x1-Fluxes
 */
 for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
      rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
#endif
      pG->U[ks][j][i].d  -= dtodx1*(rsf*x1Flux[j][i+1].d  - lsf*x1Flux[j][i].d );
      pG->U[ks][j][i].M1 -= dtodx1*(rsf*x1Flux[j][i+1].Mx - lsf*x1Flux[j][i].Mx);
      pG->U[ks][j][i].M2 -= dtodx1*(SQR(rsf)*x1Flux[j][i+1].My - SQR(lsf)*x1Flux[j][i].My);
      pG->U[ks][j][i].M3 -= dtodx1*(rsf*x1Flux[j][i+1].Mz - lsf*x1Flux[j][i].Mz);
#ifndef BAROTROPIC
      pG->U[ks][j][i].E  -= dtodx1*(rsf*x1Flux[j][i+1].E  - lsf*x1Flux[j][i].E );
#endif /* BAROTROPIC */
#ifdef MHD
      pG->U[ks][j][i].B2c -= dtodx1*(x1Flux[j][i+1].By - x1Flux[j][i].By);
      pG->U[ks][j][i].B3c -= dtodx1*(rsf*x1Flux[j][i+1].Bz - lsf*x1Flux[j][i].Bz);
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        pG->U[ks][j][i].s[n] -= dtodx1*(rsf*x1Flux[j][i+1].s[n] 
                                      - lsf*x1Flux[j][i  ].s[n]);
#endif
    }
  }

 

/*--- Step 12b -----------------------------------------------------------------
 * Update cell-centered variables in pG using 3D x2-Fluxes
 */

  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
      dtodx2 = pG->dt/(r[i]*pG->dx2);
#endif
      pG->U[ks][j][i].d  -= dtodx2*(x2Flux[j+1][i].d  - x2Flux[j][i].d );
      pG->U[ks][j][i].M1 -= dtodx2*(x2Flux[j+1][i].Mz - x2Flux[j][i].Mz);
      pG->U[ks][j][i].M2 -= dtodx2*(x2Flux[j+1][i].Mx - x2Flux[j][i].Mx);
      pG->U[ks][j][i].M3 -= dtodx2*(x2Flux[j+1][i].My - x2Flux[j][i].My);
#ifndef BAROTROPIC
      pG->U[ks][j][i].E  -= dtodx2*(x2Flux[j+1][i].E  - x2Flux[j][i].E );
#endif /* BAROTROPIC */
#ifdef MHD
      pG->U[ks][j][i].B3c -= dtodx2*(x2Flux[j+1][i].By - x2Flux[j][i].By);
      pG->U[ks][j][i].B1c -= dtodx2*(x2Flux[j+1][i].Bz - x2Flux[j][i].Bz);
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        pG->U[ks][j][i].s[n] -= dtodx2*(x2Flux[j+1][i].s[n] 
                                         - x2Flux[j  ][i].s[n]);
#endif
    }
  }

/*--- Step 12c -----------------------------------------------------------------
 * Update cell-centered variables in pG using 3D x3-Fluxes
 */

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

/*static mesh refinement part goes here*/
#ifdef STATIC_MESH_REFINEMENT
/*--- Step 12e -----------------------------------------------------------------
 * With SMR, store fluxes at boundaries of child and parent grids.
 */

  for (ncg=0; ncg<pG->NCGrid; ncg++) {

/* x1-boundaries of child Grids (interior to THIS Grid) */

    for (dim=0; dim<2; dim++){
      if (pG->CGrid[ncg].myFlx[dim] != NULL) {

        if (dim==0) i = pG->CGrid[ncg].ijks[0];
        if (dim==1) i = pG->CGrid[ncg].ijke[0] + 1;
        jcs = pG->CGrid[ncg].ijks[1];
        jce = pG->CGrid[ncg].ijke[1];

        for (j=jcs, jj=0; j<=jce; j++, jj++){
          pG->CGrid[ncg].myFlx[dim][ks][jj].d  = x1Flux[j][i].d; 
          pG->CGrid[ncg].myFlx[dim][ks][jj].M1 = x1Flux[j][i].Mx; 
          pG->CGrid[ncg].myFlx[dim][ks][jj].M2 = x1Flux[j][i].My;
          pG->CGrid[ncg].myFlx[dim][ks][jj].M3 = x1Flux[j][i].Mz; 
#ifndef BAROTROPIC
          pG->CGrid[ncg].myFlx[dim][ks][jj].E  = x1Flux[j][i].E; 
#endif /* BAROTROPIC */
#ifdef MHD
          pG->CGrid[ncg].myFlx[dim][ks][jj].B1c = 0.0;
          pG->CGrid[ncg].myFlx[dim][ks][jj].B2c = x1Flux[j][i].By; 
          pG->CGrid[ncg].myFlx[dim][ks][jj].B3c = x1Flux[j][i].Bz; 
#endif /* MHD */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++)
            pG->CGrid[ncg].myFlx[dim][ks][jj].s[n]  = x1Flux[j][i].s[n]; 
#endif
        }
#ifdef MHD
        for (j=jcs, jj=0; j<=jce+1; j++, jj++){
          pG->CGrid[ncg].myEMF3[dim][ks][jj] = emf3[j][i];
        }
#endif /* MHD */
      }
    }

/* x2-boundaries of child Grids (interior to THIS Grid) */

    for (dim=2; dim<4; dim++){
      if (pG->CGrid[ncg].myFlx[dim] != NULL) {

        ics = pG->CGrid[ncg].ijks[0];
        ice = pG->CGrid[ncg].ijke[0];
        if (dim==2) j = pG->CGrid[ncg].ijks[1];
        if (dim==3) j = pG->CGrid[ncg].ijke[1] + 1;

        for (i=ics, ii=0; i<=ice; i++, ii++){
          pG->CGrid[ncg].myFlx[dim][ks][ii].d  = x2Flux[j][i].d; 
          pG->CGrid[ncg].myFlx[dim][ks][ii].M1 = x2Flux[j][i].Mz; 
          pG->CGrid[ncg].myFlx[dim][ks][ii].M2 = x2Flux[j][i].Mx;
          pG->CGrid[ncg].myFlx[dim][ks][ii].M3 = x2Flux[j][i].My; 
#ifndef BAROTROPIC
          pG->CGrid[ncg].myFlx[dim][ks][ii].E  = x2Flux[j][i].E; 
#endif /* BAROTROPIC */
#ifdef MHD
          pG->CGrid[ncg].myFlx[dim][ks][ii].B1c = x2Flux[j][i].Bz; 
          pG->CGrid[ncg].myFlx[dim][ks][ii].B2c = 0.0;
          pG->CGrid[ncg].myFlx[dim][ks][ii].B3c = x2Flux[j][i].By; 
#endif /* MHD */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++)
            pG->CGrid[ncg].myFlx[dim][ks][ii].s[n]  = x2Flux[j][i].s[n]; 
#endif
        }
#ifdef MHD
        for (i=ics, ii=0; i<=ice+1; i++, ii++){
          pG->CGrid[ncg].myEMF3[dim][ks][ii] = emf3[j][i];
        }
#endif /* MHD */
      }
    }
  }

  for (npg=0; npg<pG->NPGrid; npg++) {

/* x1-boundaries of parent Grids (at boundaries of THIS Grid)  */

    for (dim=0; dim<2; dim++){
      if (pG->PGrid[npg].myFlx[dim] != NULL) {

        if (dim==0) i = pG->PGrid[npg].ijks[0];
        if (dim==1) i = pG->PGrid[npg].ijke[0] + 1;
        jps = pG->PGrid[npg].ijks[1];
        jpe = pG->PGrid[npg].ijke[1];

        for (j=jps, jj=0; j<=jpe; j++, jj++){
          pG->PGrid[npg].myFlx[dim][ks][jj].d  = x1Flux[j][i].d; 
          pG->PGrid[npg].myFlx[dim][ks][jj].M1 = x1Flux[j][i].Mx; 
          pG->PGrid[npg].myFlx[dim][ks][jj].M2 = x1Flux[j][i].My;
          pG->PGrid[npg].myFlx[dim][ks][jj].M3 = x1Flux[j][i].Mz; 
#ifndef BAROTROPIC
          pG->PGrid[npg].myFlx[dim][ks][jj].E  = x1Flux[j][i].E; 
#endif /* BAROTROPIC */
#ifdef MHD
          pG->PGrid[npg].myFlx[dim][ks][jj].B1c = 0.0;
          pG->PGrid[npg].myFlx[dim][ks][jj].B2c = x1Flux[j][i].By; 
          pG->PGrid[npg].myFlx[dim][ks][jj].B3c = x1Flux[j][i].Bz; 
#endif /* MHD */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++)
            pG->PGrid[npg].myFlx[dim][ks][jj].s[n]  = x1Flux[j][i].s[n]; 
#endif
        }
#ifdef MHD
        for (j=jps, jj=0; j<=jpe+1; j++, jj++){
          pG->PGrid[npg].myEMF3[dim][ks][jj] = emf3[j][i];
        }
#endif /* MHD */
      }
    }

/* x2-boundaries of parent Grids (at boundaries of THIS Grid)  */

    for (dim=2; dim<4; dim++){
      if (pG->PGrid[npg].myFlx[dim] != NULL) {

        ips = pG->PGrid[npg].ijks[0];
        ipe = pG->PGrid[npg].ijke[0];
        if (dim==2) j = pG->PGrid[npg].ijks[1];
        if (dim==3) j = pG->PGrid[npg].ijke[1] + 1;

        for (i=ips, ii=0; i<=ipe; i++, ii++){
          pG->PGrid[npg].myFlx[dim][ks][ii].d  = x2Flux[j][i].d; 
          pG->PGrid[npg].myFlx[dim][ks][ii].M1 = x2Flux[j][i].Mz; 
          pG->PGrid[npg].myFlx[dim][ks][ii].M2 = x2Flux[j][i].Mx;
          pG->PGrid[npg].myFlx[dim][ks][ii].M3 = x2Flux[j][i].My; 
#ifndef BAROTROPIC
          pG->PGrid[npg].myFlx[dim][ks][ii].E  = x2Flux[j][i].E; 
#endif /* BAROTROPIC */
#ifdef MHD
          pG->PGrid[npg].myFlx[dim][ks][ii].B1c = x2Flux[j][i].Bz; 
          pG->PGrid[npg].myFlx[dim][ks][ii].B2c = 0.0;
          pG->PGrid[npg].myFlx[dim][ks][ii].B3c = x2Flux[j][i].By; 
#endif /* MHD */
#if (NSCALARS > 0)
          for (n=0; n<NSCALARS; n++)
            pG->PGrid[npg].myFlx[dim][ks][ii].s[n]  = x2Flux[j][i].s[n]; 
#endif
        }
#ifdef MHD
        for (i=ips, ii=0; i<=ipe+1; i++, ii++){
          pG->PGrid[npg].myEMF3[dim][ks][ii] = emf3[j][i];
        }
#endif /* MHD */
      }
    }
  }

#endif /* STATIC_MESH_REFINEMENT */
}

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_init_2d(MeshS *pM)
 *  \brief Allocate temporary integration arrays 
*/
void integrate_init_1d(MeshS *pM)
{
  int nmax,size1=0,size2=0,size3=0,nl,nd;

/* Cycle over all Grids on this processor to find maximum Nx1, Nx2, Nx3 */
  for (nl=0; nl<(pM->NLevels); nl++){
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++){
      if (pM->Domain[nl][nd].Grid != NULL) {
        if (pM->Domain[nl][nd].Grid->Nx[0] > size1){
          size1 = pM->Domain[nl][nd].Grid->Nx[0];
        }
        if (pM->Domain[nl][nd].Grid->Nx[1] > size2){
          size2 = pM->Domain[nl][nd].Grid->Nx[1];
        }
        if (pM->Domain[nl][nd].Grid->Nx[2] > size3){
          size3 = pM->Domain[nl][nd].Grid->Nx[2];
        }
      }
    }
  }

  size1 = size1 + 2*nghost;
  size2 = size2 + 2*nghost;
  size3 = size3 + 2*nghost;
  nmax = MAX((MAX(size1,size2)),size3);

/*refer to material  integrate_2d_ctu.c*/
  if ((Bxc = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;
  if ((Bxi = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;


  if ((U1d= (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((W  = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;

  if ((x1Flux   =(Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS)))==NULL)
    goto on_error;
  if ((x2Flux   =(Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS)))==NULL)
    goto on_error;



  return;

  on_error:
    integrate_destruct();
    ath_error("[integrate_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_destruct_2d(void)
 *  \brief Free temporary integration arrays 
 */
void integrate_destruct_1d(void)
{
/*refer to material  integrate_3d_ctu.c*/

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


#endif /* SMAUG_INTEGRATOR */
