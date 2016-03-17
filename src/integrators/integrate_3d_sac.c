#include "../copyright.h"
/*============================================================================*/
/*! \file integrate_3d_sac.c
 *  \Compute fluxes using sac . 
 *
 * PURPOSE: Integrate MHD equations using 3D version of the directionally
 *   .  The variables updated are:
 *   -  U.[d,M1,M2,M3,E,B1c,B2c,B3c,s] -- where U is of type ConsS
 *   -  B1i, B2i, B3i  -- interface magnetic field
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
 * - integrate_3d_sac()
 * - integrate_init_3d()
 * - integrate_destruct_3d() */

/*
 * PROGRESS
 * Initial boiler plate based on athena integrate_3d_ctu.c

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

#if defined(CTU_INTEGRATOR)
#ifdef SPECIAL_RELATIVITY
#error : The CTU integrator cannot be used for special relativity.
#endif /* SPECIAL_RELATIVITY */

/* The L/R states of conserved variables and fluxes at each cell face */
static Cons1DS ***Ul_x1Face=NULL, ***Ur_x1Face=NULL;
static Cons1DS ***Ul_x2Face=NULL, ***Ur_x2Face=NULL;
static Cons1DS ***Ul_x3Face=NULL, ***Ur_x3Face=NULL;
Cons1DS ***x1Flux=NULL, ***x2Flux=NULL, ***x3Flux=NULL;

/* The interface magnetic fields and emfs */
#ifdef MHD
static Real ***B1_x1Face=NULL, ***B2_x2Face=NULL, ***B3_x3Face=NULL;
Real ***emf1=NULL, ***emf2=NULL, ***emf3=NULL;
static Real ***emf1_cc=NULL, ***emf2_cc=NULL, ***emf3_cc=NULL;
#endif /* MHD */

/* 1D scratch vectors used by lr_states and flux functions */
static Real *Bxc=NULL, *Bxi=NULL;
static Prim1DS *W=NULL, *Wl=NULL, *Wr=NULL;
static Cons1DS *U1d=NULL;

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

/* variables need for cylindrical coordinates */
#ifdef CYLINDRICAL
static Real ***geom_src=NULL;
#endif

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES: 
 *   integrate_emf1_corner() - the upwind CT method in GS05, for emf1
 *   integrate_emf2_corner() - the upwind CT method in GS05, for emf2
 *   integrate_emf3_corner() - the upwind CT method in GS05, for emf3
 *============================================================================*/

#ifdef MHD
static void integrate_emf1_corner(const GridS *pG);
static void integrate_emf2_corner(const GridS *pG);
static void integrate_emf3_corner(const GridS *pG);
#endif /* MHD */

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void integrate_3d_ctu(DomainS *pD)
 *  \brief 3D CTU integrator for MHD using 6-solve method */

void integrate_3d_sac(DomainS *pD)
{
  GridS *pG=(pD->Grid);
  Real dtodx1=pG->dt/pG->dx1, dtodx2=pG->dt/pG->dx2, dtodx3=pG->dt/pG->dx3;
  Real hdt = 0.5*pG->dt, dx2=pG->dx2;
  Real q1 = 0.5*dtodx1, q2 = 0.5*dtodx2, q3 = 0.5*dtodx3;
  int i,il,iu, is = pG->is, ie = pG->ie;
  int j,jl,ju, js = pG->js, je = pG->je;
  int k,kl,ku, ks = pG->ks, ke = pG->ke;
  Real x1,x2,x3,phicl,phicr,phifc,phil,phir,phic,M1h,M2h,M3h,Bx=0.0;
#ifndef BAROTROPIC
  Real coolfl,coolfr,coolf,Eh=0.0;
#endif
#ifdef MHD
  Real MHD_src_By,MHD_src_Bz,mdb1,mdb2,mdb3;
  Real db1,db2,db3,l1,l2,l3,B1,B2,B3,V1,V2,V3;
  Real B1ch,B2ch,B3ch;
#endif
// #if defined(MHD) || defined(SELF_GRAVITY)
  Real dx1i=1.0/pG->dx1, dx2i=1.0/pG->dx2, dx3i=1.0/pG->dx3;
// #endif
#ifdef H_CORRECTION
  Real cfr,cfl,lambdar,lambdal;
#endif
#if (NSCALARS > 0)
  int n;
#endif
#ifdef SELF_GRAVITY
  Real gxl,gxr,gyl,gyr,gzl,gzr,flx_m1l,flx_m1r,flx_m2l,flx_m2r,flx_m3l,flx_m3r;
#endif
#ifdef FEEDBACK
  Real dt1 = 1.0/pG->dt;
#endif
#ifdef SHEARING_BOX
  int my_iproc,my_jproc,my_kproc;
  Real M1n, dM2n; /* M1, dM2=(My+d*q*Omega_0*x) at time n */
  Real M1e, dM2e; /* M1, dM2 evolved by dt/2 */
  Real flx1_dM2, frx1_dM2, flx2_dM2, frx2_dM2, flx3_dM2, frx3_dM2;
  Real fact, qom, om_dt = Omega_0*pG->dt;
#endif /* SHEARING_BOX */

#ifdef ROTATING_FRAME
#ifdef FARGO
#error: Fargo cannot be used in rotating frame.
#endif
#ifndef CYLINDRICAL
#error: ROTATING_FRAME has to be in CYLINDRICAL coordinates.
#endif
  Real tmp_M1, tmp_M2;
#endif /* ROTATING_FRAME */

#ifdef STATIC_MESH_REFINEMENT
  int ncg,npg,dim;
  int ii,ics,ice,jj,jcs,jce,kk,kcs,kce,ips,ipe,jps,jpe,kps,kpe;
#endif

  /* cylindrical variables */
#ifdef CYLINDRICAL
#ifndef ISOTHERMAL
  Real Pavgh;
#endif
  Real rinv, geom_src_d, geom_src_Vx, geom_src_Vy, geom_src_P, geom_src_By, geom_src_Bz;
  const Real *r=pG->r, *ri=pG->ri;
#ifdef FARGO
  Real Om, qshear, Mrn, Mpn, Mre, Mpe, Mrav, Mpav;
#endif
#endif /* CYLINDRICAL */
  Real g,gl,gr;
  Real lsf=1.0, rsf=1.0;

#if defined(CYLINDRICAL) && defined(FARGO)
  if (OrbitalProfile==NULL || ShearProfile==NULL)
    ath_error("[integrate_3d_ctu]:  OrbitalProfile() and ShearProfile() *must* be defined.\n");
#endif

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

  for (k=kl; k<=ku; k++) {
    for (j=jl; j<=ju; j++) {
      for (i=is-nghost; i<=ie+nghost; i++) {
        U1d[i].d  = pG->U[k][j][i].d;
        U1d[i].Mx = pG->U[k][j][i].M1;
        U1d[i].My = pG->U[k][j][i].M2;
        U1d[i].Mz = pG->U[k][j][i].M3;
#ifndef BAROTROPIC
        U1d[i].E  = pG->U[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
        U1d[i].By = pG->U[k][j][i].B2c;
        U1d[i].Bz = pG->U[k][j][i].B3c;
        Bxc[i] = pG->U[k][j][i].B1c;
        Bxi[i] = pG->B1i[k][j][i];
        B1_x1Face[k][j][i] = pG->B1i[k][j][i];
#endif /* MHD */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) U1d[i].s[n] = pG->U[k][j][i].s[n];
#endif
      }

/*--- Step 1b ------------------------------------------------------------------
 * Compute L and R states at X1-interfaces, add "MHD source terms" for 0.5*dt
 */


/*--- Step 1c ------------------------------------------------------------------
 * Add source terms from static gravitational potential for 0.5*dt to L/R states
 */

      if (StaticGravPot != NULL){
        for (i=il+1; i<=iu; i++) {
          cc_pos(pG,i,j,k,&x1,&x2,&x3);

          phicr = (*StaticGravPot)( x1             ,x2,x3);
          phicl = (*StaticGravPot)((x1-    pG->dx1),x2,x3);
          phifc = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

          gl = 2.0*(phifc - phicl)*dx1i;
          gr = 2.0*(phicr - phifc)*dx1i;
#if defined(CYLINDRICAL) && defined(FARGO)
          gl -= r[i-1]*SQR((*OrbitalProfile)(r[i-1]));
          gr -= r[i  ]*SQR((*OrbitalProfile)(r[i  ]));
#endif

          Wl[i].Vx -= hdt*gl;
          Wr[i].Vx -= hdt*gr;
        }
      }



/*--- Step 1c (cont) -----------------------------------------------------------
 * Add the geometric source-terms now using cell-centered primitive
 * variables at time t^n
 */



/*--- Step 1d ------------------------------------------------------------------
 * Compute 1D fluxes in x1-direction, storing into 3D array
 */


/*=== STEP 2: Compute L/R x2-interface states and 1D x2-Fluxes ===============*/

/*--- Step 2a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M2, M3, M1, E, B3c, B1c, s[n])
 */



/*--- Step 2b ------------------------------------------------------------------
 * Compute L and R states at X2-interfaces, add "MHD source terms" for 0.5*dt
 */



/*--- Step 2c ------------------------------------------------------------------
 * Add source terms from static gravitational potential for 0.5*dt to L/R states
 */

      if (StaticGravPot != NULL){
        for (j=jl+1; j<=ju; j++) {
          cc_pos(pG,i,j,k,&x1,&x2,&x3);
          phicr = (*StaticGravPot)(x1, x2             ,x3);
          phicl = (*StaticGravPot)(x1,(x2-    pG->dx2),x3);
          phifc = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

          Wl[j].Vx -= dtodx2*(phifc - phicl);
          Wr[j].Vx -= dtodx2*(phicr - phifc);
        }
      }



/*--- Step 2d ------------------------------------------------------------------
 * Compute 1D fluxes in x2-direction, storing into 3D array
 */



/*--- Step 3c ------------------------------------------------------------------
 * Add source terms from static gravitational potential for 0.5*dt to L/R states
 */

      if (StaticGravPot != NULL){
        for (k=kl+1; k<=ku; k++) {
          cc_pos(pG,i,j,k,&x1,&x2,&x3);
          phicr = (*StaticGravPot)(x1,x2, x3             );
          phicl = (*StaticGravPot)(x1,x2,(x3-    pG->dx3));
          phifc = (*StaticGravPot)(x1,x2,(x3-0.5*pG->dx3));

          Wl[k].Vx -= dtodx3*(phifc - phicl);
          Wr[k].Vx -= dtodx3*(phicr - phifc);
        }
      }



/*--- Step 3d ------------------------------------------------------------------
 * Compute 1D fluxes in x3-direction, storing into 3D array
 */

      for (k=kl+1; k<=ku; k++) {
        Ul_x3Face[k][j][i] = Prim1D_to_Cons1D(&Wl[k],&Bxi[k]);
        Ur_x3Face[k][j][i] = Prim1D_to_Cons1D(&Wr[k],&Bxi[k]);

#ifdef MHD
        Bx = B3_x3Face[k][j][i];
#endif
        fluxes(Ul_x3Face[k][j][i],Ur_x3Face[k][j][i],Wl[k],Wr[k],Bx,
          &x3Flux[k][j][i]);
      }
    }
  }

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

  if (StaticGravPot != NULL){
  for (k=kl+1; k<=ku-1; k++) {
    for (j=jl+1; j<=ju-1; j++) {
      for (i=il+1; i<=iu; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        phic = (*StaticGravPot)(x1, x2             ,x3);
        phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

        /* correct right states; x2 and x3 gradients */
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

        /* correct left states; x2 and x3 gradients */
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
  }}




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

  if (StaticGravPot != NULL){
  for (k=kl+1; k<=ku-1; k++) {
    for (j=jl+1; j<=ju; j++) {
      for (i=il+1; i<=iu-1; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);

        /* correct right states; x1 and x3 gradients */
        phic = (*StaticGravPot)((x1            ),x2,x3);
        phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

#ifdef CYLINDRICAL
        rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
#endif
        g = (phir-phil)/pG->dx1;
#if defined(CYLINDRICAL) && defined(FARGO)
        g -= r[i]*SQR((*OrbitalProfile)(r[i])); 
#endif
        Ur_x2Face[k][j][i].Mz -= hdt*pG->U[k][j][i].d*g;
#ifndef BAROTROPIC
        Ur_x2Face[k][j][i].E -= q1*(lsf*x1Flux[k][j  ][i  ].d*(phic - phil)
                                  + rsf*x1Flux[k][j  ][i+1].d*(phir - phic));
        #ifdef ROTATING_FRAME
                Ur_x2Face[k][j][i].E -= hdt * 0.5*(x1Flux[k][j  ][i  ].d+x1Flux[k][j  ][i+1].d)
                                        * SQR(Omega_0)*Rc*cos(x2);
	#endif
#endif

        phir = (*StaticGravPot)(x1,x2,(x3+0.5*pG->dx3));
        phil = (*StaticGravPot)(x1,x2,(x3-0.5*pG->dx3));

        Ur_x2Face[k][j][i].My -= q3*(phir-phil)*pG->U[k][j][i].d;
#ifndef BAROTROPIC
        Ur_x2Face[k][j][i].E -= q3*(x3Flux[k  ][j  ][i].d*(phic - phil)
                                  + x3Flux[k+1][j  ][i].d*(phir - phic));
#endif

        /* correct left states; x1 and x3 gradients */
        phic = (*StaticGravPot)((x1            ),(x2-pG->dx2),x3);
        phir = (*StaticGravPot)((x1+0.5*pG->dx1),(x2-pG->dx2),x3);
        phil = (*StaticGravPot)((x1-0.5*pG->dx1),(x2-pG->dx2),x3);

        g = (phir-phil)/pG->dx1;
#if defined(CYLINDRICAL) && defined(FARGO)
        g -= r[i]*SQR((*OrbitalProfile)(r[i])); 
#endif
        Ul_x2Face[k][j][i].Mz -= hdt*pG->U[k][j-1][i].d*g;
#ifndef BAROTROPIC
        Ul_x2Face[k][j][i].E -= q1*(lsf*x1Flux[k][j-1][i  ].d*(phic - phil)
                                  + rsf*x1Flux[k][j-1][i+1].d*(phir - phic));
        #ifdef ROTATING_FRAME
                Ul_x2Face[k][j][i].E -= hdt * 0.5*(x1Flux[k][j-1][i  ].d+x1Flux[k][j-1][i+1].d)
                                        * SQR(Omega_0)*Rc*cos(x2-pG->dx2);
	#endif
#endif
        phir = (*StaticGravPot)(x1,(x2-pG->dx2),(x3+0.5*pG->dx3));
        phil = (*StaticGravPot)(x1,(x2-pG->dx2),(x3-0.5*pG->dx3));

        Ul_x2Face[k][j][i].My -= q3*(phir-phil)*pG->U[k][j-1][i].d;
#ifndef BAROTROPIC
        Ul_x2Face[k][j][i].E -= q3*(x3Flux[k  ][j-1][i].d*(phic - phil)
                                  + x3Flux[k+1][j-1][i].d*(phir - phic));
#endif
      }
    }
  }}





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

  if (StaticGravPot != NULL){
  for (k=kl+1; k<=ku; k++) {
    for (j=jl+1; j<=ju-1; j++) {
      for (i=il+1; i<=iu-1; i++) {
        cc_pos(pG,i,j,k,&x1,&x2,&x3);

        /* correct right states; x1 and x2 gradients */
        phic = (*StaticGravPot)((x1            ),x2,x3);
        phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,x3);
        phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);

#ifdef CYLINDRICAL
        rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
#endif
        g = (phir-phil)/pG->dx1;
#if defined(CYLINDRICAL) && defined(FARGO)
        g -= r[i]*SQR((*OrbitalProfile)(r[i])); 
#endif
        Ur_x3Face[k][j][i].My -= hdt*pG->U[k][j][i].d*g;
#ifndef BAROTROPIC
        Ur_x3Face[k][j][i].E -= q1*(lsf*x1Flux[k  ][j][i  ].d*(phic - phil)
                                  + rsf*x1Flux[k  ][j][i+1].d*(phir - phic));
        #ifdef ROTATING_FRAME
                Ur_x3Face[k][j][i].E -= hdt * 0.5*(x1Flux[k][j][i  ].d+x1Flux[k][j][i+1].d)
                                        * SQR(Omega_0)*Rc*cos(x2);
	#endif
#endif

        phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),x3);
        phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

        Ur_x3Face[k][j][i].Mz -= q2*(phir-phil)*pG->U[k][j][i].d;
#ifndef BAROTROPIC
        Ur_x3Face[k][j][i].E -= q2*(x2Flux[k  ][j  ][i].d*(phic - phil)
                                  + x2Flux[k  ][j+1][i].d*(phir - phic));
        #ifdef ROTATING_FRAME
                Ur_x3Face[k][j][i].E += hdt * 0.5*(x2Flux[k][j  ][i].d*sin(x2-0.5*pG->dx2) +
                                                   x2Flux[k][j+1][i].d*sin(x2+0.5*pG->dx2)) *SQR(Omega_0)*Rc;
        #endif
#endif

        /* correct left states; x1 and x2 gradients */
        phic = (*StaticGravPot)((x1            ),x2,(x3-pG->dx3));
        phir = (*StaticGravPot)((x1+0.5*pG->dx1),x2,(x3-pG->dx3));
        phil = (*StaticGravPot)((x1-0.5*pG->dx1),x2,(x3-pG->dx3));

        g = (phir-phil)/pG->dx1;
#if defined(CYLINDRICAL) && defined(FARGO)
        g -= r[i]*SQR((*OrbitalProfile)(r[i])); 
#endif
        Ul_x3Face[k][j][i].My -= hdt*pG->U[k-1][j][i].d*g;
#ifndef BAROTROPIC
        Ul_x3Face[k][j][i].E -= q1*(lsf*x1Flux[k-1][j][i  ].d*(phic - phil)
                                  + rsf*x1Flux[k-1][j][i+1].d*(phir - phic));
        #ifdef ROTATING_FRAME
                Ul_x3Face[k][j][i].E -= hdt * 0.5*(x1Flux[k-1][j][i  ].d+x1Flux[k-1][j][i+1].d)
                                         * SQR(Omega_0)*Rc*cos(x2);
	#endif
#endif

        phir = (*StaticGravPot)(x1,(x2+0.5*pG->dx2),(x3-pG->dx3));
        phil = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),(x3-pG->dx3));

        Ul_x3Face[k][j][i].Mz -= q2*(phir-phil)*pG->U[k-1][j][i].d;
#ifndef BAROTROPIC
        Ul_x3Face[k][j][i].E -= q2*(x2Flux[k-1][j  ][i].d*(phic - phil)
                                  + x2Flux[k-1][j+1][i].d*(phir - phic));
        #ifdef ROTATING_FRAME
                Ul_x3Face[k][j][i].E += hdt * 0.5*(x2Flux[k-1][j  ][i].d*sin(x2-0.5*pG->dx2) +
                                                x2Flux[k-1][j+1][i].d*sin(x2+0.5*pG->dx2)) *SQR(Omega_0)*Rc;
        #endif
#endif
      }
    }
  }}


/*--- Step 7e ------------------------------------------------------------------
 * Apply density floor
 */


/*=== STEP 8: Compute cell-centered values at n+1/2 ==========================*/


/*=== STEP 9: Compute 3D x1-Flux, x2-Flux, x3-Flux ===========================*/

/*--- Step 9a ------------------------------------------------------------------
 * Compute maximum wavespeeds in multidimensions (eta in eq. 10 from Sanders et
 *  al. (1998)) for H-correction
 */

#ifdef H_CORRECTION
  for (k=ks-1; k<=ke+1; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=is-1; i<=ie+2; i++) {
#ifdef MHD
        Bx = B1_x1Face[k][j][i];
#endif
        cfr = cfast(&(Ur_x1Face[k][j][i]),&Bx);
        cfl = cfast(&(Ul_x1Face[k][j][i]),&Bx);
        lambdar = Ur_x1Face[k][j][i].Mx/Ur_x1Face[k][j][i].d + cfr;
        lambdal = Ul_x1Face[k][j][i].Mx/Ul_x1Face[k][j][i].d - cfl;
        eta1[k][j][i] = 0.5*fabs(lambdar - lambdal);
      }
    }
  }

  for (k=ks-1; k<=ke+1; k++) {
    for (j=js-1; j<=je+2; j++) {
      for (i=is-1; i<=ie+1; i++) {
#ifdef MHD
        Bx = B2_x2Face[k][j][i];
#endif
        cfr = cfast(&(Ur_x2Face[k][j][i]),&Bx);
        cfl = cfast(&(Ul_x2Face[k][j][i]),&Bx);
        lambdar = Ur_x2Face[k][j][i].Mx/Ur_x2Face[k][j][i].d + cfr;
        lambdal = Ul_x2Face[k][j][i].Mx/Ul_x2Face[k][j][i].d - cfl;
        eta2[k][j][i] = 0.5*fabs(lambdar - lambdal);
      }
    }
  }

  for (k=ks-1; k<=ke+2; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=is-1; i<=ie+1; i++) {
#ifdef MHD
        Bx = B3_x3Face[k][j][i];
#endif
        cfr = cfast(&(Ur_x3Face[k][j][i]),&Bx);
        cfl = cfast(&(Ul_x3Face[k][j][i]),&Bx);
        lambdar = Ur_x3Face[k][j][i].Mx/Ur_x3Face[k][j][i].d + cfr;
        lambdal = Ul_x3Face[k][j][i].Mx/Ul_x3Face[k][j][i].d - cfl;
        eta3[k][j][i] = 0.5*fabs(lambdar - lambdal);
      }
    }
  }
#endif /* H_CORRECTION */

/*--- Step 9b ------------------------------------------------------------------
 * Compute 3D x1-fluxes from corrected L/R states.
 */

  for (k=ks-1; k<=ke+1; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=is; i<=ie+1; i++) {
#ifdef H_CORRECTION
        etah = MAX(eta2[k][j][i-1],eta2[k][j][i]);
        etah = MAX(etah,eta2[k][j+1][i-1]);
        etah = MAX(etah,eta2[k][j+1][i  ]);

        etah = MAX(etah,eta3[k  ][j][i-1]);
        etah = MAX(etah,eta3[k  ][j][i  ]);
        etah = MAX(etah,eta3[k+1][j][i-1]);
        etah = MAX(etah,eta3[k+1][j][i  ]);

        etah = MAX(etah,eta1[k  ][j][i  ]);
#endif /* H_CORRECTION */
#ifdef MHD
        Bx = B1_x1Face[k][j][i];
#endif
        Wl[i] = Cons1D_to_Prim1D(&Ul_x1Face[k][j][i],&Bx);
        Wr[i] = Cons1D_to_Prim1D(&Ur_x1Face[k][j][i],&Bx);

        fluxes(Ul_x1Face[k][j][i],Ur_x1Face[k][j][i],Wl[i],Wr[i],Bx,
               &x1Flux[k][j][i]);
      }
    }
  }

/*--- Step 9c ------------------------------------------------------------------
 * Compute 3D x2-fluxes from corrected L/R states.
 */

  for (k=ks-1; k<=ke+1; k++) {
    for (j=js; j<=je+1; j++) {
      for (i=is-1; i<=ie+1; i++) {
#ifdef H_CORRECTION
        etah = MAX(eta1[k][j-1][i],eta1[k][j][i]);
        etah = MAX(etah,eta1[k][j-1][i+1]);
        etah = MAX(etah,eta1[k][j  ][i+1]);

        etah = MAX(etah,eta3[k  ][j-1][i]);
        etah = MAX(etah,eta3[k  ][j  ][i]);
        etah = MAX(etah,eta3[k+1][j-1][i]);
        etah = MAX(etah,eta3[k+1][j  ][i]);

        etah = MAX(etah,eta2[k  ][j  ][i]);
#endif /* H_CORRECTION */
#ifdef MHD
        Bx = B2_x2Face[k][j][i];
#endif
        Wl[i] = Cons1D_to_Prim1D(&Ul_x2Face[k][j][i],&Bx);
        Wr[i] = Cons1D_to_Prim1D(&Ur_x2Face[k][j][i],&Bx);

        fluxes(Ul_x2Face[k][j][i],Ur_x2Face[k][j][i],Wl[i],Wr[i],Bx,
               &x2Flux[k][j][i]);
      }
    }
  }

/*--- Step 9d ------------------------------------------------------------------
 * Compute 3D x3-fluxes from corrected L/R states.
 */

  for (k=ks; k<=ke+1; k++) {
    for (j=js-1; j<=je+1; j++) {
      for (i=is-1; i<=ie+1; i++) {
#ifdef H_CORRECTION
        etah = MAX(eta1[k-1][j][i],eta1[k][j][i]);
        etah = MAX(etah,eta1[k-1][j][i+1]);
        etah = MAX(etah,eta1[k][j  ][i+1]);

        etah = MAX(etah,eta2[k-1][j  ][i]);
        etah = MAX(etah,eta2[k  ][j  ][i]);
        etah = MAX(etah,eta2[k-1][j+1][i]);
        etah = MAX(etah,eta2[k  ][j+1][i]);

        etah = MAX(etah,eta3[k  ][j  ][i]);
#endif /* H_CORRECTION */
#ifdef MHD
        Bx = B3_x3Face[k][j][i];
#endif
        Wl[i] = Cons1D_to_Prim1D(&Ul_x3Face[k][j][i],&Bx);
        Wr[i] = Cons1D_to_Prim1D(&Ur_x3Face[k][j][i],&Bx);

        fluxes(Ul_x3Face[k][j][i],Ur_x3Face[k][j][i],Wl[i],Wr[i],Bx,
               &x3Flux[k][j][i]);
      }
    }
  }

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

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
        rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
#endif
        pG->U[k][j][i].d  -= dtodx1*(rsf*x1Flux[k][j][i+1].d -lsf*x1Flux[k][j][i].d );
        pG->U[k][j][i].M1 -= dtodx1*(rsf*x1Flux[k][j][i+1].Mx-lsf*x1Flux[k][j][i].Mx);
        pG->U[k][j][i].M2 -= dtodx1*(SQR(rsf)*x1Flux[k][j][i+1].My-SQR(lsf)*x1Flux[k][j][i].My);
        pG->U[k][j][i].M3 -= dtodx1*(rsf*x1Flux[k][j][i+1].Mz-lsf*x1Flux[k][j][i].Mz);
#ifndef BAROTROPIC
        pG->U[k][j][i].E  -= dtodx1*(rsf*x1Flux[k][j][i+1].E -lsf*x1Flux[k][j][i].E );
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++)
          pG->U[k][j][i].s[n] -= dtodx1*(rsf*x1Flux[k][j][i+1].s[n]
                                       - lsf*x1Flux[k][j][i  ].s[n]);
#endif
      }
    }
  }

/*--- Step 12b -----------------------------------------------------------------
 * Update cell-centered variables in pG using 3D x2-Fluxes
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
        dtodx2 = pG->dt/(r[i]*pG->dx2);
#endif
        pG->U[k][j][i].d  -= dtodx2*(x2Flux[k][j+1][i].d -x2Flux[k][j][i].d );
        pG->U[k][j][i].M1 -= dtodx2*(x2Flux[k][j+1][i].Mz-x2Flux[k][j][i].Mz);
        pG->U[k][j][i].M2 -= dtodx2*(x2Flux[k][j+1][i].Mx-x2Flux[k][j][i].Mx);
        pG->U[k][j][i].M3 -= dtodx2*(x2Flux[k][j+1][i].My-x2Flux[k][j][i].My);
#ifndef BAROTROPIC
        pG->U[k][j][i].E -=dtodx2*(x2Flux[k][j+1][i].E -x2Flux[k][j][i].E );
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++)
          pG->U[k][j][i].s[n] -= dtodx2*(x2Flux[k][j+1][i].s[n]
                                       - x2Flux[k][j  ][i].s[n]);
#endif
      }
    }
  }

/*--- Step 12c -----------------------------------------------------------------
 * Update cell-centered variables in pG using 3D x3-Fluxes
 */

  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pG->U[k][j][i].d  -= dtodx3*(x3Flux[k+1][j][i].d -x3Flux[k][j][i].d );
        pG->U[k][j][i].M1 -= dtodx3*(x3Flux[k+1][j][i].My-x3Flux[k][j][i].My);
        pG->U[k][j][i].M2 -= dtodx3*(x3Flux[k+1][j][i].Mz-x3Flux[k][j][i].Mz);
        pG->U[k][j][i].M3 -= dtodx3*(x3Flux[k+1][j][i].Mx-x3Flux[k][j][i].Mx);
#ifndef BAROTROPIC
        pG->U[k][j][i].E  -= dtodx3*(x3Flux[k+1][j][i].E -x3Flux[k][j][i].E );
#endif /* BAROTROPIC */
#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++)
          pG->U[k][j][i].s[n] -= dtodx3*(x3Flux[k+1][j][i].s[n]
                                       - x3Flux[k  ][j][i].s[n]);
#endif
      }
    }
  }

/*--- Step 12d -----------------------------------------------------------------
 * Set cell centered magnetic fields to average of updated face centered fields.
 */

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
#ifdef CYLINDRICAL
        rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
#endif
        pG->U[k][j][i].B1c = 0.5*(lsf*pG->B1i[k][j][i] + rsf*pG->B1i[k][j][i+1]);
        pG->U[k][j][i].B2c = 0.5*(    pG->B2i[k][j][i] +     pG->B2i[k][j+1][i]);
        pG->U[k][j][i].B3c = 0.5*(    pG->B3i[k][j][i] +     pG->B3i[k+1][j][i]);
      }
    }
  }
#endif /* MHD */

/*static mesh refinement part goes here*/


  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_init_3d(MeshS *pM)
 *  \brief Allocate temporary integration arrays 
*/
void integrate_init_3d(MeshS *pM)
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

/*refer to material  integrate_3d_ctu.c*/


  return;

  on_error:
    integrate_destruct();
    ath_error("[integrate_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_destruct_3d(void)
 *  \brief Free temporary integration arrays 
 */
void integrate_destruct_3d(void)
{
/*refer to material  integrate_3d_ctu.c*/

  return;
}

/*=========================== PRIVATE FUNCTIONS ==============================*/

/*----------------------------------------------------------------------------*/
/*! \fn static void integrate_emf1_corner(const GridS *pG)
 *  \brief Integrates face centered B-fluxes to compute corner EMFs.  
 *
 *  Note: 
 * - x1Flux.By = VxBy - BxVy = v1*b2-b1*v2 = -EMFZ
 * - x1Flux.Bz = VxBz - BxVz = v1*b3-b1*v3 = EMFY
 * - x2Flux.By = VxBy - BxVy = v2*b3-b2*v3 = -EMFX
 * - x2Flux.Bz = VxBz - BxVz = v2*b1-b2*v1 = EMFZ
 * - x3Flux.By = VxBy - BxVy = v3*b1-b3*v1 = -EMFY
 * - x3Flux.Bz = VxBz - BxVz = v3*b2-b3*v2 = EMFX 
 */
#ifdef MHD
static void integrate_emf1_corner(const GridS *pG)
{
  
/*refer to material  integrate_3d_ctu.c*/



  return;
}

/*! \fn static void integrate_emf2_corner(const GridS *pG)
 *  \brief Integrates face centered B-fluxes to compute corner EMFs.  
 *
 *  Note: 
 * - x1Flux.By = VxBy - BxVy = v1*b2-b1*v2 = -EMFZ
 * - x1Flux.Bz = VxBz - BxVz = v1*b3-b1*v3 = EMFY
 * - x2Flux.By = VxBy - BxVy = v2*b3-b2*v3 = -EMFX
 * - x2Flux.Bz = VxBz - BxVz = v2*b1-b2*v1 = EMFZ
 * - x3Flux.By = VxBy - BxVy = v3*b1-b3*v1 = -EMFY
 * - x3Flux.Bz = VxBz - BxVz = v3*b2-b3*v2 = EMFX 
 */
static void integrate_emf2_corner(const GridS *pG)
{

/*refer to material  integrate_3d_ctu.c*/

  return;
}

/*! \fn static void integrate_emf3_corner(const GridS *pG)
 *  \brief Integrates face centered B-fluxes to compute corner EMFs.  
 *
 *  Note: 
 * - x1Flux.By = VxBy - BxVy = v1*b2-b1*v2 = -EMFZ
 * - x1Flux.Bz = VxBz - BxVz = v1*b3-b1*v3 = EMFY
 * - x2Flux.By = VxBy - BxVy = v2*b3-b2*v3 = -EMFX
 * - x2Flux.Bz = VxBz - BxVz = v2*b1-b2*v1 = EMFZ
 * - x3Flux.By = VxBy - BxVy = v3*b1-b3*v1 = -EMFY
 * - x3Flux.Bz = VxBz - BxVz = v3*b2-b3*v2 = EMFX 
 */
static void integrate_emf3_corner(const GridS *pG)
{
 
/*refer to material  integrate_3d_ctu.c*/

  return;
}
#endif /* MHD */

#endif /* CTU_INTEGRATOR */
