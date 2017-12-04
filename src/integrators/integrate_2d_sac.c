#include "../copyright.h"
/*============================================================================*/
/*! \file integrate_2d_sac.c
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
 * - integrate_2d_sac()
 * - integrate_init_2d()
 * - integrate_destruct_2d() */

/*
 * PROGRESS
 * Initial boiler plate based on athena integrate_2d_ctu.c

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
static Cons1DS **Uc_x1=NULL;
static Cons1DS **Uc_x2=NULL;
static Cons1DS **x1Flux=NULL, **x2Flux=NULL;


/* The interface magnetic fields and emfs */
/*not needed for sac*/
/*#ifdef MHD
static Real **B1_x1Face=NULL, **B2_x2Face=NULL;
static Real **emf3=NULL, **emf3_cc=NULL;
#endif*/ /* MHD */





/* 1D scratch vectors used by lr_states and flux functions */
static Real *Bxc=NULL,  *Bxb=NULL;
static Prim1DS *W=NULL;/*, not used for sac *Wl=NULL, *Wr=NULL;*/
static Cons1DS *U1d=NULL;/*, not used for sac *Ul=NULL, *Ur=NULL;*/



static Cons1DS **Uc=NULL;


/* The interface magnetic fields and emfs */
#ifdef MHD
static Real **B1_x1=NULL, **B2_x2=NULL;

#endif /* MHD */

/* density and Pressure at t^{n+1/2} needed by MHD, cooling, and gravity */
static Real **dhalf = NULL, **phalf=NULL;

/* variables needed for H-correction of Sanders et al (1998) */
extern Real etah;
#ifdef H_CORRECTION
static Real **eta1=NULL, **eta2=NULL;
#endif

/* variables needed for cylindrical coordinates */
#ifdef CYLINDRICAL
static Real **geom_src=NULL;
#endif



/*following not needed for sac*/

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
/*! \fn void integrate_2d_sac(DomainS *pD)
 *  \brief CTU integrator in 2D.
 *
 *   The numbering of steps follows the numbering in the 3D version.
 *   NOT ALL STEPS ARE NEEDED IN 2D.
 */



void integrate_2d_sac(DomainS *pD)
{
  
  GridS *pG=(pD->Grid);
  Real dtodx1 = pG->dt/pG->dx1, dtodx2 = pG->dt/pG->dx2;
  Real hdtodx1 = 0.5*dtodx1, hdtodx2 = 0.5*dtodx2;
  Real hdt = 0.5*pG->dt, dx2 = pG->dx2;
  int i,il,iu,is=pG->is, ie=pG->ie;
  int j,jl,ju,js=pG->js, je=pG->je;
  int ks=pG->ks;
  Real x1,x2,x3,phicl,phicr,phifc,phil,phir,phic,M1h,M2h,M3h,Bx=0.0;
#ifndef BAROTROPIC
  Real coolfl,coolfr,coolf,coolfc,Eh=0.0;
#endif
#ifdef MHD
  Real MHD_src,dbx,dby,B1,B2,B3,V3;
  Real B1ch, B2ch, B3ch;
#endif
  Real dx1i=1.0/pG->dx1, dx2i=1.0/pG->dx2;
#ifdef H_CORRECTION
  Real cfr,cfl,lambdar,lambdal;
#endif
#if (NSCALARS > 0)
  int n;
#endif
#ifdef SELF_GRAVITY
  Real gxl,gxr,gyl,gyr,flux_m1l,flux_m1r,flux_m2l,flux_m2r;
#endif
#ifdef FEEDBACK
  Real dt1 = 1.0/pG->dt;
#endif
#ifdef SHEARING_BOX
/* in XY shearing sheet 2=phi; in XZ shearing sheet 2=Z and 3=phi */
  Real Vphi, Mphi;
  Real M1n, dM2n, dM3n;   /* M1, dM2/3=(My+d*q*Omega_0*x) at time n */
  Real M1e, dM2e, dM3e;   /* M1, dM2/3 evolved by dt/2  */
  Real flx1_dM2, frx1_dM2, flx2_dM2, frx2_dM2;
  Real flx1_dM3, frx1_dM3, flx2_dM3, frx2_dM3;
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
  int ii,ics,ice,jj,jcs,jce,ips,ipe,jps,jpe;
#endif

  /* VARIABLES NEEDED FOR CYLINDRICAL COORDINATES */
#ifdef CYLINDRICAL
#ifndef ISOTHERMAL
  Real Pavgh;
#endif
  Real rinv,geom_src_d,geom_src_Vx,geom_src_Vy,geom_src_P,geom_src_By,geom_src_Bz;
  const Real *r=pG->r, *ri=pG->ri;
#ifdef FARGO
  Real Om, qshear, Mrn, Mpn, Mre, Mpe, Mrav, Mpav;
#endif
#endif /* CYLINDRICAL */
  Real g,gl,gr;
  Real lsf=1.0, rsf=1.0;

#if defined(CYLINDRICAL) && defined(FARGO)
  if (OrbitalProfile==NULL || ShearProfile==NULL)
    ath_error("[integrate_2d_sac]:  OrbitalProfile() and ShearProfile() *must* be defined.\n");
#endif

/* With particles, one more ghost cell must be updated in predict step */
#ifdef PARTICLES
  Real d1;
  il = is - 3;
  iu = ie + 3;
  jl = js - 3;
  ju = je + 3;
#else
  il = is - 2;
  iu = ie + 2;
  jl = js - 2;
  ju = je + 2;
#endif

/* Set etah=0 so first calls to flux functions do not use H-correction */
  etah = 0.0;

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
      Bxb[i] = pG->U[ks][j][i].B1cb;
      
      //B1_x1[j][i] = pG->B1i[ks][j][i]; 
#endif /* MHD */


#ifdef SAC_INTEGRATOR
        U1d[i].db  = pG->U[ks][j][i].db;

#ifdef MHD
        U1d[i].Byb = pG->U[ks][js][i].B2cb;
        U1d[i].Bzb = pG->U[ks][js][i].B3cb;
        Bxb[i] = pG->U[ks][js][i].B1cb;
#endif /* MHD */

#endif


#ifdef SMAUG_INTEGRATOR
        U1d[i].db  = pG->U[ks][j][i].db;

#ifdef MHD
        U1d[i].Byb = pG->U[ks][js][i].B2cb;
        U1d[i].Bzb = pG->U[ks][js][i].B3cb;
        Bxb[i] = pG->U[ks][js][i].B1cb;
#endif /* MHD */

#endif








#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) U1d[i].s[n] = pG->U[ks][j][i].s[n];
#endif
    }


/*--- Step 1b ------------------------------------------------------------------
 * Compute L and R states at X1-interfaces, add "MHD source terms" for 0.5*dt
 */

    for (i=is-nghost; i<=ie+nghost; i++) {
      W[i] = Cons1D_to_Prim1D(&U1d[i],&Bxc[i],&Bxb[i]);

      /* Calculate the cell-centered geometric source vector now using U^{n}
      * This will be used at the end of the integration step as a source term
      * for the cell-centered conserved variables (steps 6D,8B) */
#ifdef CYLINDRICAL 
      geom_src[j][i]  = (W[i].d+W[i].db)*SQR(W[i].Vy);   /*to here add background terms for B fields*/
#ifdef MHD
      geom_src[j][i] += 0.5*(SQR(Bxc[i]+Bxb[i]) - SQR(W[i].By+W[i].Byb) + SQR(W[i].Bz+W[i].Bzb));
#endif
#ifdef ISOTHERMAL
      geom_src[j][i] += Iso_csound2*(W[i].d+W[i].db);
#else
      geom_src[j][i] += W[i].P;
#endif
      geom_src[j][i] /= r[i];
#endif /* CYLINDRICAL */
    }

    /*lr_states(pG,W,Bxc,pG->dt,pG->dx1,il+1,iu-1,Wl,Wr,1);*/

/* Apply density floor */
   /* for (i=il+1; i<=iu; i++){
      if (Wl[i].d < d_MIN) {
        Wl[i].d = d_MIN;
      }
      if (Wr[i].d < d_MIN) {
        Wr[i].d = d_MIN;
      }
    }*/

/*#ifdef MHD
    for (i=il+1; i<=iu; i++) {
#ifdef CYLINDRICAL
      rsf = ri[i]/r[i-1];  lsf = ri[i-1]/r[i-1];
#endif
      MHD_src = (pG->U[ks][j][i-1].M2/pG->U[ks][j][i-1].d)*
                (rsf*pG->B1i[ks][j][i] - lsf*pG->B1i[ks][j][i-1])*dx1i;
      Wl[i].By += hdt*MHD_src;

#ifdef CYLINDRICAL
      rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
#endif
      MHD_src = (pG->U[ks][j][i].M2/pG->U[ks][j][i].d)*
               (rsf*pG->B1i[ks][j][i+1] - lsf*pG->B1i[ks][j][i])*dx1i;
      Wr[i].By += hdt*MHD_src;
    }
#endif*/ /* MHD */






/*--- Step 1c ------------------------------------------------------------------
 * Add source terms from static gravitational potential for 0.5*dt to L/R states
 */

  if (StaticGravPot != NULL){
    for (i=il+1; i<=iu; i++) {
      cc_pos(pG,i,j,ks,&x1,&x2,&x3);
// #ifdef CYLINDRICAL
//       gl = (*x1GravAcc)(x1vc(pG,i-1),x2,x3);
//       gr = (*x1GravAcc)(x1vc(pG,i),x2,x3);
//       gl = (*x1GravAcc)(x1-pG->dx1,x2,x3);
//       gr = (*x1GravAcc)(x1,x2,x3);
      /* APPLY GRAV. SOURCE TERMS TO V1 USING ACCELERATION FOR (dt/2) */
//       Wl[i].Vx -= hdt*gl;
//       Wr[i].Vx -= hdt*gr;
// #else
      phicr = (*StaticGravPot)( x1             ,x2,x3);
      phicl = (*StaticGravPot)((x1-    pG->dx1),x2,x3);
     /* phifc = (*StaticGravPot)((x1-0.5*pG->dx1),x2,x3);*/

      W[i].Vx -= dtodx1*(phicr -phicl);

     /* Wl[i].Vx -= dtodx1*(phifc - phicl);
      Wr[i].Vx -= dtodx1*(phicr - phifc);*/
// #endif /* CYLINDRICAL */
    }
  }


/*--- Step 1c (cont) -----------------------------------------------------------
 * Add source terms for self-gravity for 0.5*dt to L/R states
 */

#ifdef SELF_GRAVITY
    for (i=il+1; i<=iu; i++) {
      //Wl[i].Vx -= hdtodx1*(pG->Phi[ks][j][i] - pG->Phi[ks][j][i-1]);
      //Wr[i].Vx -= hdtodx1*(pG->Phi[ks][j][i] - pG->Phi[ks][j][i-1]);
      W[i].Vx -= hdtodx1*(pG->Phi[ks][j][i] - pG->Phi[ks][j][i-1]);
    }
#endif








/*--- Step 1c (cont) -----------------------------------------------------------
 * Add source terms from optically-thin cooling for 0.5*dt to L/R states
 */

#ifndef BAROTROPIC
    if (CoolingFunc != NULL){
      for (i=il+1; i<=iu; i++) {
        //coolfl = (*CoolingFunc)(Wl[i].d,Wl[i].P,(0.5*pG->dt));
        //coolfr = (*CoolingFunc)(Wr[i].d,Wr[i].P,(0.5*pG->dt));
        //coolfr = (*CoolingFunc)(Wr[i].d,Wr[i].P,(0.5*pG->dt));
        //W[i].P -= 0.5*pG->dt*Gamma_1*coolfl;
        //Wr[i].P -= 0.5*pG->dt*Gamma_1*coolfr;

        /*check cooling function*/
        coolf = (*CoolingFunc)(W[i].d+W[i].db,W[i].P,(pG->dt));
        W[i].P -= pG->dt*Gamma_1*coolf;


      }
    }
#endif /* BAROTROPIC */

/*--- Step 1c (cont) -----------------------------------------------------------
 * Add source terms for shearing box (Coriolis forces) for 0.5*dt to L/R states
 * starting with tidal gravity terms added through the ShearingBoxPot
 *    Vx source term = (dt/2)*( 2 Omega_0 Vy)
 *    Vy source term = (dt/2)*(-2 Omega_0 Vx)
 *    Vy source term = (dt/2)*((q-2) Omega_0 Vx) (with FARGO)
 * (x1,x2,x3) in code = (X,Z,Y) in 2D shearing sheet
 */

#ifdef SHEARING_BOX
    if (ShearingBoxPot != NULL){
      for (i=il+1; i<=iu; i++) {
        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        phicr = (*ShearingBoxPot)( x1             ,x2,x3);
        phicl = (*ShearingBoxPot)((x1-    pG->dx1),x2,x3);
        //phifc = (*ShearingBoxPot)((x1-0.5*pG->dx1),x2,x3);

        //Wl[i].Vx -= dtodx1*(phifc - phicl);
        //Wr[i].Vx -= dtodx1*(phicr - phifc);

        W[i].Vx -= dtodx1*(phicr - phicl);
      }
    }

    if (ShBoxCoord == xz){
      for (i=il+1; i<=iu; i++) {
        //Wl[i].Vx += pG->dt*Omega_0*W[i-1].Vz;
        //Wr[i].Vx += pG->dt*Omega_0*W[i].Vz;
        W[i].Vx += 0.5*pG->dt*Omega_0*(W[i].Vz+W[i-1].Vz);
#ifdef FARGO
        //Wl[i].Vz += hdt*(qshear-2.)*Omega_0*W[i-1].Vx;
        //Wr[i].Vz += hdt*(qshear-2.)*Omega_0*W[i].Vx;
        W[i].Vz += 0.5*hdt*(qshear-2.)*Omega_0*(W[i].Vx+W[i-1].Vx);
#else
        //Wl[i].Vz -= pG->dt*Omega_0*W[i-1].Vx;
        //Wr[i].Vz -= pG->dt*Omega_0*W[i].Vx;
        W[i].Vz -= 0.5*pG->dt*Omega_0*(W[i].Vx+W[i-1].Vx);
#endif
      }
    }

    if (ShBoxCoord == xy) {
      for (i=il+1; i<=iu; i++) {
        //Wl[i].Vx += pG->dt*Omega_0*W[i-1].Vy;
        //Wr[i].Vx += pG->dt*Omega_0*W[i].Vy;
        W[i].Vx += 0.5*pG->dt*Omega_0*(W[i].Vy+W[i-1].Vy);
#ifdef FARGO
        //Wl[i].Vy += hdt*(qshear-2.)*Omega_0*W[i-1].Vx;
        //Wr[i].Vy += hdt*(qshear-2.)*Omega_0*W[i].Vx;
        W[i].Vy += 0.5*hdt*(qshear-2.)*Omega_0*(W[i].Vx+W[i-1].Vx);
#else
        //Wl[i].Vy -= pG->dt*Omega_0*W[i-1].Vx;
        //Wr[i].Vy -= pG->dt*Omega_0*W[i].Vx;
        W[i].Vy -= 0.5*pG->dt*Omega_0*(W[i].Vx+W[i-1].Vx);
#endif
      }
    }
#endif /* SHEARING_BOX */

#if defined(CYLINDRICAL) && defined(FARGO)
    for (i=il+1; i<=iu; i++) {
      Om = (*OrbitalProfile)(r[i-1]);
      qshear = (*ShearProfile)(r[i-1]);
      W[i].Vx += 0.5*(pG->dt)*Om*W[i-1].Vy;
      W[i].Vy += 0.5*hdt*(qshear - 2.0)*Om*W[i-1].Vx;

      Om = (*OrbitalProfile)(r[i]);
      qshear = (*ShearProfile)(r[i]);
      W[i].Vx += 0.5*(pG->dt)*Om*W[i].Vy;
      W[i].Vy += 0.5*hdt*(qshear - 2.0)*Om*W[i].Vx;
    }
#endif

/*--- Step 1c (cont) -----------------------------------------------------------
 * Add source terms for Rotating Frame (Coriolis forces + Centrifugal force from 
 * Center of Mass) for 0.5*dt to L/R states
 *    Vx source term = (dt/2)*( 2 Omega_0 Vy)
 *    Vy source term = (dt/2)*(-2 Omega_0 Vx)
 *    (x1,x2,x3) in code = (X,Z,Y) in 2D shearing sheet
 */
#ifdef ROTATING_FRAME
      for (i=il+1; i<=iu; i++) {
        //cc_pos(pG,i-1,j,ks,&x1,&x2,&x3);
        //Wl[i].Vx += pG->dt*Omega_0*W[i-1].Vy - 0.5*pG->dt*SQR(Omega_0)*Rc*cos(x2);
        //Wl[i].Vy -= pG->dt*Omega_0*W[i-1].Vx - 0.5*pG->dt*SQR(Omega_0)*Rc*sin(x2);

        cc_pos(pG,i,j,ks,&x1,&x2,&x3);
        W[i].Vx += pG->dt*Omega_0*W[i].Vy - 0.5*pG->dt*SQR(Omega_0)*Rc*cos(x2);
        W[i].Vy -= pG->dt*Omega_0*W[i].Vx - 0.5*pG->dt*SQR(Omega_0)*Rc*sin(x2);
      }
#endif /*ROTATING_FRAME*/


/*--- Step 1c (cont) -----------------------------------------------------------
 * Add source terms for particle feedback for 0.5*dt to L/R states
 */

#ifdef FEEDBACK
    for (i=il+1; i<=iu; i++) {
      //d1 = 1.0/W[i-1].d;
      //Wl[i].Vx -= pG->Coup[ks][j][i-1].fb1*d1;
      //Wl[i].Vy -= pG->Coup[ks][j][i-1].fb2*d1;
      //Wl[i].Vz -= pG->Coup[ks][j][i-1].fb3*d1;

      d1 = 1.0/(W[i].d+W[i].db);
      W[i].Vx -= pG->Coup[ks][j][i].fb1*d1;
      W[i].Vy -= pG->Coup[ks][j][i].fb2*d1;
      W[i].Vz -= pG->Coup[ks][j][i].fb3*d1;

#ifndef BAROTROPIC
      //Wl[i].P += pG->Coup[ks][j][i-1].Eloss*Gamma_1;
      W[i].P += pG->Coup[ks][j][i].Eloss*Gamma_1;
#endif
    }
#endif /* FEEDBACK */

/*--- Step 1c (cont) -----------------------------------------------------------
 * Add the geometric source-terms now using cell-centered primitive
 * variables at time t^n
 */
#ifdef CYLINDRICAL
      for (i=il+1; i<=iu; i++) {


        /* right state geometric source term (uses W[i]) */
        rinv = 1.0/r[i];
        geom_src_d  = -(W[i].d+W[i].db)*W[i].Vx*rinv;
        geom_src_Vx =  SQR(W[i].Vy);
        geom_src_Vy = -W[i].Vx*W[i].Vy;
#ifdef MHD
        geom_src_Vx -= SQR(W[i].By)/(W[i].d+W[i].db);
        geom_src_Vy += Bxc[i]*W[i].By/(W[i].d+W[i].db);
        geom_src_By =  -W[i].Vy*Bxc[i]*rinv;
        geom_src_Bz =  -W[i].Vx*W[i].Bz*rinv;
#endif /* MHD */
        geom_src_Vx *= rinv;
        geom_src_Vy *= rinv;
#ifndef ISOTHERMAL
        geom_src_P  = -Gamma*W[i].P*W[i].Vx*rinv;
#endif /* ISOTHERMAL */

        /* add source term to right state */
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




/*--- Step 1d ------------------------------------------------------------------
 * Compute 1D fluxes in x1-direction, storing into 3D array
 */
    for (i=il+1; i<=iu; i++) {
      Uc_x1[j][i] = Prim1D_to_Cons1D(&W[i],&Bxc[i],&Bxb[i]);
      
/*not needed used for computing field on face*/
/*#ifdef MHD
      Bx = B1_x1[j][i];
      Bxb=0.0;//?????????????????????????
#endif*/
      fluxes(Uc_x1[j][i],Uc_x1[j][i],W[i],W[i],Bxc[i],Bxb[i],&x1Flux[j][i]);
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
      Bxb[j] = pG->U[ks][j][i].B2cb;
      //B2_x2[j][i] = pG->B2i[ks][j][i];
#endif /* MHD */


#ifdef SAC_INTEGRATOR
        U1d[j].db  = pG->U[ks][j][i].db;

#ifdef MHD
        U1d[j].Byb = pG->U[ks][j][i].B3cb;
        U1d[j].Bzb = pG->U[ks][j][i].B1cb;
        Bxb[j] = pG->U[ks][j][i].B2cb;
#endif /* MHD */

#endif


#ifdef SMAUG_INTEGRATOR
        U1d[j].db  = pG->U[ks][j][i].db;

#ifdef MHD
        U1d[j].Byb = pG->U[ks][j][i].B3cb;
        U1d[j].Bzb = pG->U[ks][j][i].B1cb;
        Bxb[j] = pG->U[ks][j][i].B2cb;
#endif /* MHD */

#endif




#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) U1d[j].s[n] = pG->U[ks][j][i].s[n];
#endif
    }

/*--- Step 2b ------------------------------------------------------------------
 * Compute L and R states at X2-interfaces, add "MHD source terms" for 0.5*dt
 */

    for (j=js-nghost; j<=je+nghost; j++) {
      W[j] = Cons1D_to_Prim1D(&U1d[j],&Bxc[j],&Bxb[j]);
    }

/*--- Step 2c ------------------------------------------------------------------
 * Add source terms from static gravitational potential for 0.5*dt to L/R states
 */

      if (StaticGravPot != NULL){
        for (j=jl+1; j<=ju; j++) {
          cc_pos(pG,i,j,ks,&x1,&x2,&x3);
          phicr = (*StaticGravPot)(x1, x2             ,x3);
          phicl = (*StaticGravPot)(x1,(x2-    pG->dx2),x3);
          //phifc = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

          W[j].Vx -= dtodx2*(phicr - phicl);
          
        }
      }



/*--- Step 2c (cont) -----------------------------------------------------------
 * Add source terms for self-gravity for 0.5*dt to L/R states
 */

#ifdef SELF_GRAVITY
    for (j=jl+1; j<=ju; j++) {
      //Wl[j].Vx -= hdtodx2*(pG->Phi[ks][j][i] - pG->Phi[ks][j-1][i]);
      //Wr[j].Vx -= hdtodx2*(pG->Phi[ks][j][i] - pG->Phi[ks][j-1][i]);
      W[j].Vx -= hdtodx2*(pG->Phi[ks][j][i] - pG->Phi[ks][j-1][i]);
    }
#endif

/*--- Step 2c (cont) -----------------------------------------------------------
 * Add source terms from optically-thin cooling for 0.5*dt to L/R states
 */

#ifndef BAROTROPIC
    if (CoolingFunc != NULL){
      for (j=jl+1; j<=ju; j++) {
        //coolfl = (*CoolingFunc)(Wl[j].d,Wl[j].P,(0.5*pG->dt));
        //coolfr = (*CoolingFunc)(Wr[j].d,Wr[j].P,(0.5*pG->dt));
        coolfc = (*CoolingFunc)(W[j].d+W[j].db,W[j].P,(0.5*pG->dt));

        //Wl[j].P -= 0.5*pG->dt*Gamma_1*coolfl;
        //Wr[j].P -= 0.5*pG->dt*Gamma_1*coolfr;
        W[j].P -= 0.5*pG->dt*Gamma_1*coolfc;
      }
    }
#endif /* BAROTROPIC */


/*--- Step 2c (cont) -----------------------------------------------------------
 * Add source terms for particle feedback for 0.5*dt to L/R states
 */

#ifdef FEEDBACK
   for (j=jl+1; j<=ju; j++) {
      //d1 = 1.0/W[j-1].d;
      //Wl[j].Vx -= pG->Coup[ks][j-1][i].fb2*d1;
      //Wl[j].Vy -= pG->Coup[ks][j-1][i].fb3*d1;
      //Wl[j].Vz -= pG->Coup[ks][j-1][i].fb1*d1;

      d1 = 1.0/(W[j].d+W[j].db);
      W[j].Vx -= pG->Coup[ks][j][i].fb2*d1;
      W[j].Vy -= pG->Coup[ks][j][i].fb3*d1;
      W[j].Vz -= pG->Coup[ks][j][i].fb1*d1;

#ifndef BAROTROPIC
      //Wl[i].P += pG->Coup[ks][j-1][i].Eloss*Gamma_1;
      W[i].P += pG->Coup[ks][j][i].Eloss*Gamma_1;
#endif
    }
#endif /* FEEDBACK */







/*--- Step 2d ------------------------------------------------------------------
 * Compute 1D fluxes in x2-direction, storing into 3D array
 */

    for (j=jl+1; j<=ju; j++) {
      Uc_x2[j][i] = Prim1D_to_Cons1D(&W[j],&Bxc[j],&Bxb[j]);
      
/*#ifdef MHD
      Bx = B2_x2[j][i];
      Bxb=0.0;//?????????????????????????
#endif*/
      fluxes(Uc_x2[j][i],Uc_x2[j][i],W[j],W[j],Bxc[j],Bxb[j],&x2Flux[j][i]);
    }
  }


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
 * Update cell-centered variables in pG using 2D x1-fluxes
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
//#ifdef MHD
//        for (j=jcs, jj=0; j<=jce+1; j++, jj++){
//          pG->CGrid[ncg].myEMF3[dim][ks][jj] = emf3[j][i];
//        }
//#endif /* MHD */
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
//#ifdef MHD
//        for (i=ics, ii=0; i<=ice+1; i++, ii++){
//          pG->CGrid[ncg].myEMF3[dim][ks][ii] = emf3[j][i];
//        }
//#endif /* MHD */
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
//#ifdef MHD
//        for (j=jps, jj=0; j<=jpe+1; j++, jj++){
//          pG->PGrid[npg].myEMF3[dim][ks][jj] = emf3[j][i];
//        }
//#endif /* MHD */
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
//#ifdef MHD
//        for (i=ips, ii=0; i<=ipe+1; i++, ii++){
//          pG->PGrid[npg].myEMF3[dim][ks][ii] = emf3[j][i];
//        }
//#endif /* MHD */
      }
    }
  }

#endif /* STATIC_MESH_REFINEMENT */


















  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_init_2d(MeshS *pM)
 *  \brief Allocate temporary integration arrays 
*/
void integrate_init_2d(MeshS *pM)
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
  //if ((Bxi = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;
  if ((Bxb = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;


  if ((U1d= (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((W  = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;

  if ((x1Flux   =(Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS)))==NULL)
    goto on_error;
  if ((x2Flux   =(Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS)))==NULL)
    goto on_error;



#ifndef CYLINDRICAL
#ifndef MHD
#ifndef PARTICLES
  if((StaticGravPot != NULL) || (CoolingFunc != NULL))
#endif
#endif
#endif
  {
  if ((dhalf = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
  if ((phalf = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL)
    goto on_error;
  }

  /* data structures for cylindrical coordinates */
#ifdef CYLINDRICAL
  if ((geom_src = (Real**)calloc_2d_array(size2, size1, sizeof(Real))) == NULL) 
    goto on_error;
#endif



  return;

  on_error:
    integrate_destruct();
    ath_error("[integrate_init]: malloc returned a NULL pointer\n");
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_destruct_2d(void)
 *  \brief Free temporary integration arrays 
 */
void integrate_destruct_2d(void)
{
/*#ifdef MHD
  if (emf3    != NULL) free_2d_array(emf3);
  if (emf3_cc != NULL) free_2d_array(emf3_cc);
#endif *//* MHD */
/*#ifdef H_CORRECTION
  if (eta1 != NULL) free_2d_array(eta1);
  if (eta2 != NULL) free_2d_array(eta2);
#endif *//* H_CORRECTION */
  if (Bxc != NULL) free(Bxc);
  //if (Bxi != NULL) free(Bxi);
  if (Bxb != NULL) free(Bxb);
/*#ifdef MHD
  if (B1_x1Face != NULL) free_2d_array(B1_x1Face);
  if (B2_x2Face != NULL) free_2d_array(B2_x2Face);
#endif *//* MHD */

  if (U1d      != NULL) free(U1d);
 /* if (Ul       != NULL) free(Ul);
  if (Ur       != NULL) free(Ur);*/
  if (W        != NULL) free(W);
 /* if (Wl       != NULL) free(Wl);
  if (Wr       != NULL) free(Wr);*/

/*  if (Ul_x1Face != NULL) free_2d_array(Ul_x1Face);
  if (Ur_x1Face != NULL) free_2d_array(Ur_x1Face);
  if (Ul_x2Face != NULL) free_2d_array(Ul_x2Face);
  if (Ur_x2Face != NULL) free_2d_array(Ur_x2Face);*/
  if (x1Flux    != NULL) free_2d_array(x1Flux);
  if (x2Flux    != NULL) free_2d_array(x2Flux);
  if (dhalf     != NULL) free_2d_array(dhalf);
  if (phalf     != NULL) free_2d_array(phalf);

  /* data structures for cylindrical coordinates */
#ifdef CYLINDRICAL
  if (geom_src  != NULL) free_2d_array(geom_src);
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
