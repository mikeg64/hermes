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
static ConsS ***Uinit=NULL; /*Uinit used to store initial fields*/

/* The interface magnetic fields and emfs */
/*not needed for sac*/
/*#ifdef MHD
static Real **B1_x1Face=NULL, **B2_x2Face=NULL;
static Real **emf3=NULL, **emf3_cc=NULL;
#endif*/ /* MHD */





/* 1D scratch vectors used by lr_states and flux functions */
static Real *Bxc=NULL,  *Bxb=NULL, *temp=NULL, *grad=NULL;
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
static void computemaxc(ConsS ***Uint, GridS *pG, int dim);


static void hyperdifviscr(int fieldi,int dim,ConsS ***Uint, GridS *pG);
static void hyperdifviscl(int fieldi,int dim,ConsS ***Uint, GridS *pG);

static void hyperdifrhosource(int dim,Real dt,ConsS ***Uint, GridS *pG);
/*static void hyperdifesource(int dim,Real dt,ConsS ***Uint, GridS *pG); 
static void hyperdifmomsource(int field,int dim,int ii,int ii0,Real dt,ConsS ***Uint, GridS *pG);
static void hyperdifmomsourcene(int field,int dim,int ii,int ii0, Real dt,ConsS ***Uint, GridS *pG);*/

#ifdef MHD
//static void hyperdifbsource(int ii,int ii0, Real dt, ConsS ***Uint, GridS *pG);
//static void hyperdifbsourcene(int ii,int ii0, Real dt, ConsS ***Uint, GridS *pG);
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
int jj1, jj, jj0;
int fieldi; /*integers map to following index rho, mom1, mom2, energy, b1, b2,energyb,rhob,b1b,b2b*/

int size1,size2;

size1=1+ie+2*nghost-is;
size2=1+je+2*nghost-js;



/*=== STEP 1: Compute L/R x1-interface states and 1D x1-Fluxes ===============*/


    int k=0;
  
    for (j=jl; j<=ju; j++) {
      for (i=il; i<=iu; i++) {

        Uinit[k][j][i].d  = pG->U[k][j][i].d;
        Uinit[k][j][i].M1 = pG->U[k][j][i].M1;
        Uinit[k][j][i].M2 = pG->U[k][j][i].M2;
        Uinit[k][j][i].M3 = pG->U[k][j][i].M3;
#ifndef BAROTROPIC
        Uinit[k][j][i].E  = pG->U[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
        Uinit[k][j][i].B1c = pG->U[k][j][i].B1c;
        Uinit[k][j][i].B2c = pG->U[k][j][i].B2c;
        Uinit[k][j][i].B3c = pG->U[k][j][i].B3c;
        //Bxc[i] = pG->U[k][j][i].B1c;
        //Bxb[i] = pG->B1cb[k][j][i];
        //B1_x1Face[k][j][i] = pG->B1i[k][j][i];
#endif /* MHD */


#ifdef SAC_INTEGRATOR
        Uinit[k][j][i].db  = pG->U[k][j][i].db;
        Uinit[k][j][i].B1cb = pG->U[k][j][i].B1cb;        
        Uinit[k][j][i].B2cb = pG->U[k][j][i].B2cb;
        Uinit[k][j][i].B3cb = pG->U[k][j][i].B3cb;
#endif
#ifdef SMAUG_INTEGRATOR
        Uinit[k][j][i].db  = pG->U[k][j][i].db;
        Uinit[k][j][i].B1cb = pG->U[k][j][i].B1cb;        
        Uinit[k][j][i].B2cb = pG->U[k][j][i].B2cb;
        Uinit[k][j][i].B3cb = pG->U[k][j][i].B3cb;
#endif


#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) Uinit[k][j][i].s[n] = pG->U[k][j][i].s[n];
#endif


}
}







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
//computemaxc(Uinit,pG);

//density contribution
for(dim=0; dim<2; dim++) //each direction
{


computemaxc(Uinit,pG,dim);
hyperdifviscr(rho,dim,Uinit, pG);
hyperdifviscl(rho,dim,Uinit, pG);
//hyperdifvisc1ir
//hyperdifvisc1il
//int dim,Real dt,ConsS ***Uint, GridS *pG
hyperdifrhosource(dim,pG->dt,Uinit, pG) ;
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
				           ii0=fieldi;  //f is field
		                  }
		                  else
		                  {
				           ii=fieldi;
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
				           ii0=fieldi;  //f is field
		                  }
		                  else
		                  {
				           ii=fieldi;
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
      }
    }
  }

  size1 = size1 + 2*nghost;
  size2 = size2 + 2*nghost;
 
  nmax = MAX(size1,size2);

/*refer to material  integrate_2d_ctu.c*/
  if ((Bxc = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;
  //if ((Bxi = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;
  if ((Bxb = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;
  if ((temp = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;
  if ((grad = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;

  if ((U1d= (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;

  if ((Uc_x1= (Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS))) == NULL) goto on_error;
  if ((Uc_x2= (Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS))) == NULL) goto on_error;



  if ((W  = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;

  if ((x1Flux   =(Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS)))==NULL)
    goto on_error;
  if ((x2Flux   =(Cons1DS**)calloc_2d_array(size2,size1,sizeof(Cons1DS)))==NULL)
    goto on_error;

  if ((Uinit   =(ConsS***)calloc_3d_array(1,size2,size1,sizeof(ConsS)))
    == NULL) goto on_error;

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
  if (temp != NULL) free(temp);
  if (grad != NULL) free(grad);
/*#ifdef MHD
  if (B1_x1Face != NULL) free_2d_array(B1_x1Face);
  if (B2_x2Face != NULL) free_2d_array(B2_x2Face);
#endif *//* MHD */

  if (U1d      != NULL) free(U1d);
 /* if (Ul       != NULL) free(Ul);
  if (Ur       != NULL) free(Ur);*/

  if (Uc_x1      != NULL) free(Uc_x1);
  if (Uc_x2      != NULL) free(Uc_x2);


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

  if (Uinit    != NULL) free_3d_array(Uinit);


  /* data structures for cylindrical coordinates */
#ifdef CYLINDRICAL
  if (geom_src  != NULL) free_2d_array(geom_src);
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
		rhototsq*=rhotot*rhotot;
		pthermal=Uinit[i3][i2][i1].E -((Uinit[i3][i2][i1].M1*Uinit[i3][i2][i1].M1+Uinit[i3][i2][i1].M2*Uinit[i3][i2][i1].M2+Uinit[i3][i2][i1].M3*Uinit[i3][i2][i1].M3)/rhotot)  ;
		pthermal+=0.5*((Uinit[i3][i2][i1].B1c*Uinit[i3][i2][i1].B1c+Uinit[i3][i2][i1].B2c*Uinit[i3][i2][i1].B2c+Uinit[i3][i2][i1].B3c*Uinit[i3][i2][i1].B3c));
		pthermal+=0.5*((Uinit[i3][i2][i1].B1c*Uinit[i3][i2][i1].B1cb+Uinit[i3][i2][i1].B2c*Uinit[i3][i2][i1].B2cb+Uinit[i3][i2][i1].B3c*Uinit[i3][i2][i1].B3cb));
		pthermal*=(Gamma_1-1);
		cs2=Gamma_1*pthermal+(Gamma_1-1)*(Uinit[i3][i2][i1].Eb-0.5*(    ((Uinit[i3][i2][i1].B1cb*Uinit[i3][i2][i1].B1cb+Uinit[i3][i2][i1].B2cb*Uinit[i3][i2][i1].B2cb+Uinit[i3][i2][i1].B3cb*Uinit[i3][i2][i1].B3cb)  )));
		cs2/=rhototsq;

		pG->Hv[i3][i2][i1].csound=sqrt(cs2);
                //cmax=MAX(cmax,pG->Hv[i3][i2][i1].csound)

		cfast2=cs2+((Uinit[i3][i2][i1].B1c+Uinit[i3][i2][i1].B1cb)*(Uinit[i3][i2][i1].B1c+Uinit[i3][i2][i1].B1cb)+
                        (Uinit[i3][i2][i1].B2c+Uinit[i3][i2][i1].B2cb)*(Uinit[i3][i2][i1].B2c+Uinit[i3][i2][i1].B2cb)+
			(Uinit[i3][i2][i1].B3c+Uinit[i3][i2][i1].B3cb)*(Uinit[i3][i2][i1].B3c+Uinit[i3][i2][i1].B3cb))/(rhotot);




		pG->Hv[i3][i2][i1].cfast=sqrt(cfast2);
		//cmax=MAX(cmax,pG->Hv[i3][i2][i1].cfast)

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
		case b1:
			fieldd[i3][i2][i1]=Uinit[i3][i2][i1].B1c;
		break;
		case b2:
			fieldd[i3][i2][i1]=Uinit[i3][i2][i1].B2c;
		break;
		case b3:
			fieldd[i3][i2][i1]=Uinit[i3][i2][i1].B3c;
		break;
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

		pG->Hv[i3][i2][i1].hdnur[dim][fieldi]=0.0;

	       if(fieldi==energy)
		wtemp1[i3][i2][i1]=fieldd[i3][i2][i1]-0.5*((Uinit[i3][i2][i1].B1c*Uinit[i3][i2][i1].B1c+Uinit[i3][i2][i1].B2c*Uinit[i3][i2][i1].B2c+Uinit[i3][i2][i1].B3c*Uinit[i3][i2][i1].B3c)
	+(Uinit[i3][i2][i1].M1*Uinit[i3][i2][i1].M1+Uinit[i3][i2][i1].M2*Uinit[i3][i2][i1].M2+Uinit[i3][i2][i1].M3*Uinit[i3][i2][i1].M3)/(Uinit[i3][i2][i1].d+Uinit[i3][i2][i1].db ));       
	       else
	       {
		wtemp1[i3][i2][i1]=fieldd[i3][i2][i1];
		if((fieldi ==mom1 || fieldi == mom2 || fieldi == mom3))
			wtemp1[i3][i2][i1]=fieldd[i3][i2][i1]/(Uinit[i3][i2][i1].d+Uinit[i3][i2][i1].db);

		}

		//comment removed below to test mpi 29/10/2013
		tmpnu[i3+1][i2+1][i1+1]=wtemp1[i3][i2][i1];


	}
}
}

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
		   d1[i3+1][i2+1][i1+1]=fabs((tmpnu[i3+(dim==2)+1][i2+(dim==1)+1][i1+(dim==0)+1] - tmpnu[i3+1][i2+1][i1+1] ));
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



        
	if (wtemp1 != NULL) free(wtemp1);
	if (wtemp2 != NULL) free(wtemp2);
	if (wtemp3 != NULL) free(wtemp3);
	if (fieldd != NULL) free(fieldd);
	if (tmpnu != NULL) free(tmpnu);
	if (d3 != NULL) free(d3);
	if (d1 != NULL) free(d1);


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
		case b1:
			fieldd[i3][i2][i1]=Uinit[i3][i2][i1].B1c;
		break;
		case b2:
			fieldd[i3][i2][i1]=Uinit[i3][i2][i1].B2c;
		break;
		case b3:
			fieldd[i3][i2][i1]=Uinit[i3][i2][i1].B3c;
		break;
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
		wtemp1[i3][i2][i1]=fieldd[i3][i2][i1]-0.5*((Uinit[i3][i2][i1].B1c*Uinit[i3][i2][i1].B1c+Uinit[i3][i2][i1].B2c*Uinit[i3][i2][i1].B2c+Uinit[i3][i2][i1].B3c*Uinit[i3][i2][i1].B3c)
	+(Uinit[i3][i2][i1].M1*Uinit[i3][i2][i1].M1+Uinit[i3][i2][i1].M2*Uinit[i3][i2][i1].M2+Uinit[i3][i2][i1].M3*Uinit[i3][i2][i1].M3)/(Uinit[i3][i2][i1].d+Uinit[i3][i2][i1].db ));       
	       else
	       {
		wtemp1[i3][i2][i1]=fieldd[i3][i2][i1];
		if((fieldi ==mom1 || fieldi == mom2 || fieldi == mom3))
			wtemp1[i3][i2][i1]=fieldd[i3][i2][i1]/(Uinit[i3][i2][i1].d+Uinit[i3][i2][i1].db);

		}

		//comment removed below to test mpi 29/10/2013
		tmpnu[i3+1][i2+1][i1+1]=wtemp1[i3][i2][i1];


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
		   d1[i3+1][i2+1][i1+1]=fabs(( tmpnu[i3+1][i2+1][i1+1] -   tmpnu[i3+(dim==2)+1][i2+(dim==1)+1][i1+(dim==0)+1]  ));
}
}
}

   /*to here*/
        maxt2=0.0;
        maxt1=0.0;
  for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {

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
	if (tmpnu != NULL) free(tmpnu);
	if (d3 != NULL) free(d3);
	if (d1 != NULL) free(d1);
	if (fieldd != NULL) free(fieldd);

	return;


}





static void hyperdifrhosource(int dim,Real dt,ConsS ***Uint, GridS *pG)
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

	fieldd = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	wtempr = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	wtempl = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	wtemp3 = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmp = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));
	tmp2 = (Real***)calloc_3d_array(n3z, n2z, n1z, sizeof(Real));



     //CALL setnu(w,rho_,idim,ixOmin1,ixOmin2,ixOmax1,ixOmax2,nuR,nuL)


 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			fieldd[i3][i2][i1]=Uinit[i3][i2][i1].d;
					}
				}
			}

     //CALL gradient1L(tmp,ixmin1,ixmin2,ixmax1,ixmax2,idim,tmp2)
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1l(fieldd[i3][i2], n1z,pG->dx1,tmp2[i3][i2]);
				}
			}
      
     

/*nur=pG->Hv[i3][i2][i1].hdnur[dim][fieldi];
nul=pG->Hv[i3][i2][i1].hdnur[dim][fieldi];*/

     /*tmpL(ixImin1:ixImax1,ixImin2:ixImax2)=(nuL(ixImin1:ixImax1,&
        ixImin2:ixImax2)+nushk(ixImin1:ixImax1,ixImin2:ixImax2,idim))&
        *tmp2(ixImin1:ixImax1,ixImin2:ixImax2) */ 

for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			wtempl[i3][i2][i1]=(pG->Hv[i3][i2][i1].hdnul[dim][fieldi])*tmp2[i3][i2][i1];
					}
				}
			}

        
     //CALL gradient1R(tmp,ixmin1,ixmin2,ixmax1,ixmax2,idim,tmp2)
 for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
			gradient1r(fieldd[i3][i2], n1z,pG->dx1,tmp2[i3][i2]);
				}
			}




     /*tmpR(ixImin1:ixImax1,ixImin2:ixImax2)=(nuR(ixImin1:ixImax1,&
        ixImin2:ixImax2)+nushk(ixImin1:ixImax1,ixImin2:ixImax2,idim))&
        *tmp2(ixImin1:ixImax1,ixImin2:ixImax2)*/
for (i3=kl; i3<=ku; i3++) {
    for (i2=jl; i2<=ju; i2++) {
    	for (i1=il; i1<=iu; i1++) {
			wtempr[i3][i2][i1]  =(pG->Hv[i3][i2][i1].hdnur[dim][fieldi])*tmp2[i3][i2][i1];
					}
				}
			}



     /*wnew(ixImin1:ixImax1,ixImin2:ixImax2,rho_)=wnew(ixImin1:ixImax1,&
        ixImin2:ixImax2,rho_)+(tmpR(ixImin1:ixImax1,ixImin2:ixImax2)&
        -tmpL(ixImin1:ixImax1,ixImin2:ixImax2))/dx(ixImin1:ixImax1,&
        ixImin2:ixImax2,idim)*qdt*/

     pG->U[ks][j][i].d  += dtodx1*(wtempr-wtempl);


	if (wtempr != NULL) free(wtempr);
	if (wtempl != NULL) free(wtempl);
	if (wtemp3 != NULL) free(wtemp3);
	if (tmp != NULL) free(tmp);
	if (tmp2 != NULL) free(tmp2);
	if (fieldd != NULL) free(fieldd);



	return;
}

/*
static void hyperdifesource(int dim,Real dt,ConsS ***Uint, GridS *pG)
{

	return;
}

 
static void hyperdifmomsource(int field,int dim,int ii,int ii0,Real dt,ConsS ***Uint, GridS *pG)
{

	return;
}


static void hyperdifmomsourcene(int field,int dim,int ii,int ii0, Real dt,ConsS ***Uint, GridS *pG)
{

	return;
}

*/

#ifdef MHD

/*static void hyperdifbsource(int ii,int ii0, Real dt,ConsS ***Uint, GridS *pG)
{

	return;
}


static void hyperdifbsourcene(int ii,int ii0, Real dt,ConsS ***Uint, GridS *pG)
{

	return;
}*/


#endif /* MHD */


#endif /* SAC_INTEGRATOR */
