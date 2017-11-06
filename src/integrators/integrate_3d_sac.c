#include "../copyright.h"
/*============================================================================*/
/*! \file integrate_3d_sac.c
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

#if defined(SAC_INTEGRATOR)
#ifdef SPECIAL_RELATIVITY
#error : The SAC integrator cannot be used for special relativity.
#endif /* SPECIAL_RELATIVITY */




/* The L/R states of conserved variables and fluxes at each cell face */
//static Cons1DS ***U1d=NULL, ***Ur_x1Face=NULL;
//static Cons1DS ***Ul_x2Face=NULL, ***Ur_x2Face=NULL;
//static Cons1DS ***Ul_x3Face=NULL, ***Ur_x3Face=NULL;
Cons1DS ***x1Flux=NULL, ***x2Flux=NULL, ***x3Flux=NULL;
static Cons1DS ***Uc_x1=NULL;
static Cons1DS ***Uc_x2=NULL;
static Cons1DS ***Uc_x3=NULL;


//static Cons1DS ***Uc=NULL;


/* The interface magnetic fields and emfs */
//#ifdef MHD
//static Real ***B1_x1Face=NULL, ***B2_x2Face=NULL, ***B3_x3Face=NULL;
//Real ***emf1=NULL, ***emf2=NULL, ***emf3=NULL;
//static Real ***emf1_cc=NULL, ***emf2_cc=NULL, ***emf3_cc=NULL;
//#endif /* MHD */

/* 1D scratch vectors used by lr_states and flux functions */
static Real *Bxc=NULL, *Bxb=NULL;
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

void integrate_3d_sac(DomainS *pD)
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
#ifndef BAROTROPIC
  Real coolfl,coolfr,coolf,coolfc,Eh=0.0;
#endif

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
  Real gxc,gyc,gzc,gxl,gxr,gyl,gyr,gzl,gzr,flx_m1l,flx_m1r,flx_m2l,flx_m2r,flx_m3l,flx_m3r;
#endif
#ifdef FEEDBACK
  Real dt1 = 1.0/pG->dt;
#endif


#ifdef STATIC_MESH_REFINEMENT
  int ncg,npg,dim;
  int ii,ics,ice,jj,jcs,jce,kk,kcs,kce,ips,ipe,jps,jpe,kps,kpe;
#endif


  Real g,gc,gl,gr;
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
        //Bxb[i] = pG->B1cb[k][j][i];
        //B1_x1Face[k][j][i] = pG->B1i[k][j][i];
#endif /* MHD */


#ifdef BKG
        U1d[i].db  = pG->U[k][j][i].db;
        Bxb[i] = pG->U[k][j][i].B1cb;
        U1d[i].Byb = pG->U[k][j][i].B2cb;
        U1d[i].Bzb = pG->U[k][j][i].B3cb;
#endif


#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) U1d[i].s[n] = pG->U[k][j][i].s[n];
#endif
      }





/*--- Step 1b ------------------------------------------------------------------
 * Compute L and R states at X1-interfaces, add "MHD source terms" for 0.5*dt
 */

   for (i=is-nghost; i<=ie+nghost; i++) {
      W[i] = Cons1D_to_Prim1D(&U1d[i],&Bxc[i],&Bxb[i]);

        /* Calculate the cell-centered geometric source vector now using U^{n}
         * This will be used at the end of the integration step as a source term
         * for the cell-centered conserved variables (steps 6D,7D,8B) */
#ifdef CYLINDRICAL
        geom_src[k][j][i]  = (W[i].d+W[i].db)*SQR(W[i].Vy);
#ifdef MHD
        geom_src[k][j][i] += 0.5*(SQR(Bxc[i]+Bxb[i]) - SQR(W[i].By+W[i].Byb) + SQR(W[i].Bz+W[i].Bzb));
#endif
#ifdef ISOTHERMAL
        geom_src[k][j][i] += Iso_csound2*(W[i].d+W[i].db);
#else
        geom_src[k][j][i] += W[i].P;
#endif
        geom_src[k][j][i] /= r[i];
#endif /* CYLINDRICAL */
     }



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

          //Wl[i].Vx -= hdt*gl;
          //Wr[i].Vx -= hdt*gr;
          W[i].Vx -= 0.5*hdt*(gl+gr);
        }
      }





/*--- Step 1c (cont) -----------------------------------------------------------
 * Add source terms for self-gravity for 0.5*dt to L/R states
 */

#ifdef SELF_GRAVITY
      for (i=il+1; i<=iu; i++) {
        //Wl[i].Vx -= q1*(pG->Phi[k][j][i] - pG->Phi[k][j][i-1]);
        //Wr[i].Vx -= q1*(pG->Phi[k][j][i] - pG->Phi[k][j][i-1]);
        W[i].Vx -= q1*(pG->Phi[k][j][i] - pG->Phi[k][j][i-1]);
      }
#endif



/*--- Step 1c (cont) -----------------------------------------------------------
 * Add source terms from optically-thin cooling for 0.5*dt to L/R states
 */



#ifndef BAROTROPIC
      if (CoolingFunc != NULL){
        for (i=il+1; i<=iu; i++) {
          /*coolfl = (*CoolingFunc)(Wl[i].d,Wl[i].P,(0.5*pG->dt));
          coolfr = (*CoolingFunc)(Wr[i].d,Wr[i].P,(0.5*pG->dt));

          Wl[i].P -= 0.5*pG->dt*Gamma_1*coolfl;
          Wr[i].P -= 0.5*pG->dt*Gamma_1*coolfr;*/

          coolfc = (*CoolingFunc)(W[i].d+W[i].db,W[i].P,(0.5*pG->dt));
          W[i].P -= 0.5*pG->dt*Gamma_1*coolfc;
        }
      }
#endif /* BAROTROPIC */



/*--- Step 1c (cont) -----------------------------------------------------------
 * Add source terms for shearing box (Coriolis forces) for 0.5*dt to L/R states
 * starting with tidal gravity terms added through the ShearingBoxPot
 *    Vx source term = (dt/2)*( 2 Omega_0 Vy)
 *    Vy source term = (dt/2)*(-2 Omega_0 Vx)
 *    Vy source term = (dt/2)*((q-2) Omega_0 Vx) (with FARGO)
 */

#ifdef SHEARING_BOX
      if (ShearingBoxPot != NULL){
        for (i=il+1; i<=iu; i++) {
          cc_pos(pG,i,j,k,&x1,&x2,&x3);
          phicr = (*ShearingBoxPot)( x1             ,x2,x3);
          phicl = (*ShearingBoxPot)((x1-    pG->dx1),x2,x3);
          phifc = (*ShearingBoxPot)((x1-0.5*pG->dx1),x2,x3);

          //Wl[i].Vx -= dtodx1*(phifc - phicl);
          //Wr[i].Vx -= dtodx1*(phicr - phifc);
          W[i].Vx -= 0.5*dtodx1*(phicr - phicl);
        }
      }

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

 /* Add Coriolis Force and Centrifugal Force in Rotating Frame */
#ifdef ROTATING_FRAME
      for (i=il+1; i<=iu; i++) {
        cc_pos(pG,i-1,j,k,&x1,&x2,&x3);
        W[i].Vx += 0.5*pG->dt*Omega_0*W[i-1].Vy - 0.25*pG->dt*SQR(Omega_0)*Rc*cos(x2);
        W[i].Vy -= 0.5*pG->dt*Omega_0*W[i-1].Vx - 0.25*pG->dt*SQR(Omega_0)*Rc*sin(x2);

        cc_pos(pG,i,j,k,&x1,&x2,&x3);
        W[i].Vx += 0.5*pG->dt*Omega_0*W[i].Vy - 0.25*pG->dt*SQR(Omega_0)*Rc*cos(x2);
        W[i].Vy -= 0.5*pG->dt*Omega_0*W[i].Vx - 0.25*pG->dt*SQR(Omega_0)*Rc*sin(x2);
      }
#endif /*ROTATING_FRAME*/

/*--- Step 1c (cont) -----------------------------------------------------------
 * Add source terms for particle feedback for 0.5*dt to L/R states
 */


#ifdef FEEDBACK
    for (i=il+1; i<=iu; i++) {
      d1 = 1.0/W[i-1].d;
      W[i].Vx -= 0.5*pG->Coup[k][j][i-1].fb1*d1;
      W[i].Vy -= 0.5*pG->Coup[k][j][i-1].fb2*d1;
      W[i].Vz -= 0.5*pG->Coup[k][j][i-1].fb3*d1;

      d1 = 1.0/W[i].d;
      W[i].Vx -= 0.5*pG->Coup[k][j][i].fb1*d1;
      W[i].Vy -= 0.5*pG->Coup[k][j][i].fb2*d1;
      W[i].Vz -= 0.5*pG->Coup[k][j][i].fb3*d1;

#ifndef BAROTROPIC
      W[i].P += 0.5*pG->Coup[k][j][i-1].Eloss*Gamma_1;
      W[i].P += 0.5*pG->Coup[k][j][i].Eloss*Gamma_1;
#endif

    }
#endif /* FEEDBACK */


/*--- Step 1c (cont) -----------------------------------------------------------
 * Add the geometric source-terms now using cell-centered primitive
 * variables at time t^n
 */
#ifdef CYLINDRICAL
      for (i=is-1; i<=ie+2; i++) {

        /* left state geometric source term (uses W[i-1]) */
        rinv = 1.0/r[i-1];
        geom_src_d  = -(W[i-1].d+W[i-1].db)*W[i-1].Vx*rinv;
        geom_src_Vx =  SQR(W[i-1].Vy);
        geom_src_Vy = -W[i-1].Vx*W[i-1].Vy;
#ifdef MHD
        geom_src_Vx -= SQR(W[i-1].By)/(W[i-1].d+W[i-1].db);
        geom_src_Vy += (Bxc[i-1]+Bxb[i-1])*(W[i-1].By+W[i-1].Byb)/(W[i-1].d+W[i-1].db);
        geom_src_By =  -W[i-1].Vy*(Bxc[i-1]+Bxb[i-1])*rinv;
        geom_src_Bz =  -W[i-1].Vx*(W[i-1].Bz+W[i-1].Bzb)*rinv;
#endif
        geom_src_Vx *= rinv;
        geom_src_Vy *= rinv;
#ifndef ISOTHERMAL
        geom_src_P  = -Gamma*W[i-1].P*W[i-1].Vx*rinv;
#endif

        /* add source term to left state */
        W[i].d  += 0.5*hdt*geom_src_d;
        W[i].Vx += 0.5*hdt*geom_src_Vx;
        W[i].Vy += 0.5*hdt*geom_src_Vy;
#ifdef MHD
        W[i].By += 0.5*hdt*geom_src_By;
        W[i].Bz += 0.5*hdt*geom_src_Bz;
#endif
#ifndef ISOTHERMAL
        W[i].P  += 0.5*hdt*geom_src_P;
#endif

        /* right state geometric source term (uses W[i]) */
        rinv = 1.0/r[i];
        geom_src_d  = -(W[i].d+W[i].db)*W[i].Vx*rinv;
        geom_src_Vx =  SQR(W[i].Vy);
        geom_src_Vy = -W[i].Vx*W[i].Vy;
#ifdef MHD
        geom_src_Vx -= SQR(W[i].By+W[i].Byb)/(W[i].d+W[i].db);
        geom_src_Vy += (Bxc[i]+Bxb[i])*(W[i].By+W[i].Byb)/(W[i].d+W[i].db);
        geom_src_By =  -W[i].Vy*(Bxc[i]+Bxb[i])*rinv;
        geom_src_Bz =  -W[i].Vx*(W[i].Bz+W[i].Bzb)*rinv;
#endif
        geom_src_Vx *= rinv;
        geom_src_Vy *= rinv;
#ifndef ISOTHERMAL
        geom_src_P  = -Gamma*W[i].P*W[i].Vx*rinv;
#endif

        /* add source term to right state */
        W[i].d  += 0.5*hdt*geom_src_d;
        W[i].Vx += 0.5*hdt*geom_src_Vx;
        W[i].Vy += 0.5*hdt*geom_src_Vy;
#ifdef MHD
        W[i].By += 0.5*hdt*geom_src_By;
        W[i].Bz += 0.5*hdt*geom_src_Bz;
#endif
#ifndef ISOTHERMAL
        W[i].P  += 0.5*hdt*geom_src_P;
#endif
      }
#endif /* CYLINDRICAL */




/*--- Step 1d ------------------------------------------------------------------
 * Compute 1D fluxes in x1-direction, storing into 3D array
 */
    for (i=il+1; i<=iu; i++) {
      Uc_x1[k][j][i] = Prim1D_to_Cons1D(&W[i],&Bxc[i],&Bxb[i]);
      

//#ifdef MHD
//      Bx = B1_x1[j][i];
//      Bxb=0.0;//?????????????????????????
//#endif
      fluxes(Uc_x1[k][j][i],Uc_x1[k][j][i],W[i],W[i],Bxc,Bxb,&x1Flux[k][j][i]);
    }
  }
}

/*=== STEP 2: Compute L/R x2-interface states and 1D x2-Fluxes ===============*/

/*--- Step 2a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M2, M3, M1, E, B3c, B1c, s[n])
 */

  for (k=kl; k<=ku; k++) {
    for (i=il; i<=iu; i++) {
#ifdef CYLINDRICAL
      rsf = ri[i+1]/r[i];  lsf = ri[i]/r[i];
      dx2 = r[i]*pG->dx2;
      dtodx2 = pG->dt/dx2;
#endif
    for (j=js-nghost; j<=je+nghost; j++) {
      U1d[j].d  = pG->U[k][j][i].d;
      U1d[j].Mx = pG->U[k][j][i].M2;
      U1d[j].My = pG->U[k][j][i].M3;
      U1d[j].Mz = pG->U[k][j][i].M1;
#ifndef BAROTROPIC
      U1d[j].E  = pG->U[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
      U1d[j].By = pG->U[k][j][i].B3c;
      U1d[j].Bz = pG->U[k][j][i].B1c;
      Bxc[j] = pG->U[k][j][i].B2c;
      Bxb[j] = pG->B2b[k][j][i];
//      B2_x2[j][i] = pG->B2i[k][j][i];
#endif /* MHD */


#ifdef BKG
        U1d[j].db  = pG->U[k][j][i].db;

#ifdef MHD
        Bxb[j] = pG->U[k][j][i].B2cb;
        U1d[j].Byb = pG->U[k][j][i].B3cb;
        U1d[j].Bzb = pG->U[k][j][i].B1cb;
#endif /* MHD */
#endif



#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++) U1d[j].s[n] = pG->U[k][j][i].s[n];
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
          cc_pos(pG,i,j,k,&x1,&x2,&x3);
          phicr = (*StaticGravPot)(x1, x2             ,x3);
          phicl = (*StaticGravPot)(x1,(x2-    pG->dx2),x3);
          //phifc = (*StaticGravPot)(x1,(x2-0.5*pG->dx2),x3);

          //Wl[j].Vx -= dtodx2*(phifc - phicl);
          W[j].Vx -= 0.5*dtodx2*(phicr - phicl);
        }
      }

/*--- Step 2c (cont) -----------------------------------------------------------
 * Add source terms for self-gravity for 0.5*dt to L/R states
 */

#ifdef SELF_GRAVITY
      for (j=jl+1; j<=ju; j++) {
        //Wl[j].Vx -= q2*(pG->Phi[k][j][i] - pG->Phi[k][j-1][i]);
        W[j].Vx -= q2*(pG->Phi[k][j][i] - pG->Phi[k][j-1][i]);
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
      d1 = 1.0/(W[j-1].d+W[j-1].db);
      W[j].Vx -= 0.5*pG->Coup[k][j-1][i].fb2*d1;
      W[j].Vy -= 0.5*pG->Coup[k][j-1][i].fb3*d1;
      W[j].Vz -= 0.5*pG->Coup[k][j-1][i].fb1*d1;

      d1 = 1.0/(W[j].d+W[j].db);
      W[j].Vx -= 0.5*pG->Coup[k][j][i].fb2*d1;
      W[j].Vy -= 0.5*pG->Coup[k][j][i].fb3*d1;
      W[j].Vz -= 0.5*pG->Coup[k][j][i].fb1*d1;

#ifndef BAROTROPIC
      W[i].P += 0.5*pG->Coup[k][j-1][i].Eloss*Gamma_1;
      W[i].P += 0.5*pG->Coup[k][j][i].Eloss*Gamma_1;
#endif

    }
#endif /* FEEDBACK */





/*--- Step 2d ------------------------------------------------------------------
 * Compute 1D fluxes in x2-direction, storing into 3D array
 */

    for (j=jl+1; j<=ju; j++) {
      Uc_x2[k][j][i] = Prim1D_to_Cons1D(&W[j],&Bxc[j],&Bxb[j]);
      
//#ifdef MHD
//      Bx = B2_x2[j][i];
//       Bxb=0.0;//?????????????????????????
//#endif
      fluxes(Uc_x2[k][j][i],Uc_x2[k][j][i],W[j],W[j],Bxc,Bxb,&x2Flux[k][j][i]);
    }
  }

}


/*=== STEP 3: Compute L/R x3-interface states and 1D x3-Fluxes ===============*/

/*--- Step 3a ------------------------------------------------------------------
 * Load 1D vector of conserved variables;
 * U1d = (d, M3, M1, M2, E, B1c, B2c, s[n])
 */

  for (j=jl; j<=ju; j++) {
    for (i=il; i<=iu; i++) {
      for (k=ks-nghost; k<=ke+nghost; k++) {
        U1d[k].d  = pG->U[k][j][i].d;
        U1d[k].Mx = pG->U[k][j][i].M3;
        U1d[k].My = pG->U[k][j][i].M1;
        U1d[k].Mz = pG->U[k][j][i].M2;
#ifndef BAROTROPIC
        U1d[k].E  = pG->U[k][j][i].E;
#endif /* BAROTROPIC */
#ifdef MHD
        U1d[k].By = pG->U[k][j][i].B1c;
        U1d[k].Bz = pG->U[k][j][i].B2c;
        Bxc[k] = pG->U[k][j][i].B3c;
        //Bxb[k] = pG->B3i[k][j][i];
        //B3_x3Face[k][j][i] = pG->B3i[k][j][i];
#endif /* MHD */




#ifdef BKG
        U1d[k].db  = pG->U[k][j][i].db;
        Bxb[k] = pG->U[k][j][i].B3cb;
        U1d[k].Byb = pG->U[k][j][i].B1cb;
        U1d[k].Bzb = pG->U[k][j][i].B2cb;
#endif



#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++) U1d[k].s[n] = pG->U[k][j][i].s[n];
#endif
      }


/*--- Step 3b ------------------------------------------------------------------
 * Compute L and R states at X3-interfaces, add "MHD source terms" for 0.5*dt
 */

      for (k=ks-nghost; k<=ke+nghost; k++) {
        W[k] = Cons1D_to_Prim1D(&U1d[k],&Bxc[k],&Bxb[k]);
      }


/*--- Step 3c ------------------------------------------------------------------
 * Add source terms from static gravitational potential for 0.5*dt to L/R states
 */

      if (StaticGravPot != NULL){
        for (k=kl+1; k<=ku; k++) {
          cc_pos(pG,i,j,k,&x1,&x2,&x3);
          phicr = (*StaticGravPot)(x1,x2, x3             );
          phicl = (*StaticGravPot)(x1,x2,(x3-    pG->dx3));
          //phifc = (*StaticGravPot)(x1,x2,(x3-0.5*pG->dx3));

          //Wl[k].Vx -= dtodx3*(phifc - phicl);
          W[k].Vx -= dtodx3*(phicr - phicl);
        }
      }


/*--- Step 3c (cont) -----------------------------------------------------------
 * Add source terms for self-gravity for 0.5*dt to L/R states
 */

#ifdef SELF_GRAVITY
      for (k=kl+1; k<=ku; k++) {
        //Wl[k].Vx -= q3*(pG->Phi[k][j][i] - pG->Phi[k-1][j][i]);
        W[k].Vx -= q3*(pG->Phi[k][j][i] - pG->Phi[k-1][j][i]);
      }
#endif


/*--- Step 3c (cont) -----------------------------------------------------------
 * Add source terms from optically-thin cooling for 0.5*dt to L/R states
 */

#ifndef BAROTROPIC
      if (CoolingFunc != NULL){
        for (k=kl+1; k<=ku; k++) {
          //coolfl = (*CoolingFunc)(Wl[k].d,Wl[k].P,(0.5*pG->dt));
          //coolfr = (*CoolingFunc)(Wr[k].d,Wr[k].P,(0.5*pG->dt));
          coolfc = (*CoolingFunc)(W[k].d+W[k].db,W[k].P,(0.5*pG->dt));

          //Wl[k].P -= 0.5*pG->dt*Gamma_1*coolfl;
          W[k].P -= 0.5*pG->dt*Gamma_1*coolfc;
        }
      }
#endif /* BAROTROPIC */



/*--- Step 3c (cont) -----------------------------------------------------------
 * Add source terms for particle feedback for 0.5*dt to L/R states
 */

#ifdef FEEDBACK
   for (k=kl+1; k<=ku; k++) {
      d1 = 0.5*(1.0/(W[k-1].d+W[k-1].db));
      W[k].Vx -= pG->Coup[k-1][j][i].fb3*d1;
      W[k].Vy -= pG->Coup[k-1][j][i].fb1*d1;
      W[k].Vz -= pG->Coup[k-1][j][i].fb2*d1;

      d1 = 0.5*(1.0/(W[k].d+W[k].db));
      W[k].Vx -= pG->Coup[k][j][i].fb3*d1;
      W[k].Vy -= pG->Coup[k][j][i].fb1*d1;
      W[k].Vz -= pG->Coup[k][j][i].fb2*d1;

#ifndef BAROTROPIC
      //Wl[i].P += pG->Coup[k-1][j][i].Eloss*Gamma_1;
      W[i].P += pG->Coup[k][j][i].Eloss*Gamma_1;
#endif
    }
#endif /* FEEDBACK */


/*--- Step 3d ------------------------------------------------------------------
 * Compute 1D fluxes in x3-direction, storing into 3D array
 */

      for (k=kl+1; k<=ku; k++) {
        Uc_x3[k][j][i] = Prim1D_to_Cons1D(&W[k],&Bxc[k], &Bxb[k]);
        //Ur_x3Face[k][j][i] = Prim1D_to_Cons1D(&Wr[k],&Bxi[k]);

//#ifdef MHD
//        Bx = B3_x3Face[k][j][i];
//#endif
        fluxes(Uc_x3[k][j][i],Uc_x3[k][j][i],W[k],W[k],Bxc,Bxb,&x3Flux[k][j][i]);

      }
    }
  }





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

#ifdef MHD
        //pG->U[k][j][i].B1c -= dtodx1*(rsf*x1Flux[k][j][i+1].Bx-lsf*x1Flux[k][j][i].Bx);
        pG->U[k][j][i].B2c -= dtodx1*(rsf*x1Flux[k][j][i+1].By-lsf*x1Flux[k][j][i].By);
        pG->U[k][j][i].B3c -= dtodx1*(rsf*x1Flux[k][j][i+1].Bz-lsf*x1Flux[k][j][i].Bz);
#endif /* MHD */



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

#ifdef MHD
        pG->U[k][j][i].B1c -= dtodx2*(rsf*x2Flux[k][j+1][i].Bz-lsf*x2Flux[k][j][i].Bz);
        //pG->U[k][j][i].B2c -= dtodx1*(rsf*x2Flux[k][j+1][i].By-lsf*x2Flux[k][j][i].By);
        pG->U[k][j][i].B3c -= dtodx2*(rsf*x2Flux[k][j+1][i].By-lsf*x2Flux[k][j][i].By);
#endif /* MHD */



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



#ifdef MHD
        pG->U[k][j][i].B1c -= dtodx3*(rsf*x3Flux[k+1][j][i].By-lsf*x3Flux[k][j][i].By);
        pG->U[k][j][i].B2c -= dtodx3*(rsf*x3Flux[k+1][j][i].Bz-lsf*x3Flux[k][j][i].Bz);
        //pG->U[k][j][i].B3c -= dtodx1*(rsf*x3Flux[k+1][j][i].Bz-lsf*x3Flux[k][j][i].Bz);
#endif /* MHD */



#if (NSCALARS > 0)
        for (n=0; n<NSCALARS; n++)
          pG->U[k][j][i].s[n] -= dtodx3*(x3Flux[k+1][j][i].s[n]
                                       - x3Flux[k  ][j][i].s[n]);
#endif
      }
    }
  }









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
}

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void integrate_init_2d(MeshS *pM)
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

/*refer to material  integrate_2d_ctu.c*/
  if ((Bxc = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;
  if ((Bxb = (Real*)malloc(nmax*sizeof(Real))) == NULL) goto on_error;


  if ((U1d= (Cons1DS*)malloc(nmax*sizeof(Cons1DS))) == NULL) goto on_error;
  if ((W  = (Prim1DS*)malloc(nmax*sizeof(Prim1DS))) == NULL) goto on_error;

  if ((Uc_x1=(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;

  if ((Uc_x2=(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;

  if ((Uc_x3=(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))== NULL) goto on_error;


  if ((x1Flux   =(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((x2Flux   =(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;
  if ((x3Flux   =(Cons1DS***)calloc_3d_array(size3,size2,size1,sizeof(Cons1DS)))
    == NULL) goto on_error;

#ifdef CYLINDRICAL
#ifndef MHD
#ifndef PARTICLES
  if((StaticGravPot != NULL) || (CoolingFunc != NULL))
#endif
#endif
#endif
  {
  if ((dhalf = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  if ((phalf = (Real***)calloc_3d_array(size3,size2,size1,sizeof(Real)))==NULL)
    goto on_error;
  }

#ifdef SHEARING_BOX
  if ((Flxiib = (ConsS**)calloc_2d_array(size3,size2,sizeof(ConsS)))==NULL)
    goto on_error;
  if ((Flxoib = (ConsS**)calloc_2d_array(size3,size2,sizeof(ConsS)))==NULL)
    goto on_error;
  if ((rFlxiib = (ConsS**)calloc_2d_array(size3,size2,sizeof(ConsS)))==NULL)
    goto on_error;
  if ((rFlxoib = (ConsS**)calloc_2d_array(size3,size2,sizeof(ConsS)))==NULL)
    goto on_error;
#endif

  /* data structures for cylindrical coordinates */
#ifdef CYLINDRICAL
  if ((geom_src = (Real***)calloc_3d_array(size3, size2, size1, sizeof(Real))) == NULL)
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
void integrate_destruct_3d(void)
{
/*refer to material  integrate_3d_ctu.c*/


  if (Bxc != NULL) free(Bxc);
  if (Bxb != NULL) free(Bxi);


  if (U1d      != NULL) free(U1d);
  if (W        != NULL) free(W);
  if (Wl       != NULL) free(Wl);
  if (Wr       != NULL) free(Wr);

  if (Uc_x1 != NULL) free_3d_array(Uc_x1);
  if (Uc_x2 != NULL) free_3d_array(Uc_x2);
  if (Uc_x3 != NULL) free_3d_array(Uc_x3);

  if (x1Flux    != NULL) free_3d_array(x1Flux);
  if (x2Flux    != NULL) free_3d_array(x2Flux);
  if (x3Flux    != NULL) free_3d_array(x3Flux);
  if (dhalf     != NULL) free_3d_array(dhalf);
  if (phalf     != NULL) free_3d_array(phalf);
#ifdef SHEARING_BOX
  if (Flxiib != NULL) free_2d_array(Flxiib);
  if (Flxoib != NULL) free_2d_array(Flxoib);
  if (rFlxiib != NULL) free_2d_array(rFlxiib);
  if (rFlxoib != NULL) free_2d_array(rFlxoib);
#endif

  /* data structures for cylindrical coordinates */
#ifdef CYLINDRICAL
  if (geom_src  != NULL) free_3d_array(geom_src);
#endif

  return;

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
