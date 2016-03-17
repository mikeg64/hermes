#include "../copyright.h"
/*============================================================================*/
/*! \file roe.c
 *  \brief Computes 1D fluxes for sac.
 *
 *
 *
 * CONTAINS PUBLIC FUNCTIONS: 
 * - fluxes() - all Riemann solvers in Athena must have this function name and
 *              use the same argument list as defined in rsolvers/prototypes.h
 */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

/* maximum wavespeed used by H-correction, value passed from integrator */
Real etah=0.0;

#ifdef ROE_FLUX


/* Test the intermediate states in the approximate Riemann solution. */
#define TEST_INTERMEDIATE_STATES

/*----------------------------------------------------------------------------*/
/*! \fn void fluxes(const Cons1DS Ul, const Cons1DS Ur,
 *            const Prim1DS Wl, const Prim1DS Wr,
 *            const Real Bxi, Cons1DS *pFlux)
 *  \brief Computes 1D fluxes
 *   Input Arguments:
 *   - Bxi = B in direction of slice at cell interface
 *   - Ul,Ur = L/R-states of CONSERVED variables at cell interface
 *   Output Arguments:
 *   - pFlux = pointer to fluxes of CONSERVED variables at cell interface
 */

void fluxes(const Cons1DS Ul, const Cons1DS Ur,
            const Prim1DS Wl, const Prim1DS Wr,
            const Real Bxi, Cons1DS *pFlux)
{
  Real sqrtdl,sqrtdr,isdlpdr,droe,v1roe,v2roe,v3roe,pbl=0.0,pbr=0.0;
#ifndef ISOTHERMAL
  Real hroe;
#endif /* ISOTHERMAL */
#ifdef MHD
  Real b2roe,b3roe,x,y;
#endif
  Real coeff[NWAVE];
  Real ev[NWAVE],rem[NWAVE][NWAVE],lem[NWAVE][NWAVE];
  Real dU[NWAVE],a[NWAVE];
#ifdef TEST_INTERMEDIATE_STATES
  Real u_inter[NWAVE];
#ifdef ADIABATIC
  Real p_inter=0.0;
#endif
#endif /* TEST_INTERMEDIATE_STATES */
  Real *pUl, *pUr, *pFl, *pFr, *pF;
  Cons1DS Fl,Fr;
  int n,m,hlle_flag;

/*--- Step 1. ------------------------------------------------------------------
 * Convert left- and right- states in conserved to primitive variables.
 */

/*
  pbl = Cons1D_to_Prim1D(&Ul,&Wl,&Bxi);
  pbr = Cons1D_to_Prim1D(&Ur,&Wr,&Bxi);
*/

/*--- Step 2. ------------------------------------------------------------------
 * Compute Roe-averaged data from left- and right-states
 */

  sqrtdl = sqrt((double)Wl.d);
  sqrtdr = sqrt((double)Wr.d);
  isdlpdr = 1.0/(sqrtdl + sqrtdr);

  droe  = sqrtdl*sqrtdr;
  v1roe = (sqrtdl*Wl.Vx + sqrtdr*Wr.Vx)*isdlpdr;
  v2roe = (sqrtdl*Wl.Vy + sqrtdr*Wr.Vy)*isdlpdr;
  v3roe = (sqrtdl*Wl.Vz + sqrtdr*Wr.Vz)*isdlpdr;

/* The Roe average of the magnetic field is defined differently  */

#ifdef MHD
  b2roe = (sqrtdr*Wl.By + sqrtdl*Wr.By)*isdlpdr;
  b3roe = (sqrtdr*Wl.Bz + sqrtdl*Wr.Bz)*isdlpdr;
  x = 0.5*(SQR(Wl.By - Wr.By) + SQR(Wl.Bz - Wr.Bz))/(SQR(sqrtdl + sqrtdr));
  y = 0.5*(Wl.d + Wr.d)/droe;
  pbl = 0.5*(SQR(Bxi) + SQR(Wl.By) + SQR(Wl.Bz));
  pbr = 0.5*(SQR(Bxi) + SQR(Wr.By) + SQR(Wr.Bz));
#endif

/*
 * Following Roe(1981), the enthalpy H=(E+P)/d is averaged for adiabatic flows,
 * rather than E or P directly.  sqrtdl*hl = sqrtdl*(el+pl)/dl = (el+pl)/sqrtdl
 */


/*--- Step 3. ------------------------------------------------------------------
 * Compute eigenvalues and eigenmatrices using Roe-averaged values
 */


/*--- Step 4. ------------------------------------------------------------------
 * Compute L/R fluxes 
 */

  Fl.d  = Ul.Mx;
  Fr.d  = Ur.Mx;

  Fl.Mx = Ul.Mx*Wl.Vx;
  Fr.Mx = Ur.Mx*Wr.Vx;

  Fl.My = Ul.Mx*Wl.Vy;
  Fr.My = Ur.Mx*Wr.Vy;

  Fl.Mz = Ul.Mx*Wl.Vz;
  Fr.Mz = Ur.Mx*Wr.Vz;

#ifdef ISOTHERMAL
  Fl.Mx += Wl.d*Iso_csound2;
  Fr.Mx += Wr.d*Iso_csound2;
#else
  Fl.Mx += Wl.P;
  Fr.Mx += Wr.P;

  Fl.E  = (Ul.E + Wl.P)*Wl.Vx;
  Fr.E  = (Ur.E + Wr.P)*Wr.Vx;
#endif /* ISOTHERMAL */

#ifdef MHD
  Fl.Mx -= 0.5*(Bxi*Bxi - SQR(Wl.By) - SQR(Wl.Bz));
  Fr.Mx -= 0.5*(Bxi*Bxi - SQR(Wr.By) - SQR(Wr.Bz));

  Fl.My -= Bxi*Wl.By;
  Fr.My -= Bxi*Wr.By;
    
  Fl.Mz -= Bxi*Wl.Bz;
  Fr.Mz -= Bxi*Wr.Bz;

#ifndef ISOTHERMAL
  Fl.E += (pbl*Wl.Vx - Bxi*(Bxi*Wl.Vx + Wl.By*Wl.Vy + Wl.Bz*Wl.Vz));
  Fr.E += (pbr*Wr.Vx - Bxi*(Bxi*Wr.Vx + Wr.By*Wr.Vy + Wr.Bz*Wr.Vz));
#endif /* ISOTHERMAL */

  Fl.By = Wl.By*Wl.Vx - Bxi*Wl.Vy;
  Fr.By = Wr.By*Wr.Vx - Bxi*Wr.Vy;

  Fl.Bz = Wl.Bz*Wl.Vx - Bxi*Wl.Vz;
  Fr.Bz = Wr.Bz*Wr.Vx - Bxi*Wr.Vz;
#endif /* MHD */

#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) {
    Fl.s[n] = Fl.d*Wl.r[n];
    Fr.s[n] = Fr.d*Wr.r[n];
  }
#endif

/*--- Step 5. ------------------------------------------------------------------
 * Return upwind flux if flow is supersonic
 */


/*--- Step 6. ------------------------------------------------------------------
 * Compute projection of dU onto L eigenvectors ("vector A")
 */


/*--- Step 7. ------------------------------------------------------------------
 * Check that the density and pressure in the intermediate states are positive.
 * If not, set hlle_flag=1 if d_inter<0; hlle_flag=2 if p_inter<0, get HLLE
 * fluxes, and return
 */

 

#endif /* TEST_INTERMEDIATE_STATES */

/*--- Step 8. ------------------------------------------------------------------
 * Compute Roe flux */


/* Fluxes of passively advected scalars, computed from density flux */


  return;
}


#endif /* ROE_FLUX */
