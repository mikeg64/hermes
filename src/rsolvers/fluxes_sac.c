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



#ifdef SAC_INTEGRATOR
//#ifdef BKG


/* Test the intermediate states in the approximate Riemann solution. */
//#define TEST_INTERMEDIATE_STATES

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

//#ifdef BKG
void fluxes(const Cons1DS Ul, const Cons1DS Ur,
            const Prim1DS Wl, const Prim1DS Wr,
            const Real Bxi, const Real Bxib, Cons1DS *pFlux)
//#else
//void fluxes(const Cons1DS Ul, const Cons1DS Ur,
//            const Prim1DS Wl, const Prim1DS Wr,
//            const Real Bxi, Cons1DS *pFlux)
//#endif
{
  Real pbl=0.0;
  Prim1DS W;
  Cons1DS Fc;   /*flux at cell centre the l values are passed in as values at cell centres*/

/*--- Step 1. ------------------------------------------------------------------
 * Convert left-  states in conserved to primitive variables.
 */

//#ifdef BKG
  W = Cons1D_to_Prim1D(&Ul,&Bxi, &Bxib);
//#else
//  pbl = Cons1D_to_Prim1D(&Ul,&Wl,&Bxi);
//#endif


/*--- Step 2. ------------------------------------------------------------------
 * Compute L fluxes 
 */


  Fc.d  = Ul.Mx;   /*computed using (rho+rhob)*velocity */
  Fc.Mx = Ul.Mx*Wl.Vx;
  Fc.My = Ul.Mx*Wl.Vy;
  Fc.Mz = Ul.Mx*Wl.Vz;

#ifdef ISOTHERMAL
  Fc.Mx += (Wl.d+Wl.db)*Iso_csound2;
#else
  Fc.Mx += Wl.P;
  Fc.E  = (Ul.E + Ul.Eb+ Wl.P)*Wl.Vx;
#endif /* ISOTHERMAL */

#ifdef MHD
  Fc.Mx += 0.5*(Bxi*Bxi + SQR(Wl.By) + SQR(Wl.Bz))+(Bxi*Bxib+Wl.By*Wl.Byb+Wl.Bz*Wl.Bzb);/*thermal pressure plus mag pressure time*/
  Fc.Mx += -(Bxi*Bxi + SQR(Wl.By) + SQR(Wl.Bz));
  Fc.Mx -= (Bxi*Bxib+Bxi*Wl.Byb+Bxi*Wl.Bzb)+(Bxib*Bxi+Bxib*Wl.By+Bxib*Wl.Bz);
  Fc.My -= (Wl.By*Bxib+Wl.By*Wl.Bzb+Bxi*Wl.Byb+Wl.Bz*Wl.Byb )-(Wl.By*Bxi+Wl.By*Wl.Bz);
  Fc.Mz -= (Wl.Bz*Bxib+Wl.Bz*Wl.Byb+Bxi*Wl.Bzb+Wl.By*Wl.Bzb )-(Wl.Bz*Bxi+Wl.Bz*Wl.By);

#ifndef ISOTHERMAL
  Fc.E += (pbl*Wl.Vx - Bxi*(Bxi*Wl.Vx + Wl.By*Wl.Vy + Wl.Bz*Wl.Vz));
#endif /* ISOTHERMAL */

  Fc.By = Wl.By*Wl.Vx - Bxi*Wl.Vy;
  Fc.Bz = Wl.Bz*Wl.Vx - Bxi*Wl.Vz;
#endif /* MHD */

/* Fluxes of passively advected scalars, computed from density flux */
#if (NSCALARS > 0)
  for (n=0; n<NSCALARS; n++) {
    Fc.s[n] = Fc.d*Wl.r[n];

  }
#endif




	*pFlux = Fc;



  return;
}

//#endif
#endif /* SAC_FLUX */
