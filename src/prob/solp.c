#include "copyright.h"
/*============================================================================*/
/*! \file rotor.c
 *  \brief Sets up 2D rotor test problem.
 *
 * PURPOSE: Sets up 2D rotor test problem.  The center of the grid is assumed to
 *   have coordinates (x1,x2) = [0,0]; the grid initialization must be
 *   consistent with this
 *
 * REFERENCE: G. Toth, "The div(B)=0 constraint in shock-capturing MHD codes",
 *   JCP, 161, 605 (2000)						      */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"



#ifndef MHD
#error : The rotor problem can only be run with MHD.
#endif
#ifdef ISOTHERMAL 
#error : The rotor problem can only be run with an ADIABATIC eos.
#endif


static void ry_bc(GridS *pGrid);
static Real grav_pot2(const Real x1, const Real x2, const Real x3);
static Real grav_pot3(const Real x1, const Real x2, const Real x3);


/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/
/* problem:  */

void problem(DomainS *pDomain)
{
  GridS *pGrid = pDomain->Grid;
  int i,j,k,is,ie,js,je,ks,ke;
  Real v0,p0,bx0,x1,x2,x3,rad,frac,r0,r1;

  Real amp,lx,ly,lz;
  Real a,v1,fac,w;
  double val3c[132][4],rho,pres,val;
  
  char st1[100],st2[100],st3[100],st4[100];
  int ntt;

#ifdef MHD
  Real b0,angle;
#endif


/* Read initial conditions from 'athinput' */

  v0 = par_getd("problem","v0");
  p0 = par_getd("problem","p0");
  bx0 = par_getd("problem","bx0");
  r0 = par_getd("problem","r0");
  r1 = par_getd("problem","r1");

/* Initialize the grid.  Note the center is always assumed to have coordinates
 * x1=0, x2=0; the grid range in the input file must be consistent with this */

  is = pGrid->is;  ie = pGrid->ie;
  js = pGrid->js;  je = pGrid->je;
  ks = pGrid->ks;  ke = pGrid->ke;

  lx = pDomain->RootMaxX[0] - pDomain->RootMinX[0];
  ly = pDomain->RootMaxX[1] - pDomain->RootMinX[1];
  lz = pDomain->RootMaxX[2] - pDomain->RootMinX[2];

/*read VALIIc data*/
/*see initialisation_user.h.spicule1_mpi in smaug_pmode/models*/
    
    FILE *fid=fopen("../tst/2D-mhd/VALMc_rho_132_test_sac_all.dat","r");
    printf("%d %d %d %d %d %d\n",is,js,ks,ie,je,ke);
    for(i=0; i<132; i++)
               {
                 //fscanf(fid, " %s %s %s %s %n", st1, st2, st3, st4,&ntt);
	         fscanf(fid, " %s %s %s %s", st1, st2, st3, st4);
		 //fscanf(fid, " %g %g %g %g", &val3c[131-i][0], &val3c[131-i][0], &val3c[131-i][0], &val3c[131-i][0]);
                 //printf("%s %s %s %s\n",st1,st2,st3,st4);
                 sscanf(st1,"%lf",&val); //height
                 val3c[131-i][0]=val;
		 sscanf(st2,"%lf",&val); //temp
		 val3c[131-i][1]=val;
		 sscanf(st3,"%lf",&val); //dens
		 val3c[131-i][2]=val;
		 sscanf(st4,"%lf",&val); //pres
		 val3c[131-i][3]=val;

                 

            
              }
    fclose(fid);

             for(i=0; i<132; i++)
                 if((i+js)<=je)
		{
		  cc_pos(pGrid,is,i+js,ks,&x1,&x2,&x3);
                  printf("%f %f %f %f %f\n",x2, val3c[i][0], val3c[i][1], val3c[i][2], val3c[i][3]);
                //if(p->ipe==1)
                 }

/*lagrange interpolation*/
/*   % t1=(xval-x(i+1))/(x(i)-x(i-1));
     % t2=(xval-x(i))/(x(i+1)-x(i));
     % y =t1*f(i)+t2*f(i+1); */  







/* 2D PROBLEM --------------------------------------------------------------- */
/* Initialize two fluids with interface at y=0.0.  Pressure scaled to give a
 * sound speed of 1 at the interface in the light (lower, d=1) fluid 
 * Perturb V2 using single (iprob=1) or multiple (iprob=2) mode 
 */

  if (pGrid->Nx[2] == 1) {
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
        rho=val3c[j-js][2];
        pres=val3c[j-js][3];
	pGrid->U[k][j][i].d = rho;
        pGrid->U[k][j][i].E = (pres)/Gamma_1;
	pGrid->U[k][j][i].M1 = 0.0;
        pGrid->U[k][j][i].M2 = 0.0;
        pGrid->U[k][j][i].M3 = 0.0;
        
	pGrid->U[k][j][i].E+=0.5*SQR(pGrid->U[k][j][i].M2)/pGrid->U[k][j][i].d;

        a=2.0e6;
        w=1.0e6;
        fac=exp(-(x2-a)*(x2-a)/(w*w));
        v1=(x2<(2.0e6))*fac-(x2>=(2.0e6))*fac;
        //pGrid->U[k][j][i].M1 = 1000*v1*rho;
	pGrid->U[k][j][i].M1 = 0;
          pGrid->U[k][j][i].M2 = (100*rho/4.0)*
            (sin(PI*x1/lx))*(sin(PI*x2/ly));
	pGrid->U[k][j][i].M2 = 0;

       
        //if(i==63)
        //   printf("i j v1=%d %d %f %f %f %f %f\n",i,j,x1,x2,lx,ly,(pGrid->U[k][j][i].M2)/rho);


#ifdef MHD
	pGrid->B1i[k][j][i] = b0;
	pGrid->U[k][j][i].B1c = b0;
        pGrid->U[k][j][i].E += 0.5*b0*b0;
#endif
      }
#ifdef MHD
    pGrid->B1i[k][j][ie+1] = b0;
#endif
    }
  }
}

//bvals_mhd_fun(pDomain, right_x2, ry_bc);


/* Enroll gravitational potential to give acceleration in y-direction for 2D
 * Use special boundary condition routines.  In 2D, gravity is in the
 * y-direction, so special boundary conditions needed for x2
*/


StaticGravPot = grav_pot2;

  return;
}

/*==============================================================================
 * PROBLEM USER FUNCTIONS:
 * problem_write_restart() - writes problem-specific user data to restart files
 * problem_read_restart()  - reads problem-specific user data from restart files
 * get_usr_expr()          - sets pointer to expression for special output data
 * get_usr_out_fun()       - returns a user defined output function pointer
 * get_usr_par_prop()      - returns a user defined particle selection function
 * Userwork_in_loop        - problem specific work IN     main loop
 * Userwork_after_loop     - problem specific work AFTER  main loop
 *----------------------------------------------------------------------------*/

void problem_write_restart(MeshS *pM, FILE *fp)
{
  return;
}

void problem_read_restart(MeshS *pM, FILE *fp)
{
  return;
}

ConsFun_t get_usr_expr(const char *expr)
{
  return NULL;
}

VOutFun_t get_usr_out_fun(const char *name){
  return NULL;
}

void Userwork_in_loop(MeshS *pM)
{




  DomainS *pDomain = (DomainS*)&(pM->Domain[0][0]);
  GridS *pGrid = pM->Domain[0][0].Grid;


int i, is=pGrid->is, ie = pGrid->ie;
  int j, js=pGrid->js, je = pGrid->je;
  int k, ks=pGrid->ks, ke = pGrid->ke;
  Real newtime;

  Real qt,tdep,s_period,AA;
  Real delta_x, delta_y, delta_z, xxmax, yymax, xxmin, yymin;
  Real exp_x,exp_y,exp_z,exp_xyz;
  Real r1,r2,xp, yp,zp;
  Real vvz;
  Real x1,x2,x3;

  Real xcz,xcx;

  int n1,n2;

  n1=2;
  n2=2;


  s_period=180.0; //Driver period
  AA=350.0;       //Driver amplitude
  //AA=1;
  xcz=0.5e6;
  xcx=2.0e6;
  delta_z=0.016e6;
  delta_x=0.016e6;
  delta_y=0.016e6;




  if (isnan(pGrid->dt)) ath_error("Time step is NaN!");


	qt=pGrid->time;

	tdep=sin(qt*2.0*PI/s_period);
        //tdep=1.0;


	if (pM->Nx[2] == 1)
	{
		cc_pos(pGrid,ie,je,ke,&x1,&x2,&x3);
		xxmax=x1;
		yymax=x3;
		cc_pos(pGrid,is,js,ks,&x1,&x2,&x3);
		xxmax=xxmax-x1;
		yymax=yymax-x3;
		xxmin=x1;
		yymin=x3;
	}

        /*printf("%d %d %d \n",is,js,ks);
        printf("%d %d %d \n",ie,je,ke);
        printf("%d %d %d \n", pGrid->Nx[0],pGrid->Nx[1], pGrid->Nx[2]);*/
	if (pGrid->Nx[2] == 1) {
	  for (k=ks; k<=ke; k++) {
	    for (j=js; j<=je; j++) {
	      for (i=is; i<=ie; i++) {
		cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

		xp=x1-xxmin;
		yp=x3-yymin;
		zp=x2;

		r2=(zp-xcz)*(zp-xcz);
                r1=(xp-xcx)*(xp-xcx);
		
                exp_y=exp(-r1/(delta_x*delta_x));
		exp_z=exp(-r2/(delta_z*delta_z));
                exp_x=exp(-r1/(delta_y*delta_y));

		exp_xyz=sin(PI*xp*(n1+1)/xxmax)*exp_z;
		//exp_xyz=exp_y*exp_z;
                //exp_xyz=exp_x*exp_z;

		vvz=AA*exp_xyz*tdep;
                //vvz=0;
                //if(j==12)
                //    printf("%d %d %d %f %f %f %f %f %f %f\n",i,j,k,xp,yp,zp,xcz,exp_x,exp_z,vvz);

//if(i>60 && i<68)
//if(i>is && i<ie)
//{

                if(j>8 && j<16 && qt<2)
                    printf("%d %d %d %g %g %g %g  \n",i,j,k,vvz,exp_x,exp_z,(pGrid->dt)*vvz*(pGrid->U[k][j][i].d));


		pGrid->U[k][j][i].M2 += (pGrid->dt)*vvz*(pGrid->U[k][j][i].d);
		pGrid->U[k][j][i].E += (pGrid->dt)*vvz*vvz*(pGrid->U[k][j][i].d)/2.0;
//}
	      }
              //printf("\n");

	    }
	  }
        }

	//for 3D model
	if (pM->Nx[2] > 1)
	{
		cc_pos(pGrid,ie,je,ke,&x1,&x2,&x3);
		xxmax=x1;
		yymax=x2;
		cc_pos(pGrid,is,js,ks,&x1,&x2,&x3);
		xxmax=xxmax-x1;
		yymax=yymax-x2;
		xxmin=x1;
		yymin=x2;
	}



	if (pGrid->Nx[2] > 1) {
	  for (k=ks; k<=ke; k++) {
	    for (j=js; j<=je; j++) {
	      for (i=is; i<=ie; i++) {
		cc_pos(pGrid,i,j,k,&x1,&x2,&x3);

		xp=x1-xxmin;
		yp=x2-yymin;
		zp=x3;

		r2=(x3-xcz)*(x3-xcz);
		
		exp_z=exp(-r2/(delta_z*delta_z));
		exp_xyz=sin(PI*xp*(n1+1)/xxmax)*sin(PI*yp*(n2+1)/yymax)*exp_z;

		vvz=AA*exp_xyz*tdep;
                //vvz=0;

		pGrid->U[k][j][i].M3 += (pGrid->dt)*vvz*(pGrid->U[k][j][i].d);
		pGrid->U[k][j][i].E += (pGrid->dt)*vvz*vvz*(pGrid->U[k][j][i].d)/2.0;
	      }

	    }
	  }
      }

	//newtime = pGrid->time + pGrid->dt;



  return;


}

void Userwork_after_loop(MeshS *pM)
{
}

/*----------------------------------------------------------------------------*/
/*! \fn static Real grav_pot2(const Real x1, const Real x2, const Real x3)
 *  \brief Gravitational potential; g = 0.1
 */

static Real grav_pot2(const Real x1, const Real x2, const Real x3)
{
  return 287*x2;
  //return 0;
}
/*! \fn static Real grav_pot3(const Real x1, const Real x2, const Real x3)
 *  \brief Gravitational potential; g = 0.1
 */
static Real grav_pot3(const Real x1, const Real x2, const Real x3)
{
  //return 287*x3;
  return 0;
}



/*----------------------------------------------------------------------------*/
/*! \fn static void ry_bc(GridS *pG)
 *  \brief  Apply boundary condition in right-y direction
 */

static void ry_bc(GridS *pGrid)
{
  int je = pGrid->je;
  int ks = pGrid->ks, ke = pGrid->ke, ku;
  int i,j,k,il,iu,jl,ju; /* i/j-lower/upper */

  iu = pGrid->ie + nghost;
  il = pGrid->is - nghost;

  if (pGrid->Nx[1] > 1){
    ju = pGrid->je + nghost;
    jl = pGrid->js - nghost;
  } else {
    ju = pGrid->je;
    jl = pGrid->js;
  }

  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->U[k][je+j][i]    =  pGrid->U[k][je-(j-1)][i];
        pGrid->U[k][je+j][i].M2 = 0; /* reflect 2-mom. */
        //pGrid->U[k][je+j][i].E -=
        //  pGrid->U[k][je-(j-1)][i].d*0.1*(2*j-1)*pGrid->dx2/Gamma_1;
      }
    }
  }

#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B1i[k][je+j][i] = pGrid->B1i[k][je-(j-1)][i];
      }
    }
  }

/* j=je+1 is not set for the interface field B2i */
  for (k=ks; k<=ke; k++) {
    for (j=2; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B2i[k][je+j][i] = pGrid->B2i[k][je-(j-2)][i];
      }
    }
  }

  if (pGrid->Nx[2] > 1) ku=ke+1; else ku=ke;
  for (k=ks; k<=ku; k++) {
    for (j=1; j<=nghost; j++) {
      for (i=il; i<=iu; i++) {
        pGrid->B3i[k][je+j][i] = pGrid->B3i[k][je-(j-1)][i];
      }
    }
  }
#endif

  return;
}



