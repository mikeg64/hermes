#include "copyright.h"
/*============================================================================*/
/*! \file solp.c
 *  \brief Problem generator for Solar physics oscillator driver.
 *
 * REFERENCE: For example, see: G. Toth,  "The div(B)=0 constraint in shock
 *   capturing MHD codes", JCP, 161, 605 (2000)				      */
/*============================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"
#include "debug_tools_cuda.h"

#ifndef MHD
#error : The problem generator solp.c only works for mhd.
#endif /* MHD */

/*----------------------------------------------------------------------------*/
/* problem:   */
#undef FIELD_LOOP
#ifdef SOLP

static Real grav_pot2(const Real x1, const Real x2);

void problem(Grid *pGrid)
{
	  //GridS *pGrid = pDomain->Grid;
	  int i,j,is,ie,js,je;
	  int nx1,nx2;
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

	  is = pGrid->is; ie = pGrid->ie;
	  js = pGrid->js; je = pGrid->je;
	  nx1 = (ie-is)+1 + 2*nghost;
	  nx2 = (je-js)+1 + 2*nghost;






	/*read VALIIc data*/
	/*see initialisation_user.h.spicule1_mpi in smaug_pmode/models*/


	    FILE *fid=fopen("test.dat","r");

	    //FILE *fid=fopen("VALMc_rho_132_test_sac_all.dat","r");
	    printf("%d %d %d %d \n",is,js,ie,je);
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
			  cc_pos(pGrid,is,i+js,&x1,&x2);
	                  printf("%f %f %f %f %f\n",x2, val3c[i][0], val3c[i][1], val3c[i][2], val3c[i][3]);
	                //if(p->ipe==1)
	                 }

	/*lagrange interpolation*/
	/*   % t1=(xval-x(i+1))/(x(i)-x(i-1));
	     % t2=(xval-x(i))/(x(i+1)-x(i));
	     % y =t1*f(i)+t2*f(i+1); */



	/* Initialize the grid.  Note the center is always assumed to have coordinates
	 * x1=0, x2=0; the grid range in the input file must be consistent with this */

	  is = pGrid->is;  ie = pGrid->ie;
	  js = pGrid->js;  je = pGrid->je;



	    for (j=js; j<=je; j++) {
	      for (i=is; i<=ie; i++) {


	        cc_pos(pGrid,i,j,&x1,&x2);
	        rho=val3c[j-js][2];
	        //en=val3c[j-js][4];
	        pres=val3c[j-js][3];
	        pGrid->U[j][i].d = rho;
	        pGrid->U[j][i].E = (pres)/Gamma_1;
	        //pGrid->U[j][i].E = en;

	        //pGrid->U[j][i].d = 1.0;
	        pGrid->U[j][i].M1 = 0.0;
	        pGrid->U[j][i].M2 = 0.0;
	        pGrid->U[j][i].M3 = 0.0;
	        pGrid->B1i[j][i] = bx0;
	        pGrid->B2i[j][i] = 0.0;
	        pGrid->U[j][i].B1c = 0.0;
	        pGrid->U[j][i].B2c = bx0;
	        pGrid->U[j][i].B3c = 0.0;







	      }
	    }



	    for (j=js; j<=je; j++) {
	      pGrid->B1i[j][ie+1] = bx0;
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

void problem_write_restart(Grid *pG, FILE *fp)
{
  return;
}

void problem_read_restart(Grid *pG, FILE *fp)
{
  return;
}

Gasfun_t get_usr_expr(const char *expr)
{
  return NULL;
}


void Userwork_in_loop(Grid *pGrid)
{
	  /*DomainS *pDomain = (DomainS*)&(pM->Domain[0][0]);
	  GridS *pGrid = pM->Domain[0][0].Grid;*/


	int i, is=pGrid->is, ie = pGrid->ie;
	  int j, js=pGrid->js, je = pGrid->je;
	  //int k, ks=pGrid->ks, ke = pGrid->ke;
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
	  //AA=0.0;
	  xcz=0.5e6;
	  xcx=2.0e6;
	  delta_z=0.064e6;
	  delta_x=0.064e6;
	  //delta_y=0.064e6;




	  if (isnan(pGrid->dt)) ath_error("Time step is NaN!");


		qt=pGrid->time;

		tdep=sin(qt*2.0*PI/s_period);
	        //tdep=1.0;



			cc_pos(pGrid,ie,je,&x1,&x2);
			xxmax=x1;
			yymax=x2;
			cc_pos(pGrid,is,js,&x1,&x2);
			xxmax=xxmax-x1;
			yymax=yymax-x2;
			xxmin=x1;
			yymin=x2;


	        /*printf("%d %d %d \n",is,js,ks);
	        printf("%d %d %d \n",ie,je,ke);
	        printf("%d %d %d \n", pGrid->Nx[0],pGrid->Nx[1], pGrid->Nx[2]);*/
		    for (j=js; j<=je; j++) {
		      for (i=is; i<=ie; i++) {
			cc_pos(pGrid,i,j,&x1,&x2);

			xp=x1-xxmin;
			zp=x2;

			r2=(zp-xcz)*(zp-xcz);
	        r1=(xp-xcx)*(xp-xcx);

			exp_z=exp(-r2/(delta_z*delta_z));
	        exp_x=exp(-r1/(delta_y*delta_y));

			//exp_xyz=sin(PI*xp*(n1+1)/xxmax)*exp_z;
			//exp_xyz=exp_y*exp_z;
	        exp_xyz=exp_x*exp_z;

			vvz=AA*exp_xyz*tdep;
	                //vvz=0;
	                //if(j==12)
	                //    printf("%d %d %d %f %f %f %f %f %f %f\n",i,j,k,xp,yp,zp,xcz,exp_x,exp_z,vvz);

	//if(i>60 && i<68)
	//if(i>is && i<ie)
	//{

	                //(j>8 && j<16 && qt<2)
	                //    printf("%d %d %d %g %g %g %g  \n",i,j,k,vvz,exp_x,exp_z,(pGrid->dt)*vvz*(pGrid->U[j][i].d));


			pGrid->U[j][i].M2 += (pGrid->dt)*vvz*(pGrid->U[j][i].d);
			pGrid->U[j][i].E += (pGrid->dt)*vvz*vvz*(pGrid->U[j][i].d)/2.0;
	//}
		      }
	              //printf("\n");

		    }






		//newtime = pGrid->time + pGrid->dt;



	  return;








}

void Userwork_after_loop(Grid *pG)
{
}


/*----------------------------------------------------------------------------*/
/*! \fn static Real grav_pot2(const Real x1, const Real x2, const Real x3)
 *  \brief Gravitational potential; g = 0.1
 */

static Real grav_pot2(const Real x1, const Real x2)
{
  return 287*x2;
  //return 0;
  //return x2;
}

#endif /* SOLP */
