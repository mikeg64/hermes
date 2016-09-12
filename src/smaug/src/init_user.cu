//17/12/2012

//ozttest
__device__ __host__
void init_user_MODID (real *w, real *wd, struct params *p,int *ii) {
                    
                    real b0=1.0;
                    real ptot=5.0/3.0;
                    real rrho=25.0/36.0/PI;
                    real rgamm1;
		    real e1,e2;
		    int i,j,k;
		    i=ii[0];
		    j=ii[1];

	#ifdef USE_SAC
                    b0=1.0/sqrt(4.0*PI);
                    ptot=5.0/12.0/PI;
		    w[fencode3_i(p,ii,rhob)]=25.0/(36.0*PI);

                    rgamm1=1.0/((p->gamma)-1);
		    
		    w[fencode3_i(p,ii,b1)]=b0*sin((4.0*PI*p->dx[1])*j);
		    w[fencode3_i(p,ii,b2)]=-b0*sin(2.0*PI*(p->dx[0])*i);


		    w[fencode3_i(p,ii,mom2)]=-w[fencode3_i(p,ii,rhob)]*sin(2.0*PI*i*(p->dx[0]));
                    w[fencode3_i(p,ii,mom1)]=w[fencode3_i(p,ii,rhob)]*sin(2.0*PI*j*(p->dx[1]));

                    e1=ptot*rgamm1+(0.5*(w[fencode3_i(p,ii,mom1)]*w[fencode3_i(p,ii,mom1)]+w[fencode3_i(p,ii,mom2)]*w[fencode3_i(p,ii,mom2)])/rrho);
                    w[fencode3_i(p,ii,energy)]+=0.5*((p->gamma)-2)*(w[fencode3_i(p,ii,b1)]*w[fencode3_i(p,ii,b1)]+w[fencode3_i(p,ii,b2)]*w[fencode3_i(p,ii,b2)]);

                    e2=0.5*(w[fencode3_i(p,ii,b1)]*w[fencode3_i(p,ii,b1)]+w[fencode3_i(p,ii,b2)]*w[fencode3_i(p,ii,b2)])/rrho;
                    w[fencode3_i(p,ii,energyb)]=(e1+e2);

                    w[fencode3_i(p,ii,energy)]=w[fencode3_i(p,ii,energyb)];
                    w[fencode3_i(p,ii,energyb)]=0.0;

                    w[fencode3_i(p,ii,rho)]=w[fencode3_i(p,ii,rhob)];
                    w[fencode3_i(p,ii,rhob)]=0.0;
      #endif
}


