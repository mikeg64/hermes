#ifndef TYPES_H_
#define TYPES_H_

#define DEFINE_PRECISION(T) \
  typedef T real;


#ifdef USE_REAL
DEFINE_PRECISION(float)
#else
DEFINE_PRECISION(double)
#endif




#undef DEFINE_PRECISION


#define NDIM 2
#define NVECDIM 3
#ifdef USE_SAC
   //#define NVAR 13
   #define NVAR 10
   //#define NDERV 19
   #define NDERV 19
   #define NTEMP 8
   #define NTEMP1 2
   #define NTEMP2 1
   #define NDIM 2
   #define NVECDIM 2
#endif
#ifdef USE_SAC_3D
   #define NVAR 13
   //#define NVAR 10
   #define NDERV 23
   //#define NDERV 15
   #define NTEMP 8
   #define NTEMP1 2
   #define NTEMP2 1
   #define NDIM 3
   #define NVECDIM 3
#endif


 #ifdef USE_VAC
   //#define NVAR 8
   #define NVAR 6
   //#define NDERV 17
   #define NDERV 13
   #define NTEMP 8
   #define NTEMP1 2
   #define NTEMP2 1
   #define NDIM 2
   #define NVECDIM 2
 #endif
 #ifdef ADIABHYDRO
   //#define NVAR 4
   #define NVAR 4
   //#define NDERV 9
   #define NDERV 5
   #define NTEMP 1
   #define NDIM 2
   #define NVECDIM 3
 #endif

#define PI 3.14159265358979
struct Meta {
   char *directory ;
   char *author;
   char *sdate;
   char *platform;
   char *desc;
   char *name;
   char *ini_file;
   char *log_file;
   char *out_file;
};

struct Iome {
    char *server;
    int port;
    int id;
};

struct params {
	int n[NDIM];
	int ng[NDIM];
        int npgp[NDIM];
        int fullgridini;

        real hdmean;
        real hdmax;
        real pmax;
        real pmin;
        real pmean;
        real emax;
        real emin;
        real emean;
        

        real xmax[NDIM];
        real xmin[NDIM];
	int nt;
        int it;
        real tmax;
        real qt;

        real boundu[NDIM][NVAR];
        real boundl[NDIM][NVAR];

        real cmax;
        real maxcourant;
        int steeringenabled;
        int finishsteering;     
	real dt;
        real dx[NDIM];



        real gamma;
/*constant used for adiabatic hydrodynamics*/
         #ifdef ADIABHYDRO
            real adiab;
        #endif
        real mu;
        real eta;
        real g[NDIM];

	int sodifon;
        int rkon;
        int moddton;
        int divbon;
        int divbfix;
        int cfgsavefrequency;
        int hyperdifmom; 
        int mode;
        int readini;
        real courant;
        real maxviscoef;
        real dtdiffvisc;
        real chyp[NVAR];
        real chyp3;
        real test;  
        int boundtype[NVAR][NDIM][2];  //boundtype=0 is periodic 1=mpi 2=mpiperiod 3=cont contcd4=4  fixed=5 symm=6 asymm=7
        
        int gpid[16];
        int npe;
        int noghost;
       #ifdef USE_MULTIGPU
		int ipe;
	
                int pnpe[NDIM];
                int pipe[NDIM];
                int mpiupperb[NDIM];
                int mpilowerb[NDIM];

                int gpudirectgroup;
                int ngpudirectgroups;

                int gpudirectgroupneighb[2][NDIM]; //gpudirect group ID for each neighbour
                
                //nearest neighbours                
                int phpe[NDIM];
                int pjpe[NDIM];
                int hpe;
                int jpe;
 		#ifdef USE_SAC_3D
                	int pkpe[NDIM];
			int kpe;
		#endif               
                int gpemin[NDIM];
                int gpemax[NDIM];

                //global value of box dimensions
		real gxmax[NDIM];
		real gxmin[NDIM];
       #endif   
};

//it   t   dt    rho m1 m2 e bx by
struct state{
	int it;
	real t;
	real dt;
	real rho;
        real m1;
	real m2;
	real m3;
	real e;
        real b1;
        real b2;
        real b3;
};

struct hydrovars{
    int numvars; //variables each vector component
	int num;   //total number of dimensions including any ghost variables
	real *w;

};

/*         #ifdef USE_SAC

         #else

         #endif*/


typedef enum mode {run,scatter,gather,init,redistribute} MODE;

//typedef enum oldvars {mom3, b3,b3b} CEVOLD;
#ifdef USE_SAC
   //typedef enum vars {rho, mom1, mom2, mom3, energy, b1, b2, b3,rhob,energyb,b1b,b2b,b3b} CEV;
   typedef enum vars {rho, mom1, mom2, energy, b1, b2,energyb,rhob,b1b,b2b} CEV;
#endif

#ifdef USE_SAC_3D
   typedef enum vars {rho, mom1, mom2, mom3, energy, b1, b2, b3,energyb,rhob,b1b,b2b,b3b} CEV;
   //typedef enum vars {rho, mom1, mom2, energy, b1, b2,energyb,rhob,b1b,b2b} CEV;
#endif

#ifdef ADIABHYDRO
   //typedef enum vars {rho, mom1, mom2, mom3, energy, b1, b2, b3} CEV;
   typedef enum vars {rho, mom1, mom2, energy, b1, b2} CEV;
#endif

#ifdef USE_SAC
	typedef enum dvars {vel1,vel2,flux,hdnur,hdnul,nushk1,nushk2,soundspeed,pressuret,pressurek,bdotv,divb,cfast,ptb,pkb,pos1,pos2,delx1,delx2} DEV;

//typedef enum dvars {vel1,vel2,soundspeed,pressuret,pressurek,current1,current2,bdotv,divb,cfast,hdnur,hdnul,ptb,pkb} DEV;
#endif
#ifdef USE_SAC_3D
	typedef enum dvars {vel1,vel2,vel3,flux,hdnur,hdnul,nushk1,nushk2,nushk3,soundspeed,pressuret,pressurek,bdotv,divb,cfast,ptb,pkb,pos1,pos2,pos3,delx1,delx2,delx3} DEV;
//typedef enum dvars {vel1,vel2,soundspeed,pressuret,pressurek,current1,current2,bdotv,divb,cfast,hdnur,hdnul,ptb,pkb} DEV;
#endif
#ifdef ADIABHYDRO
	typedef enum dvars {vel1,vel2,flux,soundspeed,pressuret,pressurek,current1,current2,bdotv,divb,cfast,hdnur,hdnul} DEV;
//typedef enum dvars {vel1,vel2,soundspeed,pressuret,pressurek,current1,current2,bdotv,divb,cfast,hdnur,hdnul} DEV;
#endif


//typedef enum tempvars {tmp1, tmp2, tmp3,tmp4,tmp5,tmp6, tmprhol, tmprhor } TEV;
typedef enum tempvars {tmp1, tmp2, tmp3,tmp4,tmp5,tmp6,tmp7,tmp8 } TEV;
typedef enum temp1vars {d1,d3 } TEV1;
typedef enum temp2vars {tmpnui,tmpnui1,tmpnui2 } TEV2; //note tmpnui1 and tmpnui2 not used on GPU



typedef struct Source source;
typedef struct Constants constants;
typedef struct Domain domain;
typedef struct Iome iome;
typedef struct Meta meta;
typedef struct Stateinfo stateinfo;
typedef struct params Params;
#endif

