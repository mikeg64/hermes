/////////////////////////////////////
// global variables and configuration section
/////////////////////////////////////
#define MSIZE1 128
#define MSIZE2 128
#define MSIZE3 128

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

#ifdef USE_SAC
   typedef enum vars {rho, mom1, mom2, energy, b1, b2,energyb,rhob,b1b,b2b} CEV;
#endif

#ifdef USE_SAC_3D
   typedef enum vars {rho, mom1, mom2, mom3, energy, b1, b2, b3,energyb,rhob,b1b,b2b,b3b} CEV;
#endif





// problem size (vector length) N
static int N = 123456;

// number of threads per block

//tnumThreadsPerBlock  used to control the shared memory size 
//and nume threads for global maxima routines


//parameters used for fermi
//static int numThreadsPerBlock = 512;
//static int tnumThreadsPerBlock = 128;


//parameters used for kepler
static int numThreadsPerBlock = 64;
static int tnumThreadsPerBlock = 64;
//static int tnumThreadsPerBlock = 128;
// device to use in case there is more than one
static int selectedDevice = 0;
static int blocksize = 512;

static int dimblock = 16;
