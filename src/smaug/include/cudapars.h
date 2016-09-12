/////////////////////////////////////
// global variables and configuration section
/////////////////////////////////////

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
