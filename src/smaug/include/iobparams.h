#ifndef BPARAMS_H_
#define BPARAMS_H_

#define MSIZE1 128
#define MSIZE2 128
#define MSIZE3 128


struct bparams {


        /*#ifdef USE_SAC_3D
          real fixed1[4*NDIM*NDIM*NVAR];
          real fixed2[4*NDIM*NDIM*NVAR];
          real fixed3[4*NDIM*NDIM*NVAR];
        #else
          real fixed1[4*NDIM*NVAR];
          real fixed2[4*NDIM*NVAR];
        #endif*/


        #ifdef USE_SAC_3D
          real fixed1[4*MSIZE2*MSIZE3*NVAR];
          real fixed2[4*MSIZE1*MSIZE3*NVAR];
          real fixed3[4*MSIZE1*MSIZE2*NVAR];
        #else
          real fixed1[4*MSIZE2*NVAR];
          real fixed2[4*MSIZE1*NVAR];
        #endif

 
};


typedef struct bparams Bparams;
#endif

