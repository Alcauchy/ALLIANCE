#ifndef ALLIANCE_ALPHA_1_0_SPACE_CONFIG_H
#define ALLIANCE_ALPHA_1_0_SPACE_CONFIG_H
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "fftw3-mpi.h"
#include "utils_mpi.h"
#include "array.h"
#include "utils_fftw.h"

void space_init();
void space_generateWaveSpace();     // generates wave vector space
void free_wavespace();
void space_generateMSpace();

extern double space_Lx;
extern double space_Ly;
extern double space_Lz;
extern double space_LperpMax;
extern double space_LperpMin;
extern double space_dx;
extern double space_dy;
extern double space_dz;
extern double space_dKx;
extern double space_dKy;
extern double space_dKz;
extern double space_kXmax;
extern double space_kYmax;
extern double space_kZmax;
extern double space_kPerpMax;
extern double *space_kPerp2;
extern double *space_kPerp;
extern double *space_kSq;
extern double *space_kx;
extern double *space_ky;
extern double *space_kz;
extern double *space_sqrtM;
extern double *space_zerosKx;
extern double *space_zerosKy;
extern double *space_zerosKz;
extern COMPLEX *space_iKx;
extern COMPLEX *space_iKy;
extern COMPLEX *space_iKz;
extern size_t *space_globalMIndex;


#endif //ALLIANCE_ALPHA_1_0_SPACE_CONFIG_H
