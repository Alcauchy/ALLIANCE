//
// Created by alcauchy on 03/02/2022.
//

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

extern double *space_kPerp2;
extern double *space_kPerp;
extern double *space_kx;
extern double *space_ky;
extern double *space_kz;
extern COMPLEX *space_iKx;
extern COMPLEX *space_iKy;
extern COMPLEX *space_iKz;
#endif //ALLIANCE_ALPHA_1_0_SPACE_CONFIG_H
