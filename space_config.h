//
// Created by alcauchy on 03/02/2022.
//

#ifndef ALLIANCE_ALPHA_1_0_SPACE_CONFIG_H
#define ALLIANCE_ALPHA_1_0_SPACE_CONFIG_H
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "array.h"
#include "mpi_utils.h"
#include "fftw_utils.h"

void space_generateWaveSpace();     // generates wave vector space
void space_xGrad(COMPLEX *data);    // computes gradient in x direction
void space_yGrad(COMPLEX *data);    // computes gradient in y direction
void space_zGrad(COMPLEX *data);    // computes gradient in z direction
#endif //ALLIANCE_ALPHA_1_0_SPACE_CONFIG_H
