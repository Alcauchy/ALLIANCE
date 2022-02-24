//
// Created by alcauchy on 03/02/2022.
//

#ifndef ALLIANCE_ALPHA_1_0_DISTRIB_H
#define ALLIANCE_ALPHA_1_0_DISTRIB_H
#include <stdlib.h>
#include <stdio.h>
#include "parameters_io.h"
#include "array.h"
#include "utils_mpi.h"
#include "fftw_utils.h"
#include "utils_hdf.h"
#include "fields.h"
#include "equation.h"
#include "space_config.h"
#include "solver.h"
#include "variables.h"

void distrib_getH(COMPLEX *h, const COMPLEX *g);
void distrib_getG(COMPLEX *g, const COMPLEX *h);
#endif //ALLIANCE_ALPHA_1_0_DISTRIB_H
