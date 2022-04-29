#ifndef ALLIANCE_ALPHA_1_0_DISTRIB_H
#define ALLIANCE_ALPHA_1_0_DISTRIB_H
#include <stdlib.h>
#include <stdio.h>
#include "parameters_io.h"
#include "array.h"
#include "utils_mpi.h"
#include "utils_fftw.h"
#include "utils_hdf.h"
#include "fields.h"
#include "equation.h"
#include "space_config.h"
#include "solver.h"
#include "variables.h"

void distrib_getH(COMPLEX *h, const COMPLEX *g);
void distrib_getG(COMPLEX *g, const COMPLEX *h);
void distrib_getXGrad(const COMPLEX *in, COMPLEX *out);    // computes gradient in x direction
void distrib_getYGrad(const COMPLEX *in, COMPLEX *out);    // computes gradient in y direction
void distrib_getZGrad(const COMPLEX *in, COMPLEX *out);    // computes gradient in z direction
#endif //ALLIANCE_ALPHA_1_0_DISTRIB_H
