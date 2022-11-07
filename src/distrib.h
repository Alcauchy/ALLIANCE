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
void distrib_enforceReality(COMPLEX *f);                   // enforces reality condition
void distrib_setZeroNHalf(COMPLEX *f);                     // sets all f(n/2) = 0 since there is no way to distinguish derivatives for them
void dealiasing23(COMPLEX *data_c);                        // performs dealiasing by employing 2/3 rule. Turns all the aliased modes to zero.
void distrib_dealiasing(COMPLEX *data_c);
#endif //ALLIANCE_ALPHA_1_0_DISTRIB_H
