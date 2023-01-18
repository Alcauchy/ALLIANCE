#ifndef ALLIANCE_ALPHA_1_0_VARIABLES_H
#define ALLIANCE_ALPHA_1_0_VARIABLES_H
#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include <math.h>
#include "fftw3-mpi.h"
#include "utils_mpi.h"
#include "array.h"
#include "utils_fftw.h"
#include "parameters_io.h"
#include "space_config.h"

void var_getJ0();
void var_getJ1();
void var_init();
void var_varInit();
size_t var_getJIndex(size_t ikx, size_t iky, size_t is);

struct phys_params{
        double beta;
        double *m;
        double *n;
        double *T;
        double *q;
        double *vT;
        double *rho;
        double k_rg;
        double *b;
        double B0;
        double mu_k;
        double mu_m;
        double mu_kz;
        double lap_k;
        double lap_kz;
        double pwr_m;
};

extern struct phys_params var_var;
extern double *var_J0;
extern double *var_J1;

#endif //ALLIANCE_ALPHA_1_0_VARIABLES_H
