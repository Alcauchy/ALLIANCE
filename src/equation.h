#ifndef ALLIANCE_ALPHA_1_0_EQUATION_H
#define ALLIANCE_ALPHA_1_0_EQUATION_H
#include <stdio.h>
#include <stdlib.h>
#include "array.h"

void equation_getLinearTerm(const COMPLEX *in, const COMPLEX *plus_boundary, const COMPLEX *minus_boundary, COMPLEX *out);
void equation_getDissipation(const COMPLEX *h, COMPLEX *rhs);
void equation_getNonlinearElectromagnetic();
void equation_getNonlinearElectrostatic();
void equation_getNonlinearProduct(double *in, double *chiAr, double *out, double sign);
void equation_getRHS(const COMPLEX *in_g, COMPLEX *in_h, COMPLEX *out);
void equation_getNonlinearTerm(const COMPLEX *h, COMPLEX *out);
void equation_init();
void equation_getForcing(const COMPLEX *h, COMPLEX *rhs);

#endif //ALLIANCE_ALPHA_1_0_EQUATION_H
