#ifndef ALLIANCE_ALPHA_1_0_EQUATION_H
#define ALLIANCE_ALPHA_1_0_EQUATION_H
#include <stdio.h>
#include <stdlib.h>
#include "array.h"

void
equation_getLinearTerm(const COMPLEX *in, const COMPLEX *plus_boundary, const COMPLEX *minus_boundary, COMPLEX *out);
void equation_getNonlinearElectromagnetic();
void equation_getNonlinearElectrostatic();
void equation_getNonlinearTerm();
void equation_getRHS(const COMPLEX *in, COMPLEX *out);

#endif //ALLIANCE_ALPHA_1_0_EQUATION_H
