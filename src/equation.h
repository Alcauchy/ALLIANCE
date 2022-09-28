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

extern int *equation_forceKxInd;
extern int *equation_forceKyInd;
extern int *equation_forceKzInd;
extern int *equation_forceKxIndGlobal;
extern int *equation_forceKxIndGathered;
extern int *equation_forceKyIndGathered;
extern int *equation_forceKzIndGathered;
extern int equation_forceKn;
extern int equation_forceKnTotal;
extern int equation_forceNorm;
extern int equation_forcedM;
extern int *equation_forceKnAr;
extern int *equation_displacements;
extern double equation_forcingCoef;
extern double *equation_forcingMM;
extern COMPLEX *equation_forcingAr;
#endif //ALLIANCE_ALPHA_1_0_EQUATION_H
