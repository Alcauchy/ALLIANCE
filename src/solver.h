#ifndef ALLIANCE_ALPHA_1_0_SOLVER_H
#define ALLIANCE_ALPHA_1_0_SOLVER_H
#include "parameters_io.h"

enum solverType {RK4};
void solver_init(); // initialize solver depending on the choice.
void solver_makeStep(COMPLEX *g);
struct solver{
    double dt;
    int Nt;
};
struct rk4{
    COMPLEX* K1;
    COMPLEX* K2;
    COMPLEX* K3;
    COMPLEX* K4;
    COMPLEX* buf;
};
extern struct solver solver;
extern struct rk4 rk4;
#endif //ALLIANCE_ALPHA_1_0_SOLVER_H
