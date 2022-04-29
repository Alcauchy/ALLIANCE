#ifndef ALLIANCE_ALPHA_1_0_SOLVER_H
#define ALLIANCE_ALPHA_1_0_SOLVER_H
#include "parameters_io.h"

enum solverType {RK4};
void solver_init(); // initialize solver depending on the choice.
void solver_makeStep(COMPLEX **g, COMPLEX *h);
struct solver{
    double dt;
    int Nt;
};
struct rk4{
    COMPLEX* K_buf;
    COMPLEX* RHS_buf;
    COMPLEX* g_buf;
};
extern struct solver solver;
extern struct rk4 rk4;
#endif //ALLIANCE_ALPHA_1_0_SOLVER_H
