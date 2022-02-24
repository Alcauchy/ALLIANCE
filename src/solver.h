//
// Created by alcauchy on 26/01/2022.
//

#ifndef ALLIANCE_ALPHA_1_0_SOLVER_H
#define ALLIANCE_ALPHA_1_0_SOLVER_H
#include "parameters_io.h"
void solver_init(); // initialize solver depending on the choice.
void solver_makeStep();
struct solver{
    double dt;
    int Nt;
};
extern struct solver solver;
#endif //ALLIANCE_ALPHA_1_0_SOLVER_H
