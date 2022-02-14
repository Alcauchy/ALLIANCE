//
// Created by alcauchy on 26/01/2022.
//

#include "solver.h"


struct solver solver;

void solver_init(){
    solver.dt = parameters.dt;
    solver.Nt = parameters.Nt;
};
void solver_makeStep(){};