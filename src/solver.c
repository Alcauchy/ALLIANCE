//
// Created by alcauchy on 26/01/2022.
//

#include "solver.h"


struct solver solver;

/***************************************
 * solver_init():
 ***************************************/
void solver_init(){
    solver.dt = parameters.dt;
    solver.Nt = parameters.Nt;
};

/***************************************
 * solver_makeStep():
 ***************************************/
void solver_makeStep(){};