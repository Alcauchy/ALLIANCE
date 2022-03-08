////////////////////////////////////////////////////////////////////////////////
// 26/01/2022 created by Gene Gorbunov
//                                   PARAMETERS
//
// solver_init
// solver_makeStep
//
// VERSION 1.0
////////////////////////////////////////////////////////////////////////////////
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