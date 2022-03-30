////////////////////////////////////////////////////////////////////////////////
// 26/01/2022 created by Gene Gorbunov
//                                   SOLVER
//
// solver_init
// solver_makeStep
//
// VERSION 1.0
////////////////////////////////////////////////////////////////////////////////
#include "solver.h"

#define SOLVERTYPE RK4
#define IORANK 0

enum solverType solverType;
struct solver solver;
struct rk4 rk4;

/***************************************
 * solver_init():
 ***************************************/
void solver_init() {
    solver.dt = parameters.dt;
    solver.dt = 0.01;
    solver.Nt = parameters.Nt;
    solverType = SOLVERTYPE;
    if (solverType == RK4) {
        if (mpi_my_rank == IORANK) printf("CHOSEN SOLVER IS RUNGE-KUTTA 4\n");
        rk4.K1 = malloc(array_local_size.total_comp * sizeof(*rk4.K1));
        rk4.K2 = malloc(array_local_size.total_comp * sizeof(*rk4.K2));
        rk4.K3 = malloc(array_local_size.total_comp * sizeof(*rk4.K3));
        rk4.K4 = malloc(array_local_size.total_comp * sizeof(*rk4.K4));
        rk4.buf = malloc(array_local_size.total_comp * sizeof(*rk4.buf));
        for (size_t i = 0; i < array_local_size.total_comp; i++){
            rk4.K1[i] = 0;
            rk4.K2[i] = 0;
            rk4.K3[i] = 0;
            rk4.K4[i] = 0;
            rk4.buf[i] = 0;
        }
    }
};

/***************************************
 * solver_makeStep(COMPLEX *g):
 ***************************************/
void solver_makeStep(COMPLEX *g) {
    switch (solverType) {
        case RK4:
            equation_getRHS(g, rk4.K1);
            for (size_t i = 0; i < array_local_size.total_comp; i++) {
                rk4.buf[i] = g[i] + 0.5 * solver.dt * rk4.K1[i];
            }
            equation_getRHS(rk4.buf, rk4.K2);
            for (size_t i = 0; i < array_local_size.total_comp; i++) {
                rk4.buf[i] = g[i] + 0.5 * solver.dt * rk4.K2[i];
            }
            equation_getRHS(rk4.buf, rk4.K3);
            for (size_t i = 0; i < array_local_size.total_comp; i++) {
                rk4.buf[i] = g[i] + solver.dt * rk4.K3[i];
            }
            equation_getRHS(rk4.buf, rk4.K4);
            for (size_t i = 0; i < array_local_size.total_comp; i++) {
                g[i] += 1. / 6. * solver.dt * (rk4.K1[i] + 2. * rk4.K2[i] + 2. * rk4.K3[i] + rk4.K4[i]);
            }
            break;

        default:
            printf("ERROR WHILE MAKING SOLVER STEP; CHECK CHOICE OF SOLVER! ABORTING... \n");
            exit(1);
    }
};