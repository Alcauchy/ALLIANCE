/**************************************
* @file solver.c
* \brief numerical solver
*
***************************************/
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
 * \fn void solver_init():
 * \brief initializes solver
 *
 * initializes solver with the <tt>solverType</tt>.
 ***************************************/
void solver_init() {
    solver.dt = parameters.dt;
    //solver.dt = 0.015;
    solver.Nt = parameters.Nt;
    solverType = SOLVERTYPE;
    if (solverType == RK4) {
        if (mpi_my_rank == IORANK) printf("CHOSEN SOLVER IS RUNGE-KUTTA 4\n");
        rk4.K_buf = calloc(array_local_size.total_comp, sizeof(*rk4.K_buf));
        rk4.RHS_buf = calloc(array_local_size.total_comp, sizeof(*rk4.RHS_buf));
        rk4.g_buf = calloc(array_local_size.total_comp, sizeof(*rk4.g_buf));

    }
};

/***************************************
 * \fn void solver_makeStep(COMPLEX *g):
 * \brief iterate solver forward
 *
 * solves one simulation time step
 ***************************************/
void solver_makeStep(COMPLEX **g, COMPLEX *h) {
    COMPLEX *g_ar = *g;
    switch (solverType) {
        case RK4:
            //computing k1
            equation_getRHS(g_ar, h, rk4.K_buf);
            // setting additional constrains on substeps of rk4:
            distrib_enforceReality(rk4.K_buf);
            distrib_setZeroNHalf(rk4.K_buf);
            for (size_t i = 0; i < array_local_size.total_comp; i++) {
                rk4.RHS_buf[i] = g_ar[i] + 0.5 * solver.dt * rk4.K_buf[i];
                rk4.g_buf[i] = g_ar[i] + solver.dt/6. * rk4.K_buf[i];
                rk4.K_buf[i] = 0.j;
            }
            // computing k2
            equation_getRHS(rk4.RHS_buf, h, rk4.K_buf);
            // setting additional constrains on substeps of rk4:
            distrib_enforceReality(rk4.K_buf);
            distrib_setZeroNHalf(rk4.K_buf);
            for (size_t i = 0; i < array_local_size.total_comp; i++) {
                rk4.RHS_buf[i] = g_ar[i] + 0.5 * solver.dt * rk4.K_buf[i];
                rk4.g_buf[i] += solver.dt/3. * rk4.K_buf[i];
                rk4.K_buf[i] = 0.j;
            }
            //computing k3
            equation_getRHS(rk4.RHS_buf, h, rk4.K_buf);
            // setting additional constrains on substeps of rk4:
            distrib_enforceReality(rk4.K_buf);
            distrib_setZeroNHalf(rk4.K_buf);
            for (size_t i = 0; i < array_local_size.total_comp; i++) {
                rk4.RHS_buf[i] = g_ar[i] + solver.dt * rk4.K_buf[i];
                rk4.g_buf[i] += solver.dt/3. * rk4.K_buf[i];
                rk4.K_buf[i] = 0.j;
            }
            //computing k4
            equation_getRHS(rk4.RHS_buf, h, rk4.K_buf);
            // setting additional constrains on substeps of rk4:
            distrib_enforceReality(rk4.K_buf);
            distrib_setZeroNHalf(rk4.K_buf);
            for (size_t i = 0; i < array_local_size.total_comp; i++) {
                rk4.g_buf[i] += solver.dt/6. * rk4.K_buf[i];
                rk4.K_buf[i] = 0.j;
            }
            COMPLEX *inter = rk4.g_buf;
            //printf("3rd = %p\n", inter);
            rk4.g_buf = *g;
            *g = inter;
            break;

        default:
            printf("ERROR WHILE MAKING SOLVER STEP; CHECK CHOICE OF SOLVER! ABORTING... \n");
            exit(1);
    }
};