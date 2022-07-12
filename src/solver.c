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
    solver.curTime = 0.;
    solver.Nt = parameters.Nt;
    solver.iter_dt = parameters.iter_dt;
    solverType = SOLVERTYPE;
    if (solverType == RK4) {
        if (mpi_my_rank == IORANK) printf("CHOSEN SOLVER IS RUNGE-KUTTA 4\n");
        rk4.K_buf = calloc(array_local_size.total_comp, sizeof(*rk4.K_buf));
        rk4.RHS_buf = calloc(array_local_size.total_comp, sizeof(*rk4.RHS_buf));
        rk4.g_buf = calloc(array_local_size.total_comp, sizeof(*rk4.g_buf));

    }
};

/***************************************
 * \fn void solver_makeStep(COMPLEX **g, COMPLEX *h):
 * \brief iterate solver forward
 * \param g: address of the 6D complex array. Modified gyrokinetic distribution function
 * \param h: 6D complex array. Gyrokinetic distribution function
 *
 * solves one simulation time step
 ***************************************/
void solver_makeStep(COMPLEX **g, COMPLEX *h, int it) {
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
            rk4.g_buf = *g;
            *g = inter;
            break;

        default:
            printf("ERROR WHILE MAKING SOLVER STEP; CHECK CHOICE OF SOLVER! ABORTING... \n");
            exit(1);
    }
};

/***************************************
 * \fn void solver_updateDt():
 * \brief updates dt
 *
 * uses CFL condition for nonlinear time step
 ***************************************/
void solver_updateDt(COMPLEX *g, COMPLEX *h, int it) {
    solver.curTime += solver.dt;
    double v_perp[2] = {0,0};
    double v_temp;
    size_t indChi;
    double *vField;
    vField = (double*) fftw_chiBuf;

    // update fields and h first
    fields_sendF(g);
    fields_getFields(f00, f10, f01);
    fields_getChi();
    distrib_getH(h, g);

    if (it % solver.iter_dt == 0) {
        switch (systemType) {
            case ELECTROSTATIC:
                // compute vx
                fields_getGradY(fftw_chiBuf);
                fftw_c2r_chi();
                for (size_t iy = 0; iy < array_local_size.ny; iy++) {
                    for (size_t ix = 0; ix < array_local_size.nx; ix++) {
                        for (size_t iz = 0; iz < array_local_size.nz + 2; iz++) {
                            for (size_t is = 0; is < array_local_size.ns; is++) {
                                indChi = getIndChiBufEL_r(ix, iy, iz, is);
                                v_temp = fabs(vField[indChi]);
                                v_perp[0] = (v_perp[0] > v_temp) ? v_perp[0] : v_temp;
                            }
                        }
                    }
                }
                // compute vy
                fields_getGradX(fftw_chiBuf);
                fftw_c2r_chi();
                for (size_t iy = 0; iy < array_local_size.ny; iy++) {
                    for (size_t ix = 0; ix < array_local_size.nx; ix++) {
                        for (size_t iz = 0; iz < array_local_size.nz + 2; iz++) {
                            for (size_t is = 0; is < array_local_size.ns; is++) {
                                indChi = getIndChiBufEL_r(ix, iy, iz, is);
                                v_temp = fabs(vField[indChi]);
                                v_perp[1] = (v_perp[1] > v_temp) ? v_perp[1] : v_temp;
                            }
                        }
                    }
                }
                break;
            case ELECTROMAGNETIC:
                // compute vx
                fields_getGradY(fftw_chiBuf);
                fftw_c2r_chi();
                for (size_t iy = 0; iy < array_local_size.ny; iy++) {
                    for (size_t ix = 0; ix < array_local_size.nx; ix++) {
                        for (size_t iz = 0; iz < array_local_size.nz + 2; iz++) {
                            for (size_t is = 0; is < array_local_size.ns; is++) {
                                for (size_t ifield = 0; ifield < 3; ifield++) {
                                    indChi = getIndChiBufEM_r(ix, iy, iz, is, ifield);
                                    v_temp = fabs(vField[indChi]);
                                    if (ifield == 2) v_temp *= sqrt((array_global_size.nm + 1.) / 2.);
                                    v_perp[0] = (v_perp[0] > v_temp) ? v_perp[0] : v_temp;
                                }
                            }
                        }
                    }
                }
                // compute vy
                fields_getGradX(fftw_chiBuf);
                fftw_c2r_chi();
                for (size_t iy = 0; iy < array_local_size.ny; iy++) {
                    for (size_t ix = 0; ix < array_local_size.nx; ix++) {
                        for (size_t iz = 0; iz < array_local_size.nz + 2; iz++) {
                            for (size_t is = 0; is < array_local_size.ns; is++) {
                                for (size_t ifield = 0; ifield < 3; ifield++) {
                                    indChi = getIndChiBufEM_r(ix, iy, iz, is, ifield);
                                    v_temp = fabs(vField[indChi]);
                                    if (ifield == 2) v_temp *= sqrt((array_global_size.nm + 1.) / 2.);
                                    v_perp[1] = (v_perp[1] > v_temp) ? v_perp[1] : v_temp;
                                }
                            }
                        }
                    }
                }
                break;
            default:
                exit(1);
        }
        MPI_Allreduce(MPI_IN_PLACE, &v_perp, 2, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        solver.dt = 0.5 / (v_perp[0] / space_dx + v_perp[1] / space_dy);
        if (mpi_my_rank == 0) printf("it = %d\t t = %f\t dt = %f\n", it, solver.curTime, solver.dt);
       // if (mpi_my_rank == 0) printf("t = %f\n", solver.curTime);
       // if (mpi_my_rank == 0) printf("dt = %f\n", solver.dt);

    }
};

/***************************************
 * \fn void solver_getLinearDt():
 * \brief computes linear dt
 *
 * computes linear dt approximation using
 * Gershgorin Discs
 ***************************************/
void solver_getLinearDt(){

}