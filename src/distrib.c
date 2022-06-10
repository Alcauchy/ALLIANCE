/**************************************
* @file distrib.c
*   \brief gyrokinetic distribution function module
*
*   everything required to perform different manipulations to distribution functions
***************************************/
////////////////////////////////////////////////////////////////////////////////
// 03/02/2022 created by Gene Gorbunov
//                                   DISTRIBUTION
//
// distrib_getH
// distrib_getG
// distrib_getXGrad
// distrib_getYGrad
// distrib_getZGrad
//
// VERSION 1.0
////////////////////////////////////////////////////////////////////////////////
#include "distrib.h"

/***************************************
 * \fn void distrib_getH(COMPLEX *h, const COMPLEX *g)
 * \brief computes h from g
 * \param h: complex array to store h
 * \param g: complex array with g
 *
 * computes gyrokinetic distribution function h
 * from modified gyrokinetic distribution function g.
 * Please note that before calling this function
 * gyrokinetic potentials must be computed
 ***************************************/
void distrib_getH(COMPLEX *h, const COMPLEX *g) {
    size_t ind6D;
    size_t ind4D;
    for (size_t i = 0; i < array_local_size.total_comp; i++) {
        h[i] = g[i];
    }
    switch (systemType) {
        case ELECTROSTATIC:
            break;
        case ELECTROMAGNETIC:
            for (size_t ix = 0; ix < array_local_size.nkx; ix++) {
                for (size_t iy = 0; iy < array_local_size.nky; iy++) {
                    for (size_t iz = 0; iz < array_local_size.nkz; iz++) {
                        for (size_t is = 0; is < array_local_size.ns; is++) {
                            ind4D = ix * array_local_size.nky * array_local_size.nkz * array_local_size.ns +
                                    iy * array_local_size.nkz * array_local_size.ns +
                                    iz * array_local_size.ns +
                                    is;
                            for (size_t im = 0; im < array_local_size.nm; im++) {
                                if (global_nm_index[im] == 0) {
                                    ind6D = get_flat_c(is, 0, 0, ix, iy, iz);
                                    h[ind6D] += var_var.q[is] / var_var.T[is] *
                                                (fields_chi.phi[ind4D] + fields_chi.B[ind4D]);
                                    ind6D = get_flat_c(is, 1, 0, ix, iy, iz);
                                    h[ind6D] += var_var.q[is] / var_var.T[is] * fields_chi.B[ind4D];

                                }
                                if (global_nm_index[im] == 1) {
                                    ind6D = get_flat_c(is, 0, 1, ix, iy, iz);
                                    h[ind6D] += var_var.q[is] / var_var.T[is] * sqrt(0.5) * fields_chi.A[ind4D];
                                }
                            }
                        }

                    }
                }
            }
            break;
        default:
            printf("[MPI process %d] error in getH function! Aborting...", mpi_my_rank);
            exit(1);
    }
};

/***************************************
 * \fn distrib_getG(COMPLEX *g, const COMPLEX *h)
 * \brief computes g from h
 * \param g: complex array to store g
 * \param h: complex array with h
 *
 * computes modified gyrokinetic distribution function g
 * from gyrokinetic distribution function h.
 * Please note that before calling this function
 * gyrokinetic potentials must be computed
 ***************************************/
void distrib_getG(COMPLEX *g, const COMPLEX *h) {
    size_t ind6D;
    size_t ind4D;
    for (size_t i = 0; i < array_local_size.total_comp; i++) {
        g[i] = h[i];
    }
    switch (systemType) {
        case ELECTROSTATIC:
            break;
        case ELECTROMAGNETIC:
            for (size_t ix = 0; ix < array_local_size.nkx; ix++) {
                for (size_t iy = 0; iy < array_local_size.nky; iy++) {
                    for (size_t iz = 0; iz < array_local_size.nkz; iz++) {
                        for (size_t is = 0; is < array_local_size.ns; is++) {
                            ind4D = ix * array_local_size.nky * array_local_size.nkz * array_local_size.ns +
                                    iy * array_local_size.nkz * array_local_size.ns +
                                    iz * array_local_size.ns +
                                    is;
                            for (size_t im = 0; im < array_local_size.nm; im++) {
                                if (global_nm_index[im] == 0) {
                                    ind6D = get_flat_c(is, 0, 0, ix, iy, iz);
                                    g[ind6D] -= var_var.q[is] / var_var.T[is] *
                                                (fields_chi.phi[ind4D] + fields_chi.B[ind4D]);
                                    ind6D = get_flat_c(is, 1, 0, ix, iy, iz);
                                    g[ind6D] -= var_var.q[is] / var_var.T[is] * fields_chi.B[ind4D];

                                }
                                if (global_nm_index[im] == 1) {
                                    ind6D = get_flat_c(is, 0, 1, ix, iy, iz);
                                    g[ind6D] -= var_var.q[is] / var_var.T[is] * sqrt(0.5) * fields_chi.A[ind4D];
                                }
                            }
                        }

                    }
                }
            }
            break;
        default:
            printf("[MPI process %d] error in getG function! Aborting...", mpi_my_rank);
            exit(1);
    }

};

/***************************************
 * \fn void distrib_getXGrad(const COMPLEX *in, COMPLEX *out):
 * \brief Computes gradient in kx direction
 * \param in: complex array. Distribution function of which gradient will be taken
 * \param out: complex array, where gradient is stored
 *
 *  Computes gradient in kx direction as following:\n
 *  grad(f) = i * kx * f
 *
 ***************************************/
void distrib_getXGrad(const COMPLEX *in, COMPLEX *out) {
    for (size_t ix = 0; ix < array_local_size.nkx; ix++) {
        for (size_t iy = 0; iy < array_local_size.nky; iy++) {
            for (size_t iz = 0; iz < array_local_size.nkz; iz++) {
                for (size_t im = 0; im < array_local_size.nm; im++) {
                    for (size_t il = 0; il < array_local_size.nl; il++) {
                        for (size_t is = 0; is < array_local_size.ns; is++) {
                            out[get_flat_c(is, il, im, ix, iy, iz)] = space_iKx[ix] * in[get_flat_c(is, il, im, ix, iy, iz)];
                        }
                    }
                }
            }
        }
    }
};

/***************************************
 * \fn void distrib_getYGrad(const COMPLEX *in, COMPLEX *out):
 * \brief Computes gradient in ky direction
 * \param in: complex array. Distribution function of which gradient will be taken
 * \param out: complex array, where gradient is stored
 *
 *  Computes gradient in ky direction as following:\n
 *  grad(f) = i * ky * f
 *
 ***************************************/
void distrib_getYGrad(const COMPLEX *in, COMPLEX *out) {
    for (size_t ix = 0; ix < array_local_size.nkx; ix++) {
        for (size_t iy = 0; iy < array_local_size.nky; iy++) {
            for (size_t iz = 0; iz < array_local_size.nkz; iz++) {
                for (size_t im = 0; im < array_local_size.nm; im++) {
                    for (size_t il = 0; il < array_local_size.nl; il++) {
                        for (size_t is = 0; is < array_local_size.ns; is++) {
                            out[get_flat_c(is, il, im, ix, iy, iz)] = space_iKy[iy] * in[get_flat_c(is, il, im, ix, iy, iz)];
                        }
                    }
                }
            }
        }
    }
};

/***************************************
 * \fn void distrib_getZGrad(const COMPLEX *in, COMPLEX *out):
 * \brief Computes gradient in kz direction
 * \param in: complex array. Distribution function of which gradient will be taken
 * \param out: complex array, where gradient is stored
 *
 *  Computes gradient in kz direction as following:\n
 *  grad(f) = i * kz * f
 *
 ***************************************/
void distrib_getZGrad(const COMPLEX *in, COMPLEX *out) {
    for (size_t ix = 0; ix < array_local_size.nkx; ix++) {
        for (size_t iy = 0; iy < array_local_size.nky; iy++) {
            for (size_t iz = 0; iz < array_local_size.nkz; iz++) {
                for (size_t im = 0; im < array_local_size.nm; im++) {
                    for (size_t il = 0; il < array_local_size.nl; il++) {
                        for (size_t is = 0; is < array_local_size.ns; is++) {
                            out[get_flat_c(is, il, im, ix, iy, iz)] = space_iKz[iz] * in[get_flat_c(is, il, im, ix, iy, iz)];
                        }
                    }
                }
            }
        }
    }
};

/***************************************
 * \fn void distrib_enforceReality(COMPLEX *f)
 * \brief enforces reality condition on distribution function array
 * \param f: complex array for which reality condition will be forced.
 *
 * Enforces reality condition f(k) = conj(f(-k)) in plane kz = 0.
 * For a given kx, it first checks where modes -kx are located using the
 * #mpi_whereIsX function:\n
 * <tt> where_neg = mpi_whereIsX[kxNeg * 2];</tt>\n
 * If -kx is stored on a different processor,
 * MPI_VECTOR with a 4D data slice f(kx,kz = 0) is sent to this processor, to the <tt>buffer</tt> array: \n
 * <tt>mpi_sendVector(&f[ind6D],buffer,where_pos,where_neg);</tt>\n
 * if the data stored on the same processor, no vector is being sent.
 * Reality condition is fulfilled in a loop over all other coordinates:\n
 * <tt>f[ind6D_neg] = conj(buffer[ind6D_pos]);</tt>
 ***************************************/
void distrib_enforceReality(COMPLEX *f){
    size_t local_kxPosInd;
    size_t local_kxNegInd;
    size_t kyNeg;
    size_t kxNeg;
    size_t ind6D;
    size_t ind6D_pos;
    size_t ind6D_neg;
    int where_pos;
    int where_neg;
    COMPLEX *buffer = malloc(array_local_size.nky * array_local_size.nm * array_local_size.ns * array_local_size.nl * sizeof(*buffer));
    for (size_t ix = 0; ix < array_global_size.nkx/2 + 1; ix++){
        where_pos = mpi_whereIsX[ix * 2];
        local_kxPosInd = mpi_whereIsX[ix * 2 + 1];
        kxNeg = (ix == 0)  ? ix : array_global_size.nkx - ix;
        where_neg = mpi_whereIsX[kxNeg * 2];
        local_kxNegInd = mpi_whereIsX[kxNeg * 2 + 1];
        //printf("%zu %d\n",ix,where_neg);
        if (where_pos == where_neg){
            for(size_t iy = 0; iy < array_local_size.nky; iy++){
                    for(size_t im = 0; im < array_local_size.nm; im++){
                        for(size_t il = 0; il < array_local_size.nl; il++){
                            for(size_t is = 0; is < array_local_size.ns; is++){
                                kyNeg = (iy == 0)  ? iy : array_local_size.nky - iy;
                                ind6D_neg = get_flat_c(is,il,im,local_kxNegInd,kyNeg,0);
                                ind6D_pos = get_flat_c(is,il,im,local_kxPosInd,iy,0);
                                f[ind6D_neg] = conj(f[ind6D_pos]);
                            }
                        }
                    }
            }
        }
        else{
            ind6D = get_flat_c(0,0,0,local_kxPosInd,0,0);
            //printf("[MPI process %d], from %d to %d ix = %zu, kx = %zu, buf sent = (%f,%f)\n",mpi_my_rank,where_pos, where_neg,ix,local_kxPosInd,creal(f[ind6D]),cimag(f[ind6D]));
             mpi_sendVector(&f[ind6D],buffer,where_pos,where_neg);
             if(mpi_my_row_rank == where_neg){
                 //printf("[MPI process %d], buf recieved = (%f,%f)\n",mpi_my_rank,creal(buffer[0]),cimag(buffer[0]));
                 for(size_t iy = 0; iy < array_local_size.nky; iy++){
                     for(size_t im = 0; im < array_local_size.nm; im++){
                         for(size_t il = 0; il < array_local_size.nl; il++){
                             for(size_t is = 0; is < array_local_size.ns; is++){
                                 kyNeg = (iy == 0)  ? iy : array_local_size.nky - iy;
                                 ind6D_neg = get_flat_c(is,il,im,local_kxNegInd,kyNeg,0);
                                 ind6D_pos = iy * array_local_size.nm * array_local_size.nl * array_local_size.ns +
                                             im * array_local_size.nl * array_local_size.ns +
                                             il * array_local_size.ns +
                                             is;
                                 f[ind6D_neg] = conj(buffer[ind6D_pos]);
                                 //if (iy == 0 && local_kxNegInd == 0 && im == 0 && il==0 && is == 0) printf("[MPI process %d]\n",mpi_my_rank);
                             }
                         }
                     }
                 }
             }
        }
    }
    free(buffer);
};

/***************************************
 * \fn void distrib_setZeroNHalf(COMPLEX *f)
 * \brief sets all Nk/2 modes to zero
 * \param f: complex array
 *
 * sets Nkx/2, Nky/2 and Nz/2 modes of distribution function to zero.
 * Due to reality condition, for kz yhe last mode should be set to zero.
 ***************************************/
void distrib_setZeroNHalf(COMPLEX *f){
    size_t ind6D;
    //setting f(nky/2) = 0;
    for(size_t ix = 0; ix < array_local_size.nkx; ix++){
        for(size_t iz = 0; iz < array_local_size.nkz; iz++){
            for(size_t im = 0; im < array_local_size.nm; im++){
                for(size_t il = 0; il < array_local_size.nl; il++){
                    for(size_t is = 0; is < array_local_size.ns; is++){
                        ind6D = get_flat_c(is,il,im,ix,array_global_size.nky/2,iz);
                        f[ind6D] = 0;
                    }
                }
            }
        }
    }
    //setting f(nz/2) = f(nkz) = 0;
    for(size_t ix = 0; ix < array_local_size.nkx; ix++){
        for(size_t iy = 0; iy < array_local_size.nky; iy++){
            for(size_t im = 0; im < array_local_size.nm; im++){
                for(size_t il = 0; il < array_local_size.nl; il++){
                    for(size_t is = 0; is < array_local_size.ns; is++){
                        ind6D = get_flat_c(is,il,im,ix,iy,array_local_size.nkz - 1);
                        f[ind6D] = 0;
                    }
                }
            }
        }
    }
    //setting f(nkx/2) = 0
    int local_kx = mpi_whereIsX[2 * array_global_size.nkx/2 + 1];
    int iproc = mpi_whereIsX[2 * array_global_size.nkx/2];
    if (mpi_my_row_rank == iproc){
        for(size_t iy = 0; iy < array_local_size.nky; iy++){
            for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                for(size_t im = 0; im < array_local_size.nm; im++){
                    for(size_t il = 0; il < array_local_size.nl; il++){
                        for(size_t is = 0; is < array_local_size.ns; is++){
                            ind6D = get_flat_c(is,il,im,local_kx,iy,iz);
                            f[ind6D] = 0;
                        }
                    }
                }
            }
        }
    }
};