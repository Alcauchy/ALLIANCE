////////////////////////////////////////////////////////////////////////////////
// 03/02/2022 created by Gene Gorbunov
//                                   SPACE CONFIGURATION
//
// space_init
// space_generateWaveSpace
// space_generateMSpace
// free_wavespace
//
// VERSION 1.0
////////////////////////////////////////////////////////////////////////////////

#include "space_config.h"
#include <complex.h>

#define MINUS_I -1.j
double space_Lx = 1.0;
double space_Ly = 1.0;
double space_Lz = 1.0;
double *space_kx;
double *space_ky;
double *space_kz;
double *space_kPerp;
double *space_kPerp2;
double *space_sqrtM;
size_t *space_globalMIndex;
COMPLEX *space_iKx;
COMPLEX *space_iKy;
COMPLEX *space_iKz;

/***************************************
 * space_init():
 * initializes wave space. Called in init_start() function.
 *
 *  args:
 *
 *  return:
 *
 ***************************************/
void space_init() {
    space_globalMIndex = malloc(array_local_size.nm * sizeof(*space_globalMIndex));
    for (size_t i = 0; i < array_local_size.nm; i++){
        space_globalMIndex[i] = array_global_size.nm / mpi_dims[0] * mpi_my_col_rank + i;
    }
    space_generateWaveSpace();
    space_generateMSpace();
}

/***************************************
 * space_generateWaveSpace():
 *  generates wave number arrays space_kx, space_ky, space_kz
 *  of lengths nkx,nky,nkz for a numerical box of size [lx, ly, lz]
 *  in kx,ky,kz directions as following:
 *  [0, pi / (nkx * lx), 2 pi / (nkx * lx), ... , (n / 2 + 1) pi / (nkx * lx), - (n / 2) pi / (nkx * lx), ... ,  - pi / (nkx * lx)]
 *
 *  generates arrays space_iKx, space_iKy, space_iKz,
 *  of lengths nkx,nky,nkz. These arrays are later used to compute gradients:
 *  space_iKx = - i * space_kx, etc.
 *
 *  args:
 *
 *  return:
 *
 ***************************************/
void space_generateWaveSpace() {
    space_kx = malloc(array_local_size.nkx *
                      sizeof(*space_kx));
    space_ky = malloc(array_local_size.nky *
                      sizeof(*space_ky));
    space_kz = malloc(array_local_size.nkz *
                      sizeof(*space_kz));

    space_iKx = malloc(array_local_size.nkx *
                       sizeof(*space_iKx));
    space_iKy = malloc(array_local_size.nky *
                       sizeof(*space_iKy));
    space_iKz = malloc(array_local_size.nkz *
                       sizeof(*space_iKz));

    space_kPerp = malloc(array_local_size.nkx * array_local_size.nky *
                         sizeof(*space_kPerp));
    space_kPerp2 = malloc(array_local_size.nkx * array_local_size.nky *
                          sizeof(*space_kPerp2));

    double deltaKx = M_PI / (array_global_size.nkx * space_Lx);
    double deltaKy = M_PI / (array_global_size.nky * space_Ly);
    double deltaKz = 2. * M_PI / (array_global_size.nkz * space_Lz);

    for (size_t ix = 0; ix < array_local_size.nkx; ix++) {
        if (global_nkx_index[ix] < array_global_size.nkx / 2 + 1) {
            space_kx[ix] = deltaKx * global_nkx_index[ix];
        } else {
            space_kx[ix] = deltaKx * ((double) global_nkx_index[ix] - (double) array_global_size.nkx);
        }
        space_iKx[ix] = -1.j * space_kx[ix];
        //printf("[MPI process %d] kx[%d] = %f\n",mpi_my_rank, ix, space_kx[ix]);
    }
    printf("[MPI process %d] kx generated \n", mpi_my_rank);

    for (size_t iy = 0; iy < array_global_size.nky; iy++) {
        if (iy < array_global_size.nky / 2 + 1) {
            space_ky[iy] = deltaKy * iy;
        } else {
            space_ky[iy] = deltaKy * ((double) iy - (double) array_global_size.nky);
        }
        space_iKy[iy] = -1.j * space_ky[iy];
        //printf("[MPI process %d] ky[%d] = %f\n",mpi_my_rank, iy, space_ky[iy]);
    }
    //space_iKy[1] = (0,0);
    printf("[MPI process %d] ky generated \n", mpi_my_rank);
    for (size_t iz = 0; iz < array_local_size.nkz; iz++) {
        space_kz[iz] = deltaKz * iz;
        space_iKz[iz] = -1.j * space_kz[iz];
        //printf("[MPI process %d] kz[%d] = %f\n ",mpi_my_rank, iz, space_kz[iz]);
    }
    printf("[MPI process %d] kz generated \n", mpi_my_rank);

    for (size_t ix = 0; ix < array_local_size.nkx; ix++) {
        for (size_t iy = 0; iy < array_local_size.nky; iy++) {
            space_kPerp2[ix * array_local_size.nky + iy] = space_kx[ix] * space_kx[ix] + space_ky[iy] * space_ky[iy];
            space_kPerp[ix * array_local_size.nky + iy] = sqrt(space_kPerp2[ix * array_local_size.nky + iy]);
        }
    }

};

/***************************************
 * space_generateMSpace()
 ***************************************/
void space_generateMSpace(){
    space_sqrtM = malloc((array_local_size.nm + 1)  * sizeof(*space_sqrtM));

    for(size_t i = 0; i < array_local_size.nm; i++){
        space_sqrtM[i] = sqrt(space_globalMIndex[i] * 0.5);
    }
    space_sqrtM[array_local_size.nm] = sqrt((space_globalMIndex[array_local_size.nm - 1] + 1.) * 0.5);
}

/***************************************
 * free_wavespace()
 ***************************************/
void free_wavespace() {
    free(space_kx);
    free(space_ky);
    free(space_kz);
    free(space_iKx);
    free(space_iKy);
    free(space_iKz);
}
