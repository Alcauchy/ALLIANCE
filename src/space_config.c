/**************************************
* @file space_config.c
* \brief space configuration module
*
* creates k and m spaces
***************************************/
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
#define VERBOSE 0

#include "space_config.h"
#include <complex.h>

#define MINUS_I -1.j
double space_Lx;
double space_Ly;
double space_Lz;
double space_LperpMax;
double space_LperpMin;
double space_dx;
double space_dy;
double space_dz;
double space_dKx;
double space_dKy;
double space_dKz;
double space_kXmax;
double space_kYmax;
double space_kZmax;
double space_kPerpMax;
double *space_kx;
double *space_ky;
double *space_kz;
double *space_kPerp;
double *space_kPerp2;
double *space_kSq;
double *space_sqrtM;
double *space_zerosKx;
double *space_zerosKy;
double *space_zerosKz;
size_t *space_globalMIndex;
COMPLEX *space_iKx;
COMPLEX *space_iKy;
COMPLEX *space_iKz;

/***************************************
 * \fn void space_init():
 * \brief initializes wave space. Called in init_start() function.
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
 * \fn void space_generateWaveSpace():
 * \brief generates wave space.
 *
 *  generates wave number arrays space_kx, space_ky, space_kz
 *  of lengths nkx,nky,nkz for a numerical box of size [lx, ly, lz]
 *  in kx,ky,kz directions as following:\n
 * <tt> [0, pi / lx, 2 pi / lx, ... , (n / 2 + 1) pi / lx, - (n / 2) pi / lx, ... ,  - pi / lx] </tt>
 *  generates arrays <tt>space_iKx</tt>, <tt>space_iKy</tt>, <tt>space_iKz</tt>,
 *  of lengths nkx,nky,nkz. These arrays are later used to compute gradients by
 *  #fields_getGradX,
 *  #fields_getGradY,
 *  #distrib_getXGrad,
 *  #distrib_getYGrad,
 *  #distrib_getZGrad.
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

    space_kSq = malloc(array_local_size.nkx * array_local_size.nky * array_local_size.nkz *
                          sizeof(*space_kSq));

    space_zerosKx = malloc(array_local_size.nkx * sizeof(*space_zerosKx));
    space_zerosKy = malloc(array_local_size.nky * sizeof(*space_zerosKy));
    space_zerosKz = malloc(array_local_size.nkz * sizeof(*space_zerosKz));

    space_Lx = parameters.Lx;
    space_Ly = parameters.Ly;
    space_Lz = parameters.Lz;
    space_LperpMax =  (space_Lx > space_Ly) ? space_Lx : space_Ly;

    space_dx = space_Lx / array_global_size.nx;
    space_dy = space_Ly / array_global_size.ny;
    space_dz = space_Lz / array_global_size.nz;
    space_LperpMin =  (space_dx < space_dy) ? space_dx : space_dy;

    space_dKx = 2. * M_PI / ( space_Lx);
    space_dKy = 2. * M_PI / ( space_Ly);
    space_dKz = 2. * M_PI / ( space_Lz);

    space_kXmax = 0.5 * space_dKx * array_global_size.nkx;
    space_kYmax = 0.5 * space_dKy * array_global_size.nky;
    space_kZmax = 0.5 * space_dKz * array_global_size.nkz;

    space_kPerpMax = (space_kXmax > space_kYmax) ? space_kXmax : space_kYmax;

    for (size_t ix = 0; ix < array_local_size.nkx; ix++) {
        if (global_nkx_index[ix] < array_global_size.nkx / 2 + 1) {
            space_kx[ix] = space_dKx * global_nkx_index[ix];
        } else {
            space_kx[ix] = space_dKx * ((double) global_nkx_index[ix] - (double) array_global_size.nkx);
        }
        space_iKx[ix] = 1.j * space_kx[ix];
        //filling space_zerosKx array either with zeros or ones - needed for dealiasing
        if(fabs(space_kx[ix]) > space_dKx * array_global_size.nkx / 3.){
            space_zerosKx[ix] = 0.;
        }
        else{
            space_zerosKx[ix] = 1.;
        }
    }

    for (size_t iy = 0; iy < array_global_size.nky; iy++) {
        if (iy < array_global_size.nky / 2 + 1) {
            space_ky[iy] = space_dKy * iy;
        } else {
            space_ky[iy] = space_dKy * ((double) iy - (double) array_global_size.nky);
        }
        space_iKy[iy] = 1.j * space_ky[iy];
        //filling space_zerosKy array either with zeros or ones - needed for dealiasing
        if(fabs(space_ky[iy]) > space_dKy * array_global_size.nky / 3.){
            space_zerosKy[iy] = 0.;
        }
        else{
            space_zerosKy[iy] = 1.;
        }
    }
    for (size_t iz = 0; iz < array_local_size.nkz; iz++) {
        space_kz[iz] = space_dKz * iz;
        space_iKz[iz] = 1.j * space_kz[iz];
        //filling space_zerosKz array either with zeros or ones - needed for dealiasing
        if(fabs(space_ky[iz]) > space_dKz * array_global_size.nz / 3.){
            space_zerosKz[iz] = 0.;
        }
        else{
            space_zerosKz[iz] = 1.;
        }
    }

    for (size_t ix = 0; ix < array_local_size.nkx; ix++) {
        for (size_t iy = 0; iy < array_local_size.nky; iy++) {
            space_kPerp2[ix * array_local_size.nky + iy] = space_kx[ix] * space_kx[ix] + space_ky[iy] * space_ky[iy];
            space_kPerp[ix * array_local_size.nky + iy] = sqrt(space_kPerp2[ix * array_local_size.nky + iy]);
        }
    }

    for (size_t ix = 0; ix < array_local_size.nkx; ix++) {
        for (size_t iy = 0; iy < array_local_size.nky; iy++) {
            for (size_t iz = 0; iz < array_local_size.nkz; iz++) {
                size_t ind3D = ix * array_local_size.nky * array_local_size.nkz
                               + iy * array_local_size.nkz
                               + iz;
                space_kSq[ind3D] = space_kx[ix] * space_kx[ix] + space_ky[iy] * space_ky[iy] + space_kz[iz] * space_kz[iz];
            }
        }
    }

};

/***************************************
 * \fn space_generateMSpace()
 * \brief generates Hermite space.
 *
 * to be added
 ***************************************/
void space_generateMSpace(){
    space_sqrtM = malloc((array_local_size.nm + 1)  * sizeof(*space_sqrtM));

    for(size_t i = 0; i < array_local_size.nm; i++){
        space_sqrtM[i] = sqrt(space_globalMIndex[i] * 0.5);
    }
    space_sqrtM[array_local_size.nm] = sqrt((space_globalMIndex[array_local_size.nm - 1] + 1.) * 0.5);
}

/***************************************
 * \fn free_wavespace()
 * \brief deallocates all the arrays.
 *
 * to be added...
 ***************************************/
void free_wavespace() {
    free(space_kx);
    free(space_ky);
    free(space_kz);
    free(space_iKx);
    free(space_iKy);
    free(space_iKz);
}
