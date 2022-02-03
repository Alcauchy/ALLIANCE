//
// Created by alcauchy on 03/02/2022.
//
#include "space_config.h"
#define MINUS_I -1.j
double space_Lx = 1.0;
double space_Ly = 1.0;
double space_Lz = 1.0;
double *space_kx;
double *space_ky;
double *space_kz;
COMPLEX *space_iKx;
COMPLEX *space_iKy;
COMPLEX *space_iKz;

/*
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
 */
void space_generateWaveSpace() {
    space_kx = malloc(array_local_size.nkx *
                      sizeof(space_kx));
    space_ky = malloc(array_local_size.nky *
                      sizeof(space_ky));
    space_kz = malloc(array_local_size.nkz *
                      sizeof(space_kz));

    space_iKx = malloc(array_local_size.nkx *
                      sizeof(space_iKx));
    space_iKy = malloc(array_local_size.nky *
                      sizeof(space_iKy));
    space_iKz = malloc(array_local_size.nkz *
                      sizeof(space_iKz));

    double deltaKx = M_PI / (array_global_size.nkx * space_Lx);
    double deltaKy = M_PI / (array_global_size.nky * space_Ly);
    double deltaKz = 2 * M_PI / (array_global_size.nkz * space_Lz);

    for (size_t ix = 0; ix < array_local_size.nkx; ix++)
    {
        if (global_nkx_index[ix] < array_global_size.nkx/2 + 1)
        {
            space_kx[ix] = deltaKx * global_nkx_index[ix];
        }
        else
        {
            space_kx[ix] = deltaKx * ((double)global_nkx_index[ix] - (double)array_global_size.nkx);
        }
        space_iKx[ix] = MINUS_I * space_kx[ix];
        printf("[MPI process %d] kx[%d] = %f\n",mpi_my_rank, ix, space_kx[ix]);
    }

    for (size_t iy = 0; iy < array_global_size.nky; iy++)
    {
        if (iy < array_global_size.nky/2 + 1)
        {
            space_ky[iy] = deltaKy * iy;
        }
        else
        {
            space_ky[iy] = deltaKy * (iy - array_global_size.nky);
        }
        space_iKy[iy] = MINUS_I * space_ky[iy];
        printf("[MPI process %d] ky[%d] = %f\n",mpi_my_rank, iy, space_ky[iy]);
    }

    for (size_t iz = 0; iz < array_local_size.nkz; iz++)
    {
        space_kz[iz] = deltaKz * iz;
        space_iKx[iz] = MINUS_I * space_kx[iz];
        printf("[MPI process %d] kz[%d] = %f\n ",mpi_my_rank, iz, space_kz[iz]);
    }

};

/*
 * space_xGrad(COMPLEX *data):
 *  Computes gradient in kx direction as following:
 *  grad(f) = -2 * pi * i * kx * f
 *  args:
 *  COMPLEX *data : 6D data array gradient of which is going to be taken
 *
 *  return:
 *
 */
void space_xGrad(COMPLEX *data) {
    for(size_t ix = 0; ix < array_local_size.nkx; ix++ )
    {
        for(size_t iy = 0; iy < array_local_size.nky; iy++)
        {
            for(size_t iz = 0; iz < array_local_size.nkz; iz++)
            {
                for(size_t im = 0; im < array_local_size.nm; im++)
                {
                    for(size_t il = 0; il < array_local_size.nl; il++)
                    {
                        for(size_t is = 0; is < array_local_size.ns; is++)
                        {
                            data[get_flat_c(is,il,im,ix,iy,iz)] *= space_iKx[ix];
                        }
                    }
                }
            }
        }
    }
};

/*
 * space_yGrad(COMPLEX *data):
 *  Computes gradient in ky direction as following:
 *  grad(f) = -2 * pi * i * ky * f
 *  args:
 *  COMPLEX *data : 6D data array gradient of which is going to be taken
 *
 *  return:
 *
 */
void space_yGrad(COMPLEX *data) {
    for(size_t ix = 0; ix < array_local_size.nkx; ix++ )
    {
        for(size_t iy = 0; iy < array_local_size.nky; iy++)
        {
            for(size_t iz = 0; iz < array_local_size.nkz; iz++)
            {
                for(size_t im = 0; im < array_local_size.nm; im++)
                {
                    for(size_t il = 0; il < array_local_size.nl; il++)
                    {
                        for(size_t is = 0; is < array_local_size.ns; is++)
                        {
                            data[get_flat_c(is,il,im,ix,iy,iz)] *= space_iKy[iy];
                        }
                    }
                }
            }
        }
    }
};

/*
 * space_zGrad(COMPLEX *data):
 *  Computes gradient in kz direction as following:
 *  grad(f) = -2 * pi * i * kz * f
 *  args:
 *  COMPLEX *data : 6D data array gradient of which is going to be taken
 *
 *  return:
 *
 */
void space_zGrad(COMPLEX *data) {
    for(size_t ix = 0; ix < array_local_size.nkx; ix++ )
    {
        for(size_t iy = 0; iy < array_local_size.nky; iy++)
        {
            for(size_t iz = 0; iz < array_local_size.nkz; iz++)
            {
                for(size_t im = 0; im < array_local_size.nm; im++)
                {
                    for(size_t il = 0; il < array_local_size.nl; il++)
                    {
                        for(size_t is = 0; is < array_local_size.ns; is++)
                        {
                            data[get_flat_c(is,il,im,ix,iy,iz)] *= space_iKz[iz];
                        }
                    }
                }
            }
        }
    }
};