//
// Created by alcauchy on 26/01/2022.
//
#include "array.h"

struct array_size array_local_size;
struct array_size array_global_size;

double* alloc_double6D(size_t nkx, size_t nky, size_t nz, size_t nm, size_t nl, size_t ns){
    return malloc(nkx * nky * nz * nl * nm * ns * sizeof(double));
}

COMPLEX* alloc_complex6D(size_t nkx, size_t nky, size_t nkz, size_t nm, size_t nl, size_t ns){
    return malloc(nkx * nky * nkz * nl * nm * ns * sizeof(COMPLEX));
}

size_t get_flat_c(size_t is, size_t il, size_t im, size_t ix, size_t iy, size_t iz){
    return ix * array_local_size.nky * array_local_size.nkz * array_local_size.nm * array_local_size.nl * array_local_size.ns
           + iy * array_local_size.nkz * array_local_size.nm * array_local_size.nl * array_local_size.ns
           + iz * array_local_size.nm * array_local_size.nl * array_local_size.ns
           + im * array_local_size.nl * array_local_size.ns
           + il * array_local_size.ns
           + is;
}

size_t get_flat_r(size_t is, size_t il, size_t im, size_t ix, size_t iy, size_t iz){
    return ix * array_local_size.nky * (array_local_size.nz+2) * array_local_size.nm * array_local_size.nl * array_local_size.ns
           + iy * (array_local_size.nz+2) * array_local_size.nm * array_local_size.nl * array_local_size.ns
           + iz * array_local_size.nm * array_local_size.nl * array_local_size.ns
           + im * array_local_size.nl * array_local_size.ns
           + il * array_local_size.ns
           + is;
}

void multiply_ar_c(COMPLEX *ar1, COMPLEX *ar2, COMPLEX *ret){
    for(size_t ikx = 0; ikx < array_local_size.nkx; ikx++) {
        for (size_t iky = 0; iky < array_local_size.nky; iky++) {
            for (size_t ikz = 0; ikz < array_local_size.nkz; ikz++) {
                for (size_t im = 0; im < array_local_size.nm; im++) {
                    for (size_t il = 0; il < array_local_size.nl; il++) {
                        for (size_t is = 0; is < array_local_size.ns; is++) {
                            ret[get_flat_c(is,il,im,ikx,iky,ikz)] = ar1[get_flat_c(is,il,im,ikx,iky,ikz)]*ar2[get_flat_c(is,il,im,ikx,iky,ikz)];
                        }
                    }
                }
            }
        }
    }
}

void multiply_ar_r(double *ar1, double *ar2, double *ret){
    for(size_t ikx = 0; ikx < array_local_size.nkx; ikx++) {
        for (size_t iky = 0; iky < array_local_size.nky; iky++) {
            for (size_t ikz = 0; ikz < array_local_size.nz+2; ikz++) {
                for (size_t im = 0; im < array_local_size.nm; im++) {
                    for (size_t il = 0; il < array_local_size.nl; il++) {
                        for (size_t is = 0; is < array_local_size.ns; is++) {
                            ret[get_flat_r(is,il,im,ikx,iky,ikz)] = ar1[get_flat_r(is,il,im,ikx,iky,ikz)]*ar2[get_flat_r(is,il,im,ikx,iky,ikz)];
                        }
                    }
                }
            }
        }
    }
}

void fill_rand(COMPLEX *ar1){
    srand(time(NULL));
    for(size_t ikx = 0; ikx < array_local_size.nkx; ikx++) {
        for (size_t iky = 0; iky < array_local_size.nky; iky++) {
            for (size_t ikz = 0; ikz < array_local_size.nkz; ikz++) {
                for (size_t im = 0; im < array_local_size.nm; im++) {
                    for (size_t il = 0; il < array_local_size.nl; il++) {
                        for (size_t is = 0; is < array_local_size.ns; is++) {
                            ar1[get_flat_c(is,il,im,ikx,iky,ikz)] = (0.5 - (double)rand()/(double)(RAND_MAX))+(0.5 - (double)rand()/(double)(RAND_MAX))*1.j;
                        }
                    }
                }
            }
        }
    }
}