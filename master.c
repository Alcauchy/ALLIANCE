//
// Created by alcauchy on 07/12/2021.
//
#include "master.h"
struct array_size array_local_size;
struct system_param parameters;

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



