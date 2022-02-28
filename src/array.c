//
// Created by alcauchy on 26/01/2022.
//
#include "array.h"

struct array_size array_local_size;
struct array_size array_global_size;
struct offset_size array_offset;
struct offset_size array_offset3D;

/***************************************
 * alloc_double6D(size_t nkx, size_t nky, size_t nz, size_t nm, size_t nl, size_t ns)
 ***************************************/
double *alloc_double6D(size_t nkx, size_t nky, size_t nz, size_t nm, size_t nl, size_t ns) {
    return malloc(nkx * nky * nz * nl * nm * ns * sizeof(double));
}

/***************************************
 * alloc_complex6D(size_t nkx, size_t nky, size_t nkz, size_t nm, size_t nl, size_t ns)
 ***************************************/
COMPLEX *alloc_complex6D(size_t nkx, size_t nky, size_t nkz, size_t nm, size_t nl, size_t ns) {
    return malloc(nkx * nky * nkz * nl * nm * ns * sizeof(COMPLEX));
}

/***************************************
 * get_flat_c(size_t is, size_t il, size_t im, size_t ix, size_t iy, size_t iz)
 ***************************************/
size_t get_flat_c(size_t is, size_t il, size_t im, size_t ix, size_t iy, size_t iz) {
    return ix * array_offset.kx +
           iy * array_offset.ky +
           iz * array_offset.kz +
           im * array_offset.m +
           il * array_offset.l +
           is;
}

/***************************************
 * get_flat_r(size_t is, size_t il, size_t im, size_t ix, size_t iy, size_t iz)
 ***************************************/
size_t get_flat_r(size_t is, size_t il, size_t im, size_t ix, size_t iy, size_t iz) {
    return ix * array_offset.x +
           iy * array_offset.y +
           iz * array_offset.z +
           im * array_offset.m +
           il * array_offset.l +
           is;
}

/***************************************
 * get_flatIndexComplex3D(size_t ix, size_t iy, size_t iz)
 ***************************************/
size_t get_flatIndexComplex3D(size_t ix, size_t iy, size_t iz){
    return ix * array_offset3D.kx +
           iy * array_offset3D.ky +
           iz;
}

/***************************************
 * multiply_ar_c(COMPLEX *ar1, COMPLEX *ar2, COMPLEX *ret)
 ***************************************/
void multiply_ar_c(COMPLEX *ar1, COMPLEX *ar2, COMPLEX *ret) {
    for (size_t i = 0; i < array_local_size.total_comp; i++)
    {
        ret[i] = ar1[i] * ar2[i];
    }
}

/***************************************
 * multiply_ar_r(const double *ar1, const double *ar2, double *ret)
 ***************************************/
void multiply_ar_r(const double *ar1, const double *ar2, double *ret){
    for (size_t i = 0; i < array_local_size.total_real; i++)
    {
        ret[i] = ar1[i] * ar2[i];
    }
}

/***************************************
 * fill_rand(COMPLEX *ar1)
 ***************************************/
void fill_rand(COMPLEX *ar1) {
    srand(time(NULL));
    for (size_t i = 0; i < array_local_size.total_comp; i++)
    {
        ar1[i] = (0.5 - (double) rand() / (double) (RAND_MAX)) + (0.5 - (double) rand() / (double) (RAND_MAX)) * 1.j;
    }
}