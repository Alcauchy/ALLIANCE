//
// Created by alcauchy on 26/01/2022.
//

#ifndef ALLIANCE_ALPHA_1_0_ARRAY_H
#define ALLIANCE_ALPHA_1_0_ARRAY_H
#include <stdio.h> // standard C lib
#include <stdlib.h> //Standard C lib
#include <complex.h>
#include <fftw3-mpi.h>
#include <time.h>

typedef fftw_complex COMPLEX;


double* alloc_double6D(size_t nkx, size_t nky, size_t nz, size_t nm, size_t nl, size_t ns);     // allocate 6d real array
COMPLEX* alloc_complex6D(size_t nkx, size_t nky, size_t nkz, size_t nm, size_t nl, size_t ns);  // allocate 6d complex array
size_t get_flat_c(size_t is, size_t il, size_t im, size_t ix, size_t iy, size_t iz);            // get correct array element from real
size_t get_flat_r(size_t is, size_t il, size_t im, size_t ix, size_t iy, size_t iz);
size_t get_flatIndexComplex3D(size_t ix, size_t iy, size_t iz);                                 //
void multiply_ar_c(COMPLEX *ar1, COMPLEX *ar2, COMPLEX *ret);                                   // multiply two complex 6d arrays, ar1 and ar2, and output the result to ret array.
void multiply_ar_r(double *ar1, double *ar2, double *ret);                                      // multiply two return 6d arrays, ar1 and ar2, and output the result to ret array.
void fill_rand(COMPLEX *ar1);                                                                   // fill 6D array ar1 with random values
double sinus(int kx, int ky, int kz, double f);


struct array_size {
    size_t nkx;
    size_t nky;
    size_t nkz;
    size_t nm;
    size_t nl;
    size_t ns;
    size_t nz;
    size_t total_real;
    size_t total_comp;
};

struct offset_size {
    size_t kx;
    size_t ky;
    size_t kz;
    size_t x;
    size_t y;
    size_t z;
    size_t m;
    size_t l;
    size_t s;
};

extern struct array_size array_local_size;
extern struct array_size array_global_size;
extern struct offset_size array_offset;
extern struct offset_size array_offset3D;
#endif //ALLIANCE_ALPHA_1_0_ARRAY_H
