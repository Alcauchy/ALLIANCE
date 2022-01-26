//
// Created by alcauchy on 26/01/2022.
//

#ifndef ALLIANCE_ALPHA_1_0_ARRAY_H
#define ALLIANCE_ALPHA_1_0_ARRAY_H
#include <stdio.h> // standard C lib
#include <stdlib.h> //Standard C lib
#include <complex.h>
#include <fftw3-mpi.h>

typedef fftw_complex COMPLEX;


double* alloc_double6D(size_t nkx, size_t nky, size_t nz, size_t nm, size_t nl, size_t ns);     // allocate 6d real array
COMPLEX* alloc_complex6D(size_t nkx, size_t nky, size_t nkz, size_t nm, size_t nl, size_t ns);  // allocate 6d complex array
size_t get_flat_c(size_t is, size_t il, size_t im, size_t ix, size_t iy, size_t iz);            // get correct array element from real
size_t get_flat_r(size_t is, size_t il, size_t im, size_t ix, size_t iy, size_t iz);
double sinus(int kx, int ky, int kz, double f);

struct system_param {
    size_t nkx;
    size_t nky;
    size_t nkz;
    size_t nm;
    size_t nl;
    size_t ns;
    size_t nz;

};
struct array_size {
    size_t nkx;
    size_t nky;
    size_t nkz;
    size_t nm;
    size_t nl;
    size_t ns;
    size_t nz;
};

extern struct array_size array_local_size;
extern struct system_param parameters;
#endif //ALLIANCE_ALPHA_1_0_ARRAY_H
