////////////////////////////////////////////////////////////////////////////////
// 26/01/2022 created by Gene Gorbunov
//                                   ARRAY
//
// alloc_double6D
// alloc_complex6D
// get_flat_c
// get_flat_r
// getIndChiBufEM_c
// getIndChiBufEM_r
// getIndChiBufEL_c
// getIndChiBufEL_r
// get_flatIndexComplex3D
// multiply_ar_c
// multiply_ar_r
// fill_rand
// fill_randM0
// fill_randSingleKM
// sinus
//
// VERSION 1.0
////////////////////////////////////////////////////////////////////////////////
#include "array.h"
#include "utils_fftw.h"
#include "space_config.h"

#define CHI_EM 3
#define CHI_EL 1
#define FFT_OFFSET 2
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
 * getIndChiBufEM_c(size_t ix,size_t iy, size_t iz, size_t is, size_t ifield)
 ***************************************/
size_t getIndChiBufEM_c(size_t ix,size_t iy, size_t iz, size_t is, size_t ifield) {
    return ix * array_local_size.nky * array_local_size.nkz * array_local_size.ns * CHI_EM +
           iy * array_local_size.nkz * array_local_size.ns * CHI_EM +
           iz * array_local_size.ns * CHI_EM +
           is * CHI_EM + ifield;
}

/***************************************
 * getIndChiBufEM_r(size_t ix,size_t iy, size_t iz, size_t is, size_t ifield)
 ***************************************/
size_t getIndChiBufEM_r(size_t ix,size_t iy, size_t iz, size_t is, size_t ifield) {
    return ix * array_local_size.nky * (array_local_size.nz + FFT_OFFSET) * array_local_size.ns * CHI_EM +
           iy * (array_local_size.nz + FFT_OFFSET) * array_local_size.ns * CHI_EM +
           iz * array_local_size.ns * CHI_EM +
           is * CHI_EM + ifield;
}


/***************************************
 * getIndChiBufEL_c(size_t ix,size_t iy, size_t iz, size_t is)
 ***************************************/
size_t getIndChiBufEL_c(size_t ix,size_t iy, size_t iz, size_t is) {
    return ix * array_local_size.nky * array_local_size.nkz * array_local_size.ns * CHI_EL +
           iy * array_local_size.nkz * array_local_size.ns * CHI_EL +
           iz * array_local_size.ns * CHI_EL +
           is * CHI_EL;
}

/***************************************
 * getIndChiBufEL_r(size_t ix,size_t iy, size_t iz, size_t is)
 ***************************************/
size_t getIndChiBufEL_r(size_t ix,size_t iy, size_t iz, size_t is) {
    return ix * array_local_size.nky * (array_local_size.nz + FFT_OFFSET) * array_local_size.ns * CHI_EL +
           iy * (array_local_size.nz + FFT_OFFSET) * array_local_size.ns * CHI_EL +
           iz * array_local_size.ns * CHI_EL +
           is * CHI_EL;
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
size_t get_flatIndexComplex3D(size_t ix, size_t iy, size_t iz) {
    return ix * array_offset3D.kx +
           iy * array_offset3D.ky +
           iz;
}

/***************************************
 * getIndChi(size_t ix,size_t iy, size_t iz, size_t is)
 ***************************************/
size_t getIndChi(size_t ix,size_t iy, size_t iz, size_t is) {
    return ix * array_local_size.nky * array_local_size.nkz * array_local_size.ns +
           iy * array_local_size.nkz * array_local_size.ns +
           iz * array_local_size.ns +
           is;
}

/***************************************
 * multiply_ar_c(COMPLEX *ar1, COMPLEX *ar2, COMPLEX *ret)
 ***************************************/
void multiply_ar_c(COMPLEX *ar1, COMPLEX *ar2, COMPLEX *ret) {
    for (size_t i = 0; i < array_local_size.total_comp; i++) {
        ret[i] = ar1[i] * ar2[i];
    }
}

/***************************************
 * multiply_ar_r(const double *ar1, const double *ar2, double *ret)
 ***************************************/
void multiply_ar_r(const double *ar1, const double *ar2, double *ret) {
    for (size_t i = 0; i < array_local_size.total_real; i++) {
        ret[i] = ar1[i] * ar2[i];
    }
}

