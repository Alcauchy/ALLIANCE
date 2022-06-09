/**************************************
* @file array.c
*   \brief array manipulation module
*
*   contains functions which are supposed
*   to make array manipulation simpler
***************************************/
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
 * \fn size_t get_flat_c(size_t is, size_t il, size_t im, size_t ix, size_t iy, size_t iz)
 * \brief returns flat index of the element of complex 6D array
 * \param is: species type
 * \param il: Laguerre moment
 * \param im: Hermite moment
 * \param ix: kx index
 * \param iy: ky index
 * \param iz: kz index
 *
 * returns flattened index of a complex array from its 6D index.
 * Flattened index then can be passed to distribution function
 * 6D array to get a required element at position (is,il,im,ix,iy,iz).
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
 * \fn size_t getIndChiBufEM_c(size_t ix,size_t iy, size_t iz, size_t is, size_t ifield)
 * \brief returns flat index of an element of electromagnetic gyrokinetic potential in FOURIER SPACE
 * \param ix: kx index
 * \param iy: ky index
 * \param iz: kz index
 * \param is: particle species index
 * \param ifield: field type
 *
 * returns flattened index of a gyrokinetic potential \f$\chi^{\phi,A,B}\f$ from its 4D index in FOURIER SPACE.
 * flattened index is then can be used to access required value of the gyrokinetic potential at position (ix,iy,iz,is).
 * Type of gyrokinetic potential is specified by ifield parameter. Use 0 is to access \f$\chi^{\phi}(\mathbf{k})\f$,
 * 1 to access \f$\chi^{A}(\mathbf{k})\f$ and 2 to access \f$\chi^{B}(\mathbf{k})\f$.
 ***************************************/
size_t getIndChiBufEM_c(size_t ix,size_t iy, size_t iz, size_t is, size_t ifield) {
    return ix * array_local_size.nky * array_local_size.nkz * array_local_size.ns * CHI_EM +
           iy * array_local_size.nkz * array_local_size.ns * CHI_EM +
           iz * array_local_size.ns * CHI_EM +
           is * CHI_EM + ifield;
}

/***************************************
 * \fn size_t getIndChiBufEM_r(size_t ix,size_t iy, size_t iz, size_t is, size_t ifield)
 * \brief returns flat index of an element of electromagnetic gyrokinetic potential in POSITION SPACE
 * \param ix: x index
 * \param iy: y index
 * \param iz: z index
 * \param is: particle species index
 * \param ifield: field type
 *
 * returns flattened index of a gyrokinetic potential \f$\chi^{\phi,A,B}(\mathbf{k})\f$ from its 4D index in POSITION SPACE.
 * flattened index is then can be used to access required value of the gyrokinetic potential at position (ix,iy,iz,is).
 * Type of gyrokinetic potential is specified by ifield parameter. Use 0 is to access \f$\chi^{\phi}(\mathbf{r})\f$,
 * 1 to access \f$\chi^{A}(\mathbf{r})\f$ and 2 to access \f$\chi^{B}(\mathbf{r})\f$.
 ***************************************/
size_t getIndChiBufEM_r(size_t ix,size_t iy, size_t iz, size_t is, size_t ifield) {
    return ix * array_local_size.nky * (array_local_size.nz + FFT_OFFSET) * array_local_size.ns * CHI_EM +
           iy * (array_local_size.nz + FFT_OFFSET) * array_local_size.ns * CHI_EM +
           iz * array_local_size.ns * CHI_EM +
           is * CHI_EM + ifield;
}


/***************************************
 * \fn size_t getIndChiBufEL_c(size_t ix,size_t iy, size_t iz, size_t is)
 * \brief returns returns flat index of an element of electrostatic gyrokinetic potential in FOURIER SPACE
 * \param ix: kx index
 * \param iy: ky index
 * \param iz: kz index
 * \param is: particle species index
 *
 * returns flattened index of a gyrokinetic potential \f$\chi^{\phi}(\mathbf{k})\f$ from its 4D index in FOURIER SPACE.
 * flattened index is then can be used to access required value of the gyrokinetic potential at position (ix,iy,iz,is).
 ***************************************/
size_t getIndChiBufEL_c(size_t ix,size_t iy, size_t iz, size_t is) {
    return ix * array_local_size.nky * array_local_size.nkz * array_local_size.ns * CHI_EL +
           iy * array_local_size.nkz * array_local_size.ns * CHI_EL +
           iz * array_local_size.ns * CHI_EL +
           is * CHI_EL;
}

/***************************************
 * \fn size_t getIndChiBufEL_r(size_t ix,size_t iy, size_t iz, size_t is)
 * \brief returns returns flat index of an element of electrostatic gyrokinetic potential in REAL SPACE
 * \param ix: x index
 * \param iy: y index
 * \param iz: z index
 * \param is: particle species index
 *
 * returns flattened index of a gyrokinetic potential \f$\chi^{\phi}(\mathbf{r})\f$ from its 4D index in REAL SPACE.
 * flattened index is then can be used to access required value of the gyrokinetic potential at position (ix,iy,iz,is).
 ***************************************/
size_t getIndChiBufEL_r(size_t ix,size_t iy, size_t iz, size_t is) {
    return ix * array_local_size.nky * (array_local_size.nz + FFT_OFFSET) * array_local_size.ns * CHI_EL +
           iy * (array_local_size.nz + FFT_OFFSET) * array_local_size.ns * CHI_EL +
           iz * array_local_size.ns * CHI_EL +
           is * CHI_EL;
}

/***************************************
 * \fn size_t get_flat_r(size_t is, size_t il, size_t im, size_t ix, size_t iy, size_t iz)
 * \brief returns flat index of the element of real 6D array
 * \param is: species type
 * \param il: Laguerre moment
 * \param im: Hermite moment
 * \param ix: x index
 * \param iy: y index
 * \param iz: z index
 *
 * returns flattened index of a real array from its 6D index.
 * Flattened index then can be passed to distribution function
 * 6D array to get a required element at position (is,il,im,ix,iy,iz).
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
 * \fn size_t get_flatIndexComplex3D(size_t ix, size_t iy, size_t iz)
 * \brief returns flat array of complex 3D array
 * \param ix: kx index
 * \param iy: ky index
 * \param iz: kz index
 *
 * returns flattened index of a complex array from its 3D position index.
 * Flattened index then can be passed to one of the fields (\f$ \phi(\mathbf{k}), A_{||}(\mathbf{k}), B_{||}(\mathbf{k}) \f$)
 * 6D array to get a required element at position (ix,iy,iz).
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

