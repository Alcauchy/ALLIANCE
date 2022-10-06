#ifndef ALLIANCE_ALPHA_1_0_UTILS_FFTW_H
#define ALLIANCE_ALPHA_1_0_UTILS_FFTW_H
#include <stdio.h> // standard C lib
#include <stdlib.h> //Standard C lib
#include <math.h> //Math constants
#include <complex.h> //Complex numbers
#include <fftw3-mpi.h>
#include <string.h>
#include "utils_mpi.h"
#include "array.h"
#include "init.h"


void fftw_init(MPI_Comm communicator);                      // initializes fftw mpi environment
void fftw_kill();                                           // kills environment
void fftw_r2c();             // real to complex transform; copies data_r array to fftw buffer fftw_hBuf, performs transform to complex and copies the resulting array to data_c array
void fftw_c2r();             // complex to real transform; copies data_c array to fftw buffer fftw_hBuf, performs transform to complex and copies the resulting array to data_r array
void fftw_r2c_chi();
void fftw_c2r_chi();
void fftw_c2r_field();
void fftw_r2c_field();
void fftw_copy_buffer_r(double *to, double *from);          // copies real array to another real array, function used in fftw_r2c/fftw_c2r to copy input array to buffer and/or back;
void fftw_copy_buffer_c(COMPLEX *to, COMPLEX *from);        // copies complex array to another complex array, function used in fftw_r2c/fftw_c2r to copy input array to buffer and/or back
void fftw_copyChiBuf_r(double *ar1, double *ar2);
void fftw_copyChiBuf_c(COMPLEX *ar1, COMPLEX *ar2);
void fftw_copyFieldBuf_r(double *to, double *from);
void fftw_copyFieldBuf_c(COMPLEX *to, COMPLEX *from);
void fftw_normalise_data(COMPLEX *data);                     // normilises real data after fftw by fftw_norm = nkx*nky*nz.
void fftw_normalise_chi_r(double *data);
void fftw_normalise_data_r(double *data);
void fftw_normalise_field_r(double *data);
void fftw_test_fill(double *ar,double f);                   // fills input real array with cos(2*pi*f*x).
void fftw_transposeToXY();
void fftw_transposeToYX();
void fftw_transposeToXY_chi();
void fftw_transposeToYX_chi();
void fftw_transposeToXY_field();
void fftw_transposeToYX_field();
double cosinus(double f,int ix);                            // returns cos(2*pi*f*ix)

extern COMPLEX *fftw_hBuf;
extern COMPLEX *fftw_chiBuf;
extern COMPLEX *fftw_field;
extern double fftw_norm;
extern int *global_nkx_index;

#endif //ALLIANCE_ALPHA_1_0_UTILS_FFTW_H
