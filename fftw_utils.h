//
// Created by alcauchy on 08/12/2021.
//

#ifndef ALLIANCE_ALPHA_1_0_FFTW_UTILS_H
#define ALLIANCE_ALPHA_1_0_FFTW_UTILS_H
#include <stdio.h> // standard C lib
#include <stdlib.h> //Standard C lib
#include <math.h> //Math constants
#include <complex.h> //Complex numbers
#include <fftw3-mpi.h>
#include <string.h>
#include "mpi_utils.h"
#include "array.h"


void fftw_init(MPI_Comm communicator);                      // initializes fftw mpi environment
void fftw_kill();                                           // kills environment
void fftw_dealiasing();                                     // dealiasing, yet to be implemented
void fftw_r2c(double *data_r, COMPLEX *data_c);             // real to complex transform; copies data_r array to fftw buffer r_d, performs transform to complex and copies the resulting array to data_c array
void fftw_c2r(COMPLEX *data_c, double *data_r);             // complex to real transform; copies data_c array to fftw buffer c_d, performs transform to complex and copies the resulting array to data_r array
void fftw_copy_buffer_r(double *ar1, double *ar2);          // copies real array to another real array, function used in fftw_r2c/fftw_c2r to copy input array to buffer and/or back;
void fftw_copy_buffer_c(COMPLEX *ar1, COMPLEX *ar2);        // copies complex array to another complex array, function used in fftw_r2c/fftw_c2r to copy input array to buffer and/or back
void fftw_normalise_data(double *data);                     // normilises real data after fftw by fftw_norm = nkx*nky*nz.
void fftw_test_fill(double *ar,double f);                   // fills input real array with cos(2*pi*f*x).
void dealiasing23(COMPLEX *data_c);                         // performs dealiasing by employing 2/3 rule. Turns all the aliased modes to zero.
double cosinus(double f,int ix);                            // returns cos(2*pi*f*ix)

extern int *global_nkx_index;
#endif //ALLIANCE_ALPHA_1_0_FFTW_UTILS_H
