//
// Created by alcauchy on 20/01/2022.
//

#ifndef ALLIANCE_ALPHA_1_0_DIAGNOSTICS_H
#define ALLIANCE_ALPHA_1_0_DIAGNOSTICS_H

#include <stdio.h> // standard C lib
#include <stdlib.h> //Standard C lib
#include <time.h>
#include "array.h"

void diag_multiply_ar(COMPLEX *ar1, COMPLEX *ar2, COMPLEX *ret);        // multiply two complex 6d arrays, ar1 and ar2, and output the result to ret array.
void diag_multiply_ar_r(double *ar1, double *ar2, double *ret);         // multiply two return 6d arrays, ar1 and ar2, and output the result to ret array.
void diag_fill_rand(COMPLEX *ar1);                                      // fill 6D array ar1 with random values
void diag_compute_spectra(COMPLEX *data, double *spec);                 // compute spectra spec from complex array data.
double diag_compute_free_energy(COMPLEX *data);                         // computes free energy from the complex data.

#endif //ALLIANCE_ALPHA_1_0_DIAGNOSTICS_H
