//
// Created by alcauchy on 20/01/2022.
//

#ifndef ALLIANCE_ALPHA_1_0_DIAGNOSTICS_H
#define ALLIANCE_ALPHA_1_0_DIAGNOSTICS_H

#include <stdio.h> // standard C lib
#include <stdlib.h> //Standard C lib
#include <time.h>
#include "array.h"

void compute_spectra(COMPLEX *data, double *spec);                 // compute spectra spec from complex array data.
double diag_compute_free_energy(COMPLEX *data);                         // computes free energy from the complex data.


#endif //ALLIANCE_ALPHA_1_0_DIAGNOSTICS_H
