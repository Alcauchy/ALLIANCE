//
// Created by alcauchy on 20/01/2022.
//

#ifndef ALLIANCE_ALPHA_1_0_DIAGNOSTICS_H
#define ALLIANCE_ALPHA_1_0_DIAGNOSTICS_H

#include <stdio.h> // standard C lib
#include <stdlib.h> //Standard C lib
#include <time.h>
#include "array.h"

void diag_computeSpectra(COMPLEX *data, double *spec);                 // compute spectra spec from complex array data.
double diag_computeFreeEnergy(COMPLEX *g,COMPLEX *h);                  // computes free energy from the complex fields g and h
double diag_computeFreeEnergyFields(COMPLEX *h, COMPLEX *fields);      // comp[utes free energy from field data fields and h distribution function
void diag_computeKSpectrum();
void diag_computeMSpectrum();
#endif //ALLIANCE_ALPHA_1_0_DIAGNOSTICS_H
