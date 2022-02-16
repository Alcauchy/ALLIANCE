//
// Created by alcauchy on 20/01/2022.
//

#ifndef ALLIANCE_ALPHA_1_0_DIAGNOSTICS_H
#define ALLIANCE_ALPHA_1_0_DIAGNOSTICS_H

#include <stdio.h> // standard C lib
#include <stdlib.h> //Standard C lib
#include <time.h>
#include "array.h"
#include "parameters_io.h"
#include "space_config.h"

void diag_computeSpectra(const COMPLEX *g, const COMPLEX *h, int timestep);                 // compute spectra spec from complex array data.
void diag_computeFreeEnergy(COMPLEX *g, COMPLEX *h, int it);                  // computes free energy from the complex fields g and h
double diag_computeFreeEnergyFields(COMPLEX *h, COMPLEX *fields);      // comp[utes free energy from field data fields and h distribution function
void diag_initSpec();
void diag_computeKSpectrum(const COMPLEX *g, const COMPLEX *h, double *spec);
void diag_computeMSpectrum(const COMPLEX *g, const COMPLEX *h, double *spec);
void diag_initSpec();
void diag_getShells();

extern double *diag_mSpec;
extern double *diag_kSpec;
extern double *diag_shells;
extern double diag_freeEnergy;
#endif //ALLIANCE_ALPHA_1_0_DIAGNOSTICS_H
