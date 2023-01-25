#ifndef ALLIANCE_ALPHA_1_0_DIAGNOSTICS_H
#define ALLIANCE_ALPHA_1_0_DIAGNOSTICS_H

#include <stdio.h> // standard C lib
#include <stdlib.h> //Standard C lib
#include <time.h>
#include "array.h"
#include "parameters_io.h"
#include "space_config.h"
#include "equation.h"

void diag_computeSpectra(const COMPLEX *g, const COMPLEX *h, int timestep);                 // compute spectra spec from complex array data.
void diag_computeFreeEnergy(COMPLEX *g, COMPLEX *h);                  // computes free energy from the complex fields g and h
double diag_computeFreeEnergyFields(COMPLEX *h, COMPLEX *fields);      // comp[utes free energy from field data fields and h distribution function
void diag_initSpec();
void diag_computeKSpectrum(const COMPLEX *g, const COMPLEX *h, double *spec);
void diag_computeMSpectrum(const COMPLEX *g, const COMPLEX *h, double *spec);
void diag_computeFieldSpectrum();
void diag_computeHSpectrum(const COMPLEX *h);
void diag_initSpec();
void diag_getShells();
void diag_compute(COMPLEX *g, COMPLEX *h, int timestep);
void diag_computeEnergyBalance(const COMPLEX *h, const COMPLEX *g);
void diag_computeEnergy(const COMPLEX *h);
void diag_print(const COMPLEX *h, const COMPLEX *g, int it);
void diag_filterK(const COMPLEX *in, COMPLEX *out, double k_c);
void diag_filterKHP(const COMPLEX *in, COMPLEX *out, double k_c);
void diag_filterKz(const COMPLEX *in, COMPLEX *out, double k_c);
void diag_filterKzHP(const COMPLEX *in, COMPLEX *out, double k_c);
void diag_computeNonlinearFlux(const COMPLEX *in);

extern double *diag_mSpec;
extern double *diag_kSpec;
extern double *diag_kSpecPhi;
extern double *diag_kSpecBperp;
extern double *diag_kSpecBpar;
extern double *diag_kSpecH;
extern double *diag_shells;
extern double *diag_shellNorm;
extern double *diag_shellCentres;
extern double *diag_MM;
extern double *diag_nonlinearFlux;
extern double *diag_nonlinearFluxInverse;
extern double *diag_nonlinearFluxForward;
extern double *diag_nonlinearNorm;

extern double diag_freeEnergy;
extern double diag_free_energy0;
extern int diag_numOfShells;
extern int diag_numOfShellBounds;
extern double *diag_energyH;
extern double diag_energyPhi;
extern double diag_energyBperp;
extern double diag_energyBpar;
extern double diag_energyTotal;
extern double diag_injected;
extern double diag_dissipated;
extern double diag_dissipated_kPerp;
extern double diag_dissipated_kZ;
extern double diag_dissipated_m;
extern double diag_etakmax;
extern double diag_LInt;
extern double diag_LDis;

#endif //ALLIANCE_ALPHA_1_0_DIAGNOSTICS_H
