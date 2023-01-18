//
// Created by alcauchy on 03/02/2022.
//

#ifndef ALLIANCE_ALPHA_1_0_INIT_H
#define ALLIANCE_ALPHA_1_0_INIT_H
#include <stdlib.h>
#include <stdio.h>
#include "parameters_io.h"
#include "array.h"
#include "utils_mpi.h"
#include "utils_fftw.h"
#include "utils_hdf.h"
#include "fields.h"
#include "equation.h"
#include "space_config.h"
#include "solver.h"
#include "variables.h"
#include "diagnostics.h"

enum electromagnetic{ELECTROSTATIC,ELECTROMAGNETIC};
enum adiabatic{NONADIABATIC,ADIABATIC};
enum initial{RANDOM, FROMFILE};
enum spectrum{UNIT, SHELL};
enum dealiasing {ALIASED, TWOTHIRDS};

void init_start(char *filename);
void init_printParameters();
void init_initFlags();
void init_conditions(COMPLEX *data);
void init_computation();
void init_physicalSystem();
double init_energySpec(double k, double m, double amp, double disp, double k_z, double disp_z);
double init_sinc(double amp, double f, double x, double y, double z, double x0, double y0, double z0);
double init_exp2(double amp, double f, double x, double y, double z, double x0, double y0, double z0);
double init_sinX(double amp, double f, double x);
double init_cosX(double amp, double f, double x);
void init_fillSinc(COMPLEX *out);
extern enum adiabatic kinetic;
extern enum electromagnetic systemType;
extern enum initial initialConditions;
extern enum spectrum spectrumType;
extern enum dealiasing dealiasingType;
#endif //ALLIANCE_ALPHA_1_0_INIT_H
