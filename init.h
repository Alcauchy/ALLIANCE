//
// Created by alcauchy on 03/02/2022.
//

#ifndef ALLIANCE_ALPHA_1_0_INIT_H
#define ALLIANCE_ALPHA_1_0_INIT_H
#include <stdlib.h>
#include <stdio.h>
#include "parameters_io.h"
#include "array.h"
#include "mpi_utils.h"
#include "fftw_utils.h"
#include "hdf_utils.h"
#include "fields.h"
#include "equation.h"
#include "space_config.h"
#include "solver.h"
#include "variables.h"

enum electromagnetic{ELECTROSTATIC,ELECTROMAGNETIC};
enum adiabatic{NONADIABATIC,ADIABATIC};
enum initial{RANDOM, FROMFILE};

void init_init(char *filename);

void init_printParameters();
void init_initEnums();
void init_conditions(COMPLEX *data);

extern enum adiabatic kinetic;
extern enum electromagnetic systemType;
extern enum initial initialConditions;
#endif //ALLIANCE_ALPHA_1_0_INIT_H
