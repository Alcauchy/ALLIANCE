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
#include "fields.h"
#include "equation.h"
#include "space_config.h"
#include "solver.h"

enum electromagnetic{ELECTROSTATIC,ELECTROMAGNETIC};
enum adiabatic{NONADIABATIC,ADIABATIC};

void init_init(char *filename);

void init_printParameters();
void init_initEnums();


#endif //ALLIANCE_ALPHA_1_0_INIT_H
