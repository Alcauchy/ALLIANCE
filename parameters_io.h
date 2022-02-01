//
// Created by alcauchy on 27/12/2021.
//

#ifndef ALLIANCE_ALPHA_1_0_PARAMETERS_IO_C_H
#define ALLIANCE_ALPHA_1_0_PARAMETERS_IO_C_H

#include <stdio.h> // standard C lib
#include <stdlib.h> //Standard C lib
#include <string.h>
#include "mpi_utils.h"
#include "array.h"

enum DEALIASING {ALIASED, TWOTHIRDS};

void read_parameters(char *filename);

struct system_param {
    size_t nkx;
    size_t nky;
    size_t nkz;
    size_t nm;
    size_t nl;
    size_t ns;
    size_t nz;
    int dealiasing;

};


#endif //ALLIANCE_ALPHA_1_0_PARAMETERS_IO_C_H
