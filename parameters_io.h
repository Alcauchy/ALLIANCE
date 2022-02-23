//
// Created by alcauchy on 27/12/2021.
//

#ifndef ALLIANCE_ALPHA_1_0_PARAMETERS_IO_C_H
#define ALLIANCE_ALPHA_1_0_PARAMETERS_IO_C_H

#include <stdio.h> // standard C lib
#include <stdlib.h> //Standard C lib
#include <string.h>
#include "hdf_utils.h"
#include "mpi_utils.h"
#include "array.h"

enum DEALIASING {ALIASED, TWOTHIRDS};

void read_parameters(char *filename);
void read_parametersFromFile(char *filename);

struct system_param {
    size_t nkx;
    size_t nky;
    size_t nkz;
    size_t nm;
    size_t nl;
    size_t ns;
    size_t nz;
    int dealiasing;
    int electromagnetic;
    int adiabatic;
    int initial;
    char from_simulationName[128];
    double beta;
    double *mass;
    double *density;
    double *temperature;
    double *charge;

    double dt;
    int Nt;

    int compute_k;
    int compute_m;
    int k_shells;
    int compute_k_every;
    int compute_m_every;
    int save_kSpec_step;
    int save_mSpec_step;
    int save_energy_step;
    int save_field_step;
    int checkpoints;
    int save_checkpoint_step;
    int save_distrib_step;
    int save_kSpec;
    int save_mSpec;
    int save_energy;
    int save_field;
    int save_checkpoint;
    int save_distrib;
    double firstShell;
    double lastShell;
};

extern struct system_param parameters;

#endif //ALLIANCE_ALPHA_1_0_PARAMETERS_IO_C_H
