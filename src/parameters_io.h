#ifndef ALLIANCE_ALPHA_1_0_PARAMETERS_IO_C_H
#define ALLIANCE_ALPHA_1_0_PARAMETERS_IO_C_H

#include <stdio.h> // standard C lib
#include <stdlib.h> //Standard C lib
#include <string.h>
#include "utils_hdf.h"
#include "utils_mpi.h"
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
    double Lx;
    double Ly;
    double Lz;

    double mu_k;
    double mu_m;

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
    double linDt;
    double dissipDt;
    int iter_dt;
    int Nt;

    double forceKmin;
    double forceKmax;
    double forcePower;

    int compute_k;
    int compute_m;
    int k_shells;
    int iter_diagnostics;
    int iter_EMfield;
    int checkpoints;
    int iter_checkpoint;
    int iter_distribution;
    int save_diagnostics;
    int save_EMfield;
    int save_checkpoint;
    int save_distrib;
    double firstShell;
    double lastShell;

    int nproc_m;
    int nproc_k;
};

extern struct system_param parameters;

#endif //ALLIANCE_ALPHA_1_0_PARAMETERS_IO_C_H
