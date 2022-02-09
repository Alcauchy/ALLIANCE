//
// Created by alcauchy on 20/01/2022.
//
#include "diagnostics.h"
#define TO_ROOT 0
#define BUFFER_SIZE 1
void diag_computeSpectra(COMPLEX *data, double *spec){};


double diag_computeFreeEnergy(COMPLEX *g,COMPLEX *h){
    COMPLEX sum = 0;
    COMPLEX freeEnergy = 0;
    for (size_t i = 0; i < array_local_size.total_comp; i++)
    {
        sum += g[i] * conj(h[i]);
    }
    MPI_Reduce(&sum, &freeEnergy, BUFFER_SIZE, MPI_DOUBLE_COMPLEX, MPI_SUM,TO_ROOT,MPI_COMM_WORLD);
    return creal(freeEnergy);
};

double diag_computeFreeEnergyFields(COMPLEX *g, COMPLEX *fields){};

