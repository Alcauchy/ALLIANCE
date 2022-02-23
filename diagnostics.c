//
// Created by alcauchy on 20/01/2022.
//
#include "diagnostics.h"
#include "parameters_io.h"

#define TO_ROOT 0
#define BUFFER_SIZE 1

double *diag_kSpec = 0;
double *diag_mSpec = 0;
double *diag_shells = 0;
double diag_freeEnergy;

void diag_computeSpectra(const COMPLEX *g, const COMPLEX *h, int timestep) {
    if (parameters.compute_k && timestep % parameters.compute_k_every == 0) {
        diag_computeKSpectrum(g, h, diag_kSpec);
    }
    if (parameters.compute_m && timestep % parameters.compute_m_every == 0) {
        diag_computeMSpectrum(g, h, diag_mSpec);
    }
};

void diag_initSpec() {
    if (parameters.compute_k) {
        diag_kSpec = malloc(parameters.k_shells * sizeof(*diag_kSpec));
        diag_getShells();
    }
    if (parameters.compute_m) {
        if (mpi_my_row_rank == 0)
        {
            diag_mSpec = malloc(array_local_size.nm * sizeof(*diag_mSpec));
        }
        else
        {
            diag_mSpec = 0;
        }

    }
};

void diag_computeFreeEnergy(COMPLEX *g, COMPLEX *h, int it) {
    if (parameters.save_energy && it % parameters.save_energy_step == 0) {
        COMPLEX sum = 0;
        COMPLEX freeEnergy = 0;
        for (size_t i = 0; i < array_local_size.total_comp; i++) {
            sum += g[i] * conj(h[i]);
        }
        MPI_Reduce(&sum, &freeEnergy, BUFFER_SIZE, MPI_DOUBLE_COMPLEX, MPI_SUM, TO_ROOT, MPI_COMM_WORLD);
        diag_freeEnergy = creal(freeEnergy);
    }

};

void diag_computeKSpectrum(const COMPLEX *g, const COMPLEX *h, double *spec) {
    COMPLEX *sum = malloc(parameters.k_shells * sizeof(*sum));
    COMPLEX *buf = malloc(parameters.k_shells * sizeof(*buf));
    double *norm = malloc(parameters.k_shells * sizeof(*norm));
    double *total_norm = malloc(parameters.k_shells * sizeof(*total_norm));
    size_t ind2D;
    size_t ind6D;

    for (size_t ishell = 1; ishell < parameters.k_shells + 1; ishell++) {
        sum[ishell - 1] = 0.j;
        norm[ishell - 1] = 0;
        for (size_t ix = 0; ix < array_local_size.nkx; ix++) {
            for (size_t iy = 0; iy < array_local_size.nky; iy++) {
                ind2D = ix * array_local_size.nky + iy;
                if (diag_shells[ishell - 1] < space_kPerp2[ind2D] && diag_shells[ishell] >= space_kPerp2[ind2D]) {
                    for (size_t iz = 0; iz < array_local_size.nkz; iz++) {
                        for (size_t im = 0; im < array_local_size.nm; im++) {
                            for (size_t il = 0; il < array_local_size.nl; il++) {
                                for (size_t is = 0; is < array_local_size.ns; is++) {
                                    norm[ishell - 1] += 1.;
                                    ind6D = get_flat_c(is, il, im, ix, iy, iz);
                                    sum[ishell - 1] += g[ind6D] * conj(h[ind6D]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    MPI_Reduce(sum, buf, parameters.k_shells, MPI_DOUBLE_COMPLEX, MPI_SUM, TO_ROOT, MPI_COMM_WORLD);
    MPI_Reduce(norm, total_norm, parameters.k_shells, MPI_DOUBLE, MPI_SUM, TO_ROOT, MPI_COMM_WORLD);
    if (mpi_my_rank == TO_ROOT)
    {
        for (size_t i = 0; i < parameters.k_shells; i++){
            if (total_norm[i])
            {
                diag_kSpec[i] = creal(buf[i]/total_norm[i]);
            }
            else
            {
                diag_kSpec[i] = 0;
            }

        }
    }
};

void diag_computeMSpectrum(const COMPLEX *g, const COMPLEX *h, double *spec) {
    COMPLEX *sum = malloc(array_local_size.nm * sizeof(*sum));
    COMPLEX *buf = malloc(array_local_size.nm * sizeof(*buf));
    size_t ind6D;
    for(size_t ix = 0; ix < array_local_size.nkx; ix++){
        for(size_t iy = 0; iy < array_local_size.nky; iy++){
            for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                for(size_t im = 0; im < array_local_size.nm; im++){
                    sum[im] = 0;
                    for(size_t il = 0; il < array_local_size.nl; il++){
                        for(size_t is = 0; is < array_local_size.ns; is++){
                            ind6D = get_flat_c(is,il,im,ix,iy,iz);
                            sum[im] += g[ind6D] * conj(h[ind6D]);
                        }
                    }
                }
            }
        }
    }
    MPI_Reduce(sum, buf, array_local_size.nm, MPI_DOUBLE_COMPLEX, MPI_SUM, TO_ROOT, mpi_row_comm);
    if(mpi_my_row_rank == TO_ROOT)
    {
        for (size_t i = 0; i < array_local_size.nm; i++){
            diag_mSpec[i] = creal(buf[i]);
        }
    }



};

void diag_getShells() {
    diag_shells = malloc((parameters.k_shells + 1) * sizeof(*diag_shells));
    for (size_t i = 0; i < parameters.k_shells + 1; i++) {
        diag_shells[i] = (parameters.lastShell - parameters.firstShell) / parameters.k_shells * i;
        printf("[MPI process %d] shell[%d] = %f\n", mpi_my_rank, i, diag_shells[i]);
    }

}

double diag_computeFreeEnergyFields(COMPLEX *g, COMPLEX *fields) {};

