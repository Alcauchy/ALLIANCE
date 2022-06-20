/**************************************
* @file diagnostics.c
*   \brief diagnostics module
*
*   different diagnostic tools are gathered in this module
***************************************/
////////////////////////////////////////////////////////////////////////////////
// 20/01/2022 created by Gene Gorbunov
//                                   DIAGNOSTICS
//
// diag_computeSpectra
// diag_computeFreeEnergy
// diag_computeFreeEnergyFields
// diag_initSpec
// diag_computeKSpectrum
// diag_computeMSpectrum
// diag_initSpec
// diag_getShells
// diag_compute
//
// VERSION 1.0
////////////////////////////////////////////////////////////////////////////////
#include "diagnostics.h"
#include "parameters_io.h"

#define TO_ROOT 0
#define BUFFER_SIZE 1
/**\var double *diag_kSpec
 * \brief used to store free energy k spectra*/
double *diag_kSpec = 0;

/**\var double *diag_mSpec
 * \brief used to store free energy m spectra*/
double *diag_mSpec = 0;

/**\var double *diag_shells
 * \brief used to store positions of k shells  required to compute k spectra*/
double *diag_shells = 0;

/**\var double diag_freeEnergy
 * \brief free energy*/
double diag_freeEnergy;

/***************************************
 * \fn void diag_computeSpectra(const COMPLEX *g, const COMPLEX *h, int timestep)
 * \brief general function to compute k or m spectra
 * \param g: gyrokinetic distribution function
 * \param h: distribution function
 * \param timestep: current time step
 *
 * function computes spectra at timestep as given in parameter file.
 * \f$k_{\perp}\f$ spectra is computed using #diag_computeKSpectrum, and
 * m spectra is computed using #diag_computeMSpectrum
 ***************************************/
void diag_computeSpectra(const COMPLEX *g, const COMPLEX *h, int timestep) {
    if (parameters.compute_k && timestep % parameters.iter_diagnostics == 0) {
        diag_computeKSpectrum(g, h, diag_kSpec);
    }
    if (parameters.compute_m && timestep % parameters.iter_diagnostics == 0) {
        diag_computeMSpectrum(g, h, diag_mSpec);
    }
};

/***************************************
 * \fn void diag_initSpec()
 * \brief initialize spectra computation
 *
 * Prepares free energy spectra computation. For spectra in k:
 * allocates diag_kSpec array used to store k spectra.
 * Allocates diag_shells array and fills it with
 * shell positions \f$k^{shells}\f$, used for binning of wave vectors
 * when computing \f$k_\perp\f$ spectra.
 * For spectra in m:
 * allocates diag_mSpec array used to store m spectra.
 * Called in #init_start function
 ***************************************/
void diag_initSpec() {
    if (parameters.compute_k) {
        diag_kSpec = malloc(parameters.k_shells * sizeof(*diag_kSpec));
        diag_getShells();
    }
    if (parameters.compute_m) {
        if (mpi_my_row_rank == 0) {
            diag_mSpec = malloc(array_local_size.nm * sizeof(*diag_mSpec));
        } else {
            diag_mSpec = 0;
        }

    }
};

/***************************************
 * \fn void diag_computeFreeEnergy(COMPLEX *g, COMPLEX *h)
 * \brief compute free energy
 * \param g: modified gyrokinetic distribution function
 * \param h: gyrokintic distribution function
 *
 * computes free energy as \f$ W = 2. \Re(\sum_{k_x,k_y, k_z>0,m,l,s} g * \bar{h})\f$,
 * taking into account reality condition.
 ***************************************/
void diag_computeFreeEnergy(COMPLEX *g, COMPLEX *h) {
    COMPLEX sum = 0;
    COMPLEX freeEnergy = 0;
    for (size_t i = 0; i < array_local_size.total_comp; i++) {
        sum += g[i] * conj(h[i]);
    }
    MPI_Reduce(&sum, &freeEnergy, BUFFER_SIZE, MPI_C_DOUBLE_COMPLEX, MPI_SUM, TO_ROOT, MPI_COMM_WORLD);
    diag_freeEnergy = 2. * creal(freeEnergy);

};

/***************************************
 * \fn void diag_computeKSpectrum(const COMPLEX *g, const COMPLEX *h, double *spec)
 * \param g: modified gyrokinetic distribution function
 * \param h: gyrokinetic distribution function
 * \param spec: spectra array
 *
 * computes free energy \f$k_{\perp}\f$ spectra
 * \f$ W(k^{shell}_i) = \frac{1}{N}\sum_{k^{shell}_{i-1}<|\mathbf{k}_{\perp}|<k^{shell}_{i}}\sum_{k_z,l,m,s} g \bar{h}\f$
 * where \f$N\f$ is a number of wave vectors between shells \f$k^{shell}_{i-1}\f$ and \f$k^{shell}_{i}\f$
 ***************************************/
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
    if (mpi_my_rank == TO_ROOT) {
        for (size_t i = 0; i < parameters.k_shells; i++) {
            if (total_norm[i]) {
                diag_kSpec[i] = creal(buf[i] / total_norm[i]);
            } else {
                diag_kSpec[i] = 0;
            }

        }
    }
};

/***************************************
 * \fn void diag_computeMSpectrum(const COMPLEX *g, const COMPLEX *h, double *spec)
 * \brief computes free energy spectra in m space
 * \param g: modified gyrokinetic distribution function
 * \param h: gyrokinetic distribution function
 * \param spec: spectra array
 *
 * computes free energy m spectra as \f$ W(m) = \sum_{k_x,k_y,k_z,l,s} g \bar{h} \f$
 ***************************************/
void diag_computeMSpectrum(const COMPLEX *g, const COMPLEX *h, double *spec) {
    COMPLEX *sum = calloc(array_local_size.nm , sizeof(*sum));
    COMPLEX *buf = malloc(array_local_size.nm * sizeof(*buf));
    size_t ind6D;
    for (size_t ix = 0; ix < array_local_size.nkx; ix++) {
        for (size_t iy = 0; iy < array_local_size.nky; iy++) {
            for (size_t iz = 0; iz < array_local_size.nkz; iz++) {
                for (size_t im = 0; im < array_local_size.nm; im++) {
                    for (size_t il = 0; il < array_local_size.nl; il++) {
                        for (size_t is = 0; is < array_local_size.ns; is++) {
                            ind6D = get_flat_c(is, il, im, ix, iy, iz);
                            sum[im] += g[ind6D] * conj(h[ind6D]);
                        }
                    }
                }
            }
        }
    }
    MPI_Reduce(sum, buf, array_local_size.nm, MPI_DOUBLE_COMPLEX, MPI_SUM, TO_ROOT, mpi_row_comm);
    if (mpi_my_row_rank == TO_ROOT) {
        for (size_t i = 0; i < array_local_size.nm; i++) {
            diag_mSpec[i] = creal(buf[i]);
        }
    }


};

/***************************************
 * \fn diag_getShells()
 * \brief computes shells from parameters
 *
 * computes positions of k_shells in between
 * last_shell and first_shell as provided
 * by user in parameter file. Position of \f$i^{th}\f$ shell is computed as
 * \f$k^{shell}_i = (last\_shell - first\_shell)/(k\_shells) \cdot i\f$
 ***************************************/
void diag_getShells() {
    diag_shells = malloc((parameters.k_shells + 1) * sizeof(*diag_shells));
    for (size_t i = 0; i < parameters.k_shells + 1; i++) {
        diag_shells[i] = (parameters.lastShell - parameters.firstShell) / parameters.k_shells * i;
        //printf("[MPI process %d] shell[%d] = %f\n", mpi_my_rank, i, diag_shells[i]);
    }

}

/***************************************
 * \fn diag_computeFreeEnergyFields(COMPLEX *g, COMPLEX *fields)
 * \brief to be done later
 * \param g:
 * \param fields:
 *
 * computes free energy from the fields and distribution function.
 ***************************************/
double diag_computeFreeEnergyFields(COMPLEX *g, COMPLEX *fields) {};

/***************************************
 * \fn void diag_compute(COMPLEX *g, COMPLEX *h, int timestep)
 * \brief computes all diagnostics
 * \param g: modified distribution function
 * \param h: distribution function
 * \param iter: current time step
 ***************************************/
void diag_compute(COMPLEX *g, COMPLEX *h, int timestep) {
    diag_computeFreeEnergy(g, h);
    hdf_saveEnergy(timestep);
    if (parameters.compute_k) {
        diag_computeKSpectrum(g, h, diag_kSpec);
        hdf_saveKSpec(timestep);
    }
    if (parameters.compute_m) {
        diag_computeMSpectrum(g, h, diag_mSpec);
        hdf_saveMSpec(timestep);
    }
}