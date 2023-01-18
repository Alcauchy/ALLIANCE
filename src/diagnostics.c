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
// diag_initNonlinearFLux
// diag_computeKSpectrum
// diag_computeMSpectrum
// diag_initSpec
// diag_getShells
// diag_compute
// diag_filterK
// diag_computeNonlinearFlux
// VERSION 1.0
////////////////////////////////////////////////////////////////////////////////
#include "diagnostics.h"
#include "parameters_io.h"

#define TO_ROOT 0
#define BUFFER_SIZE 1
#define IO_RANK 0
/**\var double *diag_kSpec
 * \brief used to store free energy k spectra*/
double *diag_kSpec = 0;

/**\var double *diag_kSpecPhi
 * \brief used to store phi energy k spectra*/
double *diag_kSpecPhi = 0;

/**\var double *diag_kSpecBperp
 * \brief used to store B_perp energy k spectra*/
double *diag_kSpecBperp = 0;

/**\var double *diag_kSpecBpar
 * \brief used to store B_par energy k spectra*/
double *diag_kSpecBpar = 0;

/**\var double *diag_kSpecH
 * \brief used to store h energy k spectra*/
double *diag_kSpecH = 0;

/**\var double *diag_LDis
 * \brief dissipation length scale*/
double diag_LDis;

/**\var double *diag_LInt
 * \brief integral length scale*/
double diag_LInt;

/**\var double *diag_mSpec
 * \brief used to store free energy m spectra*/
double *diag_mSpec = 0;

/**\var double *diag_shells
 * \brief used to store positions of k shells  required to compute k spectra*/
double *diag_shells = 0;

/**\var double diag_shellCentres
 * \brief centers of shells*/
double *diag_shellCentres = 0;

/**\var double diag_shellCentres
 * \brief centers of shells*/
double *diag_shellNorm = 0;

/**\var double diag_nonlinearFlux
 * \brief used to store nonlinear flux */
double *diag_nonlinearFlux = 0;

/**\var double diag_nonlinearFluxInverse
 * \brief used to store inverse nonlinear flux */
double *diag_nonlinearFluxInverse = 0;

/**\var double diag_nonlinearNorm
* \brief normalization factor for nonlinear flux */
double *diag_nonlinearNorm = 0;

/**\var double diag_freeEnergy
 * \brief free energy*/
double diag_freeEnergy;

/**\var double diag_freeEnergy
 * \brief free energy at initial timestep*/
double diag_free_energy0;

/**\var double diag_fib
 * \brief root of golden ratio scale shells*/
double diag_sqrtGoldenRatio;


/**\var double diag_fib
 * \brief number of shells*/
int diag_numOfShells;

/**\var double diag_fib
 * \brief number of shells bounds*/
int diag_numOfShellBounds;

/**\var double diag_energyH
 * \brief h contribution to free energy*/
double diag_energyH;

/**\var double diag_energyPhi
 * \brief phi contribution to free energy*/
double diag_energyPhi;

/**\var double diag_energyBperp
 * \brief Bperp contribution to free energy*/
double diag_energyBperp;

/**\var double diag_energyBpar
 * \brief Bpar contribution to free energy*/
double diag_energyBpar;

/**\var double diag_energyTotal
 * \brief total amount of free energy*/
double diag_energyTotal;

/**\var double diag_injected
 * \brief diag_injected amount of energy due to forcing*/
double diag_injected;

/**\var double diag_dissipated
 * \brief diag_dissipated amount of energy*/
double diag_dissipated;

/**\var double diag_dissipated_kPerp
 * \brief diag_dissipated amount of energy in k_perp*/
double diag_dissipated_kPerp;

/**\var double diag_dissipated_kZ
 * \brief diag_dissipated amount of energy in k_z*/
double diag_dissipated_kZ;

/**\var double diag_dissipated_m
 * \brief diag_dissipated amount of energy in m*/
double diag_dissipated_m;

/**\var double diag_etakmax
 * \brief eta*kmax*/
double diag_etakmax;

/**\var double *diag_realityConditionNorm
 * \brief equals to 1 when z = 0, and 2 otherwise*/
double *diag_MM;

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
        diag_computeFieldSpectrum();
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
    diag_MM = malloc((array_local_size.nkz) * sizeof(*diag_MM));
    for (size_t iz = 0; iz < array_local_size.nkz; iz++ ) {
        if (iz == 0 || iz == array_local_size.nkz - 1) {
            diag_MM[iz] = 1.0;
        }
        else {
            diag_MM[iz] = 2.0;
        }
    }
        if (parameters.compute_k || parameters.compute_nonlinear){
            diag_getShells();
            diag_kSpec = malloc((diag_numOfShells) * sizeof(*diag_kSpec));
            if(parameters.compute_k){
                if (systemType == ELECTROMAGNETIC){
                    diag_kSpecPhi = calloc((diag_numOfShells), sizeof(*diag_kSpecPhi));
                    diag_kSpecBperp = calloc((diag_numOfShells), sizeof(*diag_kSpecBperp));
                    diag_kSpecBpar = calloc((diag_numOfShells), sizeof(*diag_kSpecBpar));
                    diag_kSpecH = calloc((diag_numOfShells) * array_local_size.nm, sizeof(*diag_kSpecH));
                }
                if (systemType == ELECTROSTATIC){
                    diag_kSpecPhi = calloc((diag_numOfShells), sizeof(*diag_kSpecPhi));
                    diag_kSpecH = calloc((diag_numOfShells) * array_local_size.nm, sizeof(*diag_kSpecH));
                }
            }
            if(parameters.compute_nonlinear){
                diag_nonlinearFlux = calloc((diag_numOfShells),sizeof(*diag_nonlinearFlux));
                diag_nonlinearFluxInverse = calloc((diag_numOfShells),sizeof(*diag_nonlinearFluxInverse));
            }
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
    size_t ind6D;
    size_t indJ0;
    for (size_t ix = 0; ix < array_local_size.nkx; ix++) {
        for (size_t iy = 0; iy < array_local_size.nky; iy++) {
            for (size_t iz = 0; iz < array_local_size.nkz; iz++) {
                for (size_t im = 0; im < array_local_size.nm; im++) {
                    for (size_t il = 0; il < array_local_size.nl; il++) {
                        for (size_t is = 0; is < array_local_size.ns; is++) {
                            ind6D = get_flat_c(is, il, im, ix, iy, iz);
                            indJ0 = var_getJIndex(ix,iy,is);
                            sum += 0.5 * var_var.T[is] * var_var.n[is] * g[ind6D] * conj(h[ind6D]) * diag_MM[iz];
                        }
                    }
                }
            }
        }
    }
    MPI_Reduce(&sum, &freeEnergy, BUFFER_SIZE, MPI_C_DOUBLE_COMPLEX, MPI_SUM, TO_ROOT, MPI_COMM_WORLD);
    diag_freeEnergy = cabs(freeEnergy);
    //printf("computed as old = %f\n",diag_freeEnergy);

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
    COMPLEX *sum = malloc((diag_numOfShells) * sizeof(*sum));
    COMPLEX *buf = malloc((diag_numOfShells) * sizeof(*buf));
    double *norm = malloc((diag_numOfShells) * sizeof(*norm));
    double *total_norm = malloc((diag_numOfShells) * sizeof(*total_norm));
    size_t ind2D;
    size_t ind6D;
    for (size_t ishell = 1; ishell < diag_numOfShellBounds; ishell++) {
        sum[ishell - 1] = 0.j;
        for (size_t ix = 0; ix < array_local_size.nkx; ix++) {
            for (size_t iy = 0; iy < array_local_size.nky; iy++) {
                ind2D = ix * array_local_size.nky + iy;
                if (space_kPerp[ind2D] > diag_shells[ishell - 1]  && space_kPerp[ind2D] <= diag_shells[ishell]) {
                    for (size_t iz = 0; iz < array_local_size.nkz; iz++) {
                        for (size_t im = 0; im < array_local_size.nm; im++) {
                            for (size_t il = 0; il < array_local_size.nl; il++) {
                                for (size_t is = 0; is < array_local_size.ns; is++) {
                                    norm[ishell - 1] += 1.;
                                    ind6D = get_flat_c(is, il, im, ix, iy, iz);
                                    sum[ishell - 1] += g[ind6D] * conj(h[ind6D]);
                                    //sum[ishell - 1] /= diag_shellNorm[ishell - 1];
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    MPI_Reduce(sum, buf, diag_numOfShells, MPI_DOUBLE_COMPLEX, MPI_SUM, TO_ROOT, MPI_COMM_WORLD);
    diag_LDis = 0;
    diag_LInt = 0;
    double divisorLDis = 0;
    double divisorLInt = 0;
    if (mpi_my_rank == TO_ROOT) {
        for (size_t i = 0; i < diag_numOfShells; i++) {
            diag_kSpec[i] = creal(buf[i])/diag_shellNorm[i];
            //computing sum(k^(2 alpha) * E(k))
            divisorLDis += pow(diag_shellCentres[i],2 * var_var.mu_k) * diag_kSpec[i];
            //computing sum (E(k))
            divisorLInt += diag_kSpec[i];
            //computing sum(k^(2 alpha - 1) * E(k))
            diag_LDis += pow(diag_shellCentres[i],2 * var_var.mu_k - 1) * diag_kSpec[i];
            //computing sum(k^(-1) * E(k))
            diag_LInt += diag_kSpec[i]/diag_shellCentres[i];
        }
        printf("%f\n", diag_LDis);
        printf("%f\n", diag_LInt);
        diag_LDis *= 2 * M_PI / divisorLDis;
        diag_LInt *= 2 * M_PI / divisorLInt;
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
    if (mpi_my_row_rank == TO_ROOT){
        for (size_t i = 0; i < array_local_size.nm; i++){
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
 * also computes norms for nonlinear flux
 ***************************************/
void diag_getShells() {
    diag_sqrtGoldenRatio = sqrt((1. + sqrt(5.)) * 0.5);
    double shell1;
    double k0;
    if (spectrumType == SHELL){
        k0 = (space_dKx > space_dKy) ? space_dKx : space_dKy;
        if (parameters.firstShell == 0){
            shell1 = 4. * k0;
        }
        else{
            shell1 = parameters.firstShell;
        }
        diag_numOfShells = (int) (log(parameters.lastShell/shell1)/log(diag_sqrtGoldenRatio) + 1) + 1;
        diag_numOfShellBounds = diag_numOfShells + 1;
        if(diag_numOfShells < 3) {
            if(mpi_my_rank == 0) {
                printf("Wrong shell boundaries with k_first = %f and k_last = %f\n"
                       "EXITING \n",
                       parameters.firstShell,
                       parameters.lastShell);
            }
            exit(1);
        }
        if(mpi_my_rank == 0) {
            printf("number of shells = %d  shell 1 = %f\n",
                   diag_numOfShells,
                   shell1);
        }
        diag_shells = malloc((diag_numOfShellBounds) * sizeof(*diag_shells));

        // shell boundaries
        diag_shells[0] = 0.;
        diag_shells[1] = shell1;
        for (size_t i = 2; i < diag_numOfShellBounds; i++) {
            diag_shells[i] = pow(diag_sqrtGoldenRatio,i - 1) * shell1;
        }
        diag_shells[diag_numOfShellBounds - 1] = parameters.lastShell;

        //shell centres and normalizations
        diag_shellCentres = malloc((diag_numOfShells) * sizeof(*diag_shellCentres));
        diag_shellNorm = calloc((diag_numOfShells), sizeof(*diag_shellNorm));
        diag_nonlinearNorm = calloc((diag_numOfShells), sizeof(*diag_nonlinearNorm));

        for (size_t ishell = 1; ishell < (diag_numOfShellBounds); ishell++){
            diag_shellCentres[ishell - 1] = 0.5 * (diag_shells[ishell - 1] + diag_shells[ishell ]);
            for(size_t ix = 0; ix < array_local_size.nkx; ix++){
                for(size_t iy = 0; iy < array_local_size.nky; iy++){
                    for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                        size_t ind2D = ix * array_local_size.nky + iy;
                        if (diag_shells[ishell - 1] < space_kPerp[ind2D] && diag_shells[ishell] >= space_kPerp[ind2D]){
                            if (iz == 0 || iz == array_local_size.nkz - 1){
                                diag_shellNorm[ishell - 1] += 1.;
                            }
                            else{
                                diag_shellNorm[ishell - 1] += 2.;
                            }
                        }
                        if (space_kPerp[ind2D] < diag_shells[ishell]){
                            if (iz == 0 || iz == array_local_size.nkz - 1){
                                diag_nonlinearNorm[ishell - 1] += 1.;
                            }
                            else{
                                diag_nonlinearNorm[ishell - 1] += 2.;
                            }
                        }
                    }
                }
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, diag_shellNorm,diag_numOfShells, MPI_DOUBLE, MPI_SUM, mpi_row_comm);
        MPI_Allreduce(MPI_IN_PLACE, diag_nonlinearNorm,diag_numOfShells, MPI_DOUBLE, MPI_SUM, mpi_row_comm);
    }

    if (spectrumType == UNIT){
        k0 = parameters.unitK;
        double maxK = ((array_global_size.nkx / 2 + 1) * space_dKx < (array_global_size.nky / 2 + 1) * space_dKy) ? (array_global_size.nkx / 2 + 1) * space_dKx : (array_global_size.nky / 2 + 1) * space_dKy;
        diag_numOfShells = (int) maxK/k0;
        diag_numOfShellBounds = diag_numOfShells + 1;
        if(mpi_my_rank == 0) {
            printf("number of k for spectrum = %d\n", diag_numOfShells);
        }

        //shell boundaries
        diag_shells = malloc((diag_numOfShellBounds) * sizeof(*diag_shells));
        for (size_t i = 0; i < diag_numOfShellBounds; i++) {
            diag_shells[i] = k0 * i;
        }

        //shell normalizations and centres
        diag_shellCentres = malloc((diag_numOfShells) * sizeof(*diag_shellCentres));
        diag_shellNorm = calloc((diag_numOfShells), sizeof(*diag_shellNorm));
        diag_nonlinearNorm = calloc((diag_numOfShells), sizeof(*diag_nonlinearNorm));
        for (size_t ishell = 1; ishell < (diag_numOfShellBounds); ishell++){
            diag_shellCentres[ishell - 1] = 0.5 * (diag_shells[ishell - 1] + diag_shells[ishell ]);
            for(size_t ix = 0; ix < array_local_size.nkx; ix++){
                for(size_t iy = 0; iy < array_local_size.nky; iy++){
                    for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                        size_t ind2D = ix * array_local_size.nky + iy;
                        if (diag_shells[ishell - 1] < space_kPerp[ind2D] && diag_shells[ishell] >= space_kPerp[ind2D]){
                            if (iz == 0 || iz == array_local_size.nkz - 1){
                                diag_shellNorm[ishell - 1] += 1.;
                            }
                            else{
                                diag_shellNorm[ishell - 1] += 2.;
                            }
                        }
                        if (space_kPerp[ind2D] <= diag_shells[ishell]){
                            if (iz == 0 || iz == array_local_size.nkz - 1){
                                diag_nonlinearNorm[ishell - 1] += 1.;
                            }
                            else{
                                diag_nonlinearNorm[ishell - 1] += 2.;
                            }
                        }
                    }
                }
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, diag_shellNorm,diag_numOfShells, MPI_DOUBLE, MPI_SUM, mpi_row_comm);
        MPI_Allreduce(MPI_IN_PLACE, diag_nonlinearNorm,diag_numOfShells, MPI_DOUBLE, MPI_SUM, mpi_row_comm);
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
    if (parameters.save_diagnostics && (timestep % parameters.iter_diagnostics == 0 || timestep == -1)){
        // update fields and h first
        fields_sendF(g);
        fields_getFields(f00, f10, f01);
        fields_getChi();
        distrib_getH(h, g);
        //compute and save free energy
        diag_computeFreeEnergy(g, h);
        diag_computeEnergy(h);
        //diag_freeEnergy = diag_energyTotal;
        if (timestep == 0) diag_free_energy0 = diag_energyTotal;
        hdf_saveEnergy(timestep);
        hdf_saveForcing();
        //compute spectra
        if (parameters.compute_k) {
            diag_computeKSpectrum(g, h, diag_kSpec);
            diag_computeFieldSpectrum();
            diag_computeHSpectrum(h);
            hdf_saveKSpec(timestep);
        }
        if (parameters.compute_m) {
            diag_computeMSpectrum(g, h, diag_mSpec);
            hdf_saveMSpec(timestep);
        }
    }
    if(parameters.compute_nonlinear && (timestep % parameters.compute_nonlinear == 0 || timestep == -1)){
        // update fields and h first
        fields_sendF(g);
        fields_getFields(f00, f10, f01);
        fields_getChi();
        distrib_getH(h, g);
        //compute and save free energy
        diag_computeFreeEnergy(g, h);
        diag_computeEnergy(h);
        diag_computeNonlinearFlux(h);
        hdf_saveNonlinearFlux(timestep);

    }
}

/***************************************
 * \fn diag_computeFieldSpectrum()
 *
 * computes field energy \f$k_{\perp}\f$ spectra
 * \f$ W(k^{shell}_i) = \frac{1}{N}\sum_{k^{shell}_{i-1}<|\mathbf{k}_{\perp}|<k^{shell}_{i}}\sum_{k_z,l,m,s} g \bar{h}\f$
 * where \f$N\f$ is a number of wave vectors between shells \f$k^{shell}_{i-1}\f$ and \f$k^{shell}_{i}\f$
 ***************************************/
void diag_computeFieldSpectrum() {
    size_t ind3D;
    size_t ind2D;
    switch(systemType){
        case ELECTROSTATIC:
            for(size_t ishell = 0; ishell < diag_numOfShells; ishell++){
                diag_kSpecPhi[ishell] = 0.;
                for (size_t ix = 0; ix < array_local_size.nkx; ix++){
                    for (size_t iy = 0; ix < array_local_size.nky; iy++){
                        for (size_t iz = 0; ix < array_local_size.nkz; iz++){
                            ind3D = ix * array_local_size.nky * array_local_size.nkz +
                                    iy * array_local_size.nkz +
                                    iz;
                            ind2D = ix * array_local_size.nky +
                                    iy;
                            if (space_kPerp[ind2D] > diag_shells[ishell]  && space_kPerp[ind2D] <= diag_shells[ishell + 1]){
                                for (size_t is = 0; is < array_local_size.ns; is++){
                                    diag_kSpecPhi[ishell] += 0.5 * var_var.q[is] * var_var.q[is] * var_var.n[is] / var_var.T[is]
                                                             * cabs(fields_fields.phi[ind3D]) * cabs(fields_fields.phi[ind3D]);
                                }
                            }
                        }
                    }
                }
                diag_kSpecPhi[ishell] /= diag_shellNorm[ishell];
            }
            MPI_Reduce(MPI_IN_PLACE, diag_kSpecPhi, diag_numOfShells, MPI_DOUBLE, MPI_SUM, TO_ROOT, MPI_COMM_WORLD);
            break;
        case ELECTROMAGNETIC:
            for(size_t ishell = 0; ishell < diag_numOfShells; ishell++){
                diag_kSpecPhi[ishell] = 0;
                diag_kSpecBpar[ishell] = 0;
                diag_kSpecPhi[ishell] = 0;
                for (size_t ix = 0; ix < array_local_size.nkx; ix++){
                    for (size_t iy = 0; iy < array_local_size.nky; iy++){
                        for (size_t iz = 0; iz < array_local_size.nkz; iz++){
                            ind3D = ix * array_local_size.nky * array_local_size.nkz +
                                    iy * array_local_size.nkz +
                                    iz;
                            ind2D = ix * array_local_size.nky +
                                    iy;
                            if (space_kPerp[ind2D] > diag_shells[ishell]  && space_kPerp[ind2D] <= diag_shells[ishell + 1]){
                                /*phi spectra*/
                                for (size_t is = 0; is < array_local_size.ns; is++){
                                    diag_kSpecPhi[ishell] += 0.5 * var_var.q[is] * var_var.q[is] * var_var.n[is] / var_var.T[is]
                                                                 * cabs(fields_fields.phi[ind3D]) * cabs(fields_fields.phi[ind3D]);
                                }
                                /* Bpar spectra */
                                diag_kSpecBpar[ishell] += cabs(fields_fields.B[ind3D]) * cabs(fields_fields.B[ind3D])/8./M_PI;
                                /* Bperp spectra */
                                diag_kSpecBperp[ishell] += space_kPerp2[ind2D] * cabs(fields_fields.A[ind3D]) * cabs(fields_fields.A[ind3D])/8./M_PI;
                            }
                        }
                    }
                }
                diag_kSpecPhi[ishell] /= diag_shellNorm[ishell];
                diag_kSpecBpar[ishell] /= diag_shellNorm[ishell];
                diag_kSpecBperp[ishell] /= diag_shellNorm[ishell];
            }
            if (mpi_my_row_rank == 0){
                MPI_Reduce(MPI_IN_PLACE, diag_kSpecPhi, diag_numOfShells, MPI_DOUBLE, MPI_SUM, TO_ROOT, mpi_row_comm);
                MPI_Reduce(MPI_IN_PLACE, diag_kSpecBpar, diag_numOfShells, MPI_DOUBLE, MPI_SUM, TO_ROOT, mpi_row_comm);
                MPI_Reduce(MPI_IN_PLACE, diag_kSpecBperp, diag_numOfShells, MPI_DOUBLE, MPI_SUM, TO_ROOT, mpi_row_comm);
            }
            else{
                MPI_Reduce(diag_kSpecPhi, diag_kSpecPhi, diag_numOfShells, MPI_DOUBLE, MPI_SUM, TO_ROOT, mpi_row_comm);
                MPI_Reduce(diag_kSpecBpar, diag_kSpecBpar, diag_numOfShells, MPI_DOUBLE, MPI_SUM, TO_ROOT, mpi_row_comm);
                MPI_Reduce(diag_kSpecBperp, diag_kSpecBperp, diag_numOfShells, MPI_DOUBLE, MPI_SUM, TO_ROOT, mpi_row_comm);
            }

            break;
        default:
            printf("ERROR COMPUTING FIELD SPECTRA! EXITING...\n");
            exit(1);
    }
};

/***************************************
 * \fn diag_computeHSpectrum(const COMPLEX *h)
 *
 * computes field energy \f$k_{\perp}\f$ spectra
 * \f$ W(k^{shell}_i) = \frac{1}{N}\sum_{k^{shell}_{i-1}<|\mathbf{k}_{\perp}|<k^{shell}_{i}}\sum_{k_z,l,m,s} g \bar{h}\f$
 * where \f$N\f$ is a number of wave vectors between shells \f$k^{shell}_{i-1}\f$ and \f$k^{shell}_{i}\f$
 ***************************************/
void diag_computeHSpectrum(const COMPLEX *h) {
    size_t ind6D;
    size_t ind2Dm;
    size_t ind2D;
    for (size_t ishell = 0; ishell < diag_numOfShells; ishell++){
        for(size_t im = 0; im < array_local_size.nm; im ++){
            ind2Dm = ishell * array_local_size.nm + im;
            diag_kSpecH[ind2Dm] = 0;
            for(size_t ix = 0; ix < array_local_size.nkx; ix ++){
                for(size_t iy = 0; iy < array_local_size.nky; iy ++){
                    ind2D = ix * array_local_size.nky + iy;
                    if (space_kPerp[ind2D] > diag_shells[ishell]  && space_kPerp[ind2D] <= diag_shells[ishell + 1]){
                        for(size_t iz = 0; iz < array_local_size.nkz; iz ++){
                            for(size_t il = 0; il < array_local_size.nl; il ++){
                                for(size_t is = 0; is < array_local_size.ns; is ++){
                                    ind6D = get_flat_c(is,il,im,ix,iy,iz);
                                    diag_kSpecH[ind2Dm] += 0.5 * var_var.T[is] * var_var.n[is] * cabs(h[ind6D]) * cabs(h[ind6D]) / diag_shellNorm[ishell];;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    if(mpi_my_row_rank == 0){
        MPI_Reduce(MPI_IN_PLACE, diag_kSpecH, (diag_numOfShells) * array_local_size.nm, MPI_DOUBLE, MPI_SUM, TO_ROOT, mpi_row_comm);
    }
    else{
        MPI_Reduce(diag_kSpecH, diag_kSpecH, (diag_numOfShells)* array_local_size.nm, MPI_DOUBLE, MPI_SUM, TO_ROOT, mpi_row_comm);
    }
}

/***************************************
 * \fn diag_computeEnergyBalance(const COMPLEX *h)
 *
 * computes field energy \f$k_{\perp}\f$ spectra
 * \f$ W(k^{shell}_i) = \frac{1}{N}\sum_{k^{shell}_{i-1}<|\mathbf{k}_{\perp}|<k^{shell}_{i}}\sum_{k_z,l,m,s} g \bar{h}\f$
 * where \f$N\f$ is a number of wave vectors between shells \f$k^{shell}_{i-1}\f$ and \f$k^{shell}_{i}\f$
 ***************************************/
void diag_computeEnergyBalance(const COMPLEX *h, const COMPLEX *g) {
    diag_injected = 0;
    diag_dissipated = 0;
    diag_dissipated_kPerp = 0;
    diag_dissipated_kZ = 0;
    diag_dissipated_m = 0;
    double dissipated_k = 0;
    double dissipated_local_k;
    int forced_proc = mpi_whereIsM[2];
    size_t ind6D;
    size_t ind2D;
    // amount of energy diag_injected
    if (mpi_my_col_rank == forced_proc){
        for(size_t i = 0; i < equation_forceKn; i++){
            for(size_t is = 0; is < array_local_size.ns; is++){
                diag_injected += equation_forcingMM[i] * equation_forcingCoef;
            }
        }
    }
    if(mpi_my_rank == 0){
        MPI_Reduce(MPI_IN_PLACE, &diag_injected, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    else{
        MPI_Reduce(&diag_injected, &diag_injected, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    // amount of energy diag_dissipated
    for (size_t ix = 0; ix < array_local_size.nkx; ix++){
        for (size_t iy = 0; iy < array_local_size.nky; iy++){
            for (size_t iz = 0; iz < array_local_size.nkz; iz++){
                for (size_t im = 0; im < array_local_size.nm; im++){
                    for (size_t il = 0; il < array_local_size.nl; il++){
                        for (size_t is = 0; is < array_local_size.ns; is++){
                            ind2D = ix * array_local_size.nky + iy;
                            ind6D = get_flat_c(is,il,im,ix,iy,iz);
                            dissipated_local_k = creal(equation_getLocalDissipation(g[ind6D], space_kPerp2[ind2D], 0.0, 0.0)
                                    * conj(h[ind6D]) * diag_MM[iz]);
                            dissipated_k += dissipated_local_k;
                            diag_dissipated_kPerp += dissipated_local_k;
                            diag_dissipated_kZ += creal(equation_getLocalDissipation(g[ind6D], 0.0, space_kz[iz], 0.0)
                                                        * conj(h[ind6D]) * diag_MM[iz]);
                            diag_dissipated_m += creal(equation_getLocalDissipation(g[ind6D], 0.0, 0.0, global_nm_index[im])
                                                        * conj(h[ind6D]) * diag_MM[iz]);
                            diag_dissipated += creal(equation_getLocalDissipation(g[ind6D], space_kPerp2[ind2D], space_kz[iz], global_nm_index[im])
                                    * conj(h[ind6D]) * diag_MM[iz]);
                        }
                    }
                }
            }
        }
    }
    if(mpi_my_rank == 0){
        MPI_Reduce(MPI_IN_PLACE, &diag_dissipated, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        //computing eta kmax
        MPI_Reduce(MPI_IN_PLACE,&dissipated_k,1,MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        diag_etakmax = pow(var_var.mu_k * var_var.mu_k * var_var.mu_k / dissipated_k, 0.25) * space_kPerpMax;
    }
    else{
        MPI_Reduce(&diag_dissipated, &diag_dissipated, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&dissipated_k,&dissipated_k,1,MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
}

/***************************************
 * \fn diag_computeEnergy(const COMPLEX *h)
 *
 * computes field energy \f$k_{\perp}\f$ spectra
 * \f$ W(k^{shell}_i) = \frac{1}{N}\sum_{k^{shell}_{i-1}<|\mathbf{k}_{\perp}|<k^{shell}_{i}}\sum_{k_z,l,m,s} g \bar{h}\f$
 * where \f$N\f$ is a number of wave vectors between shells \f$k^{shell}_{i-1}\f$ and \f$k^{shell}_{i}\f$
 ***************************************/
void diag_computeEnergy(const COMPLEX *h){
    diag_energyH = 0;
    diag_energyPhi = 0;
    diag_energyBperp = 0;
    diag_energyBpar = 0;
    diag_energyTotal = 0;
    size_t ind6D;
    size_t ind3D;
    size_t ind2D;
    switch(systemType){
        case ELECTROSTATIC:
            for (size_t ix = 0; ix < array_local_size.nkx; ix++){
                for (size_t iy = 0; iy < array_local_size.nky; iy++){
                    for (size_t iz = 0; iz < array_local_size.nkz; iz++){
                        ind3D = get_flatIndexComplex3D(ix,iy,iz);
                        for (size_t is = 0; is < array_local_size.ns; is++){
                            diag_energyPhi += 0.5 * var_var.q[is] * var_var.q[is] * var_var.n[is] / var_var.T[is]
                                              * cabs(fields_fields.phi[ind3D]) * cabs(fields_fields.phi[ind3D]) * diag_MM[iz];
                        }
                        for (size_t im = 0; im < array_local_size.nm; im++){
                            for (size_t il = 0; il < array_local_size.nl; il++){
                                for (size_t is = 0; is < array_local_size.ns; is++){
                                    ind6D = get_flat_c(is,il,im,ix,iy,iz);
                                    diag_energyH += 0.5 * var_var.T[is] * var_var.n[is] * cabs(h[ind6D]) * cabs(h[ind6D]) * diag_MM[iz] ;

                                }
                            }
                        }
                    }
                }
            }
            MPI_Allreduce(MPI_IN_PLACE, &diag_energyH,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &diag_energyPhi,1,MPI_DOUBLE,MPI_SUM,mpi_row_comm);
            diag_energyTotal = diag_energyH - diag_energyPhi;
            break;
        case ELECTROMAGNETIC:
            for (size_t ix = 0; ix < array_local_size.nkx; ix++){
                for (size_t iy = 0; iy < array_local_size.nky; iy++){
                    ind2D = ix * array_local_size.nky  + iy;
                    for (size_t iz = 0; iz < array_local_size.nkz; iz++){
                        ind3D = get_flatIndexComplex3D(ix,iy,iz);
                        diag_energyBpar += cabs(fields_fields.B[ind3D]) * cabs(fields_fields.B[ind3D]) * diag_MM[iz];// / 8. / M_PI
                        diag_energyBperp += space_kPerp2[ind2D] * cabs(fields_fields.A[ind3D]) * cabs(fields_fields.A[ind3D])* diag_MM[iz];// / 8. / M_PI
                        for (size_t is = 0; is < array_local_size.ns; is++){
                            diag_energyPhi += 0.5 * var_var.q[is] * var_var.q[is] * var_var.n[is] / var_var.T[is]
                                              * cabs(fields_fields.phi[ind3D]) * cabs(fields_fields.phi[ind3D]) * diag_MM[iz];
                        }
                        for (size_t im = 0; im < array_local_size.nm; im++){
                            for (size_t il = 0; il < array_local_size.nl; il++){
                                for (size_t is = 0; is < array_local_size.ns; is++){
                                    ind6D = get_flat_c(is,il,im,ix,iy,iz);
                                    diag_energyH += 0.5 * var_var.T[is] * var_var.n[is] * cabs(h[ind6D]) * cabs(h[ind6D]) * diag_MM[iz];
                                }
                            }
                        }
                    }
                }
            }
            MPI_Allreduce(MPI_IN_PLACE, &diag_energyH,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
            MPI_Allreduce(MPI_IN_PLACE, &diag_energyPhi,1,MPI_DOUBLE,MPI_SUM,mpi_row_comm);
            MPI_Allreduce(MPI_IN_PLACE, &diag_energyBpar,1,MPI_DOUBLE,MPI_SUM,mpi_row_comm);
            MPI_Allreduce(MPI_IN_PLACE, &diag_energyBperp,1,MPI_DOUBLE,MPI_SUM,mpi_row_comm);
            diag_energyTotal = diag_energyH - diag_energyPhi + diag_energyBpar + diag_energyBperp;
            break;
        default:
            exit(1);
    }
}
/***************************************
 * \fn diag_computeEnergy(const COMPLEX *h)
 *
 * computes field energy \f$k_{\perp}\f$ spectra
 * \f$ W(k^{shell}_i) = \frac{1}{N}\sum_{k^{shell}_{i-1}<|\mathbf{k}_{\perp}|<k^{shell}_{i}}\sum_{k_z,l,m,s} g \bar{h}\f$
 * where \f$N\f$ is a number of wave vectors between shells \f$k^{shell}_{i-1}\f$ and \f$k^{shell}_{i}\f$
 ***************************************/
void diag_print(const COMPLEX *h, const COMPLEX *g, int it) {
    if (it%parameters.iter_diagnostics == 0){
        diag_computeEnergyBalance(h, g);
        diag_computeEnergy(h);
        if(mpi_my_rank == 0){
            if (mpi_my_rank == 0) printf("============TIMESTEP = %d============\n",it);
            printf("W/W0 = %.16f\n", diag_energyTotal / diag_free_energy0);
            printf("old W/W0 = %.16f\n", diag_freeEnergy / diag_free_energy0);
            printf("total energy = %f\n", diag_energyTotal);
            printf("old total energy = %f\n", diag_freeEnergy);
            printf("diag_injected = %f\n", diag_injected);
            printf("dissipated / injected = %f\n", diag_dissipated / diag_injected);
            printf("dissipated_kP / dissipated = %f\n",diag_dissipated_kPerp / diag_dissipated);
            printf("dissipated_kZ / dissipated = %f\n", diag_dissipated_kZ / diag_dissipated);
            printf("dissipated_m / dissipated = %f\n", diag_dissipated_m / diag_dissipated);
            printf("h energy / total = %f\n", diag_energyH / diag_energyTotal);
            printf("phi energy / total = %f\n", diag_energyPhi / diag_energyTotal);
            printf("B_par energy / total = %f\n", diag_energyBpar / diag_energyTotal);
            printf("B_perp energy / total = %f\n", diag_energyBperp / diag_energyTotal);
            printf("eta * kmax = %f\n", diag_etakmax);
            printf("L_dis/L_min = %f\n", diag_LDis/space_LperpMin);
            printf("L_int/L_max = %f\n", diag_LInt/space_LperpMax);

        }
    }
}

/***************************************
 * \fn diag_filterK(const COMPLEX *in, COMPLEX *out, k_c)
 * \brief low-pass sharp filter in \f$ k_\perp \f$
 * \param in: input 6D array
 * \param out: output 6D array
 * \param k_c: cutoff wave number
 ***************************************/
void diag_filterK(const COMPLEX *in, COMPLEX *out, double k_c){
    COMPLEX *maskKx = calloc(array_local_size.nkx, sizeof(*maskKx));
    COMPLEX *maskKy = calloc(array_local_size.nky, sizeof(*maskKx));
    size_t ind2D;
    size_t ind6D;

    for (size_t ix = 0; ix < array_local_size.nkx; ix++){
        for(size_t iy = 0; iy < array_local_size.nky; iy++){
            ind2D = ix * array_local_size.nky + iy;
            if (space_kPerp[ind2D] <= k_c){
                maskKx[ix] = 1.0;
                maskKy[iy] = 1.0;
            }
        }
    }
    for (size_t ix = 0; ix < array_local_size.nkx; ix++){
        for(size_t iy = 0; iy < array_local_size.nky; iy++){
            for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                for(size_t im = 0; im < array_local_size.nm; im++){
                    for(size_t il = 0; il < array_local_size.nl; il++){
                        for(size_t is = 0; is < array_local_size.ns; is++){
                            ind6D = get_flat_c(is,il,im,ix,iy,iz);
                            out[ind6D] =  in[ind6D] * maskKx[ix] * maskKy[iy];
                        }
                    }
                }
            }
        }
    }
    free(maskKx);
    free(maskKy);
}

/***************************************
 * \fn diag_filterKHP(const COMPLEX *in, COMPLEX *out, k_c)
 * \brief high-pass sharp filter in \f$ k_\perp \f$
 * \param in: input 6D array
 * \param out: output 6D array
 * \param k_c: cutoff wave number
 ***************************************/
void diag_filterKHP(const COMPLEX *in, COMPLEX *out, double k_c){
    COMPLEX *maskKx = calloc(array_local_size.nkx, sizeof(*maskKx));
    COMPLEX *maskKy = calloc(array_local_size.nky, sizeof(*maskKx));
    size_t ind2D;
    size_t ind6D;

    for (size_t ix = 0; ix < array_local_size.nkx; ix++){
        for(size_t iy = 0; iy < array_local_size.nky; iy++){
            ind2D = ix * array_local_size.nky + iy;
            if (space_kPerp[ind2D] <= k_c){
                maskKx[ix] = 1.0;
                maskKy[iy] = 1.0;
            }
        }
    }
    for (size_t ix = 0; ix < array_local_size.nkx; ix++){
        for(size_t iy = 0; iy < array_local_size.nky; iy++){
            for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                for(size_t im = 0; im < array_local_size.nm; im++){
                    for(size_t il = 0; il < array_local_size.nl; il++){
                        for(size_t is = 0; is < array_local_size.ns; is++){
                            ind6D = get_flat_c(is,il,im,ix,iy,iz);
                            out[ind6D] =  in[ind6D] * (1.0 - maskKx[ix] * maskKy[iy]);
                        }
                    }
                }
            }
        }
    }
    free(maskKx);
    free(maskKy);
}

/***************************************
 * \fn diag_computeNonlinearFlux(const COMPLEX *in)
 * \brief computes nonlinear free energy flux
 * \param in: input 6D array
 ***************************************/
void diag_computeNonlinearFlux(const COMPLEX *in){
    if (mpi_my_rank == IO_RANK) printf("computing nonlinear flux\n");
    size_t ind6D;
    COMPLEX *nonlinearTermK = calloc(array_local_size.total_comp, sizeof(*nonlinearTermK));
    equation_getNonlinearTerm(in,nonlinearTermK);
    fftw_copy_buffer_c(fftw_hBuf,nonlinearTermK);
    fftw_c2r();
    free(nonlinearTermK);
    double *nonlinearTermR = calloc(array_local_size.total_real, sizeof(*nonlinearTermR));
    fftw_copy_buffer_r(nonlinearTermR, fftw_hBuf);
    // forward cascade
    for (size_t ii = 0; ii < diag_numOfShells; ii++){
        diag_nonlinearFlux[ii] = 0;
        double kc = diag_shells[ii + 1]; //we don't want to compute flux through k = 0
        diag_filterK(in,fftw_hBuf,kc);
        fftw_c2r();
        double *h_r = fftw_hBuf;
        for(size_t iy = 0; iy < array_local_size.ny; iy++){
            for(size_t ix = 0; ix < array_local_size.nx; ix++){
                for(size_t iz = 0; iz < array_local_size.nz; iz++){
                    for(size_t im = 0; im < array_local_size.nm; im++){
                        for(size_t il = 0; il < array_local_size.nl; il++){
                            for(size_t is = 0; is < array_local_size.ns; is++){
                                ind6D = get_flat_r(is,il,im,ix,iy,iz);
                                diag_nonlinearFlux[ii] += h_r[ind6D] * nonlinearTermR[ind6D] * var_var.T[is] * var_var.n[is];
                            }
                        }
                    }
                }
            }
        }
    }

    //inverse cascade
    for (size_t ii = 0; ii < diag_numOfShells; ii++){
        diag_nonlinearFluxInverse[ii] = 0;
        double kc = diag_shells[ii + 1]; //we don't want to compute flux through k = 0
        diag_filterKHP(in,fftw_hBuf,kc);
        fftw_c2r();
        double *h_r = fftw_hBuf;
        for(size_t iy = 0; iy < array_local_size.ny; iy++){
            for(size_t ix = 0; ix < array_local_size.nx; ix++){
                for(size_t iz = 0; iz < array_local_size.nz; iz++){
                    for(size_t im = 0; im < array_local_size.nm; im++){
                        for(size_t il = 0; il < array_local_size.nl; il++){
                            for(size_t is = 0; is < array_local_size.ns; is++){
                                ind6D = get_flat_r(is,il,im,ix,iy,iz);
                                diag_nonlinearFluxInverse[ii] += h_r[ind6D] * nonlinearTermR[ind6D] * var_var.T[is] * var_var.n[is];
                            }
                        }
                    }
                }
            }
        }
    }



    free(nonlinearTermR);
    MPI_Allreduce(MPI_IN_PLACE, diag_nonlinearFlux, diag_numOfShells, MPI_DOUBLE, MPI_SUM, mpi_row_comm);
    MPI_Allreduce(MPI_IN_PLACE, diag_nonlinearFlux, diag_numOfShells, MPI_DOUBLE, MPI_SUM, mpi_col_comm);
    MPI_Allreduce(MPI_IN_PLACE, diag_nonlinearFluxInverse, diag_numOfShells, MPI_DOUBLE, MPI_SUM, mpi_row_comm);
    MPI_Allreduce(MPI_IN_PLACE, diag_nonlinearFluxInverse, diag_numOfShells, MPI_DOUBLE, MPI_SUM, mpi_col_comm);
    if (mpi_my_rank == IO_RANK){ printf("nonlinear flux computed\n");}
}