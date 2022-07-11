/**************************************
* @file init.c
*   \brief initialization module for alliance.
*
*   all the inititalization routines are here.
***************************************/
////////////////////////////////////////////////////////////////////////////////
// 03/02/2022 created by Gene Gorbunov
//                                   INITIALIZATION
//
// init_start
// init_printParameters
// init_initEnums
// init_conditions
//
// VERSION 1.0
////////////////////////////////////////////////////////////////////////////////
#include "init.h"
#include "distrib.h"

/**************************************
* \def RANK_IO
*   defines rank of the processor
*   used to output information to console
***************************************/
#define RANK_IO 0

enum adiabatic kinetic;
enum electromagnetic systemType;
enum initial initialConditions;

/****************************************
* \fn void init_start(char *filename)
* \brief initialization of ALLIANCE
* \param filename: specifies parameter filename
*
* initializes all the modules required for ALLIANCE to work.
 ***************************************/
void init_start(char *filename){
    mpi_init();

    read_parameters(filename);
    init_computation();
    init_physicalSystem();
    hdf_createFiles();
    init_printParameters();
};

/****************************************
* \fn void init_physicalSystem()
* \brief initialization of physical system and parameters
 ***************************************/
void init_physicalSystem(){
    space_init();
    diag_initSpec();
    var_init();
    fields_init();
    solver_init();
};

/****************************************
* \fn void init_computation()
* \brief initialize hdf, fftw and mpi
 ***************************************/
void init_computation(){
    srand(0);
    init_initEnums();
    mpi_generateTopology();
    mpi_initMExchange();
    fftw_init(mpi_row_comm);
    hdf_init();
};

/***************************************
 * \fn void init_printParameters()
 * \brief parameter output
 *
 * prints parameters of the simulation
 ***************************************/
void init_printParameters(){
    if (mpi_my_rank == RANK_IO)
    {
        printf("created %d x %d communicator\n", mpi_dims[0], mpi_dims[1]);
        switch(kinetic)
        {
            case NONADIABATIC:
                printf("system is fully kinetic\n");
                break;
            case ADIABATIC:
                printf("system is adiabatic\n");
                break;
        }
        switch(systemType)
        {
            case ELECTROMAGNETIC:
                printf("system is electromagnetic\n");
                break;
            case ELECTROSTATIC:
                printf("system is electrostatic\n");
                break;
        }
    }

};

/***************************************
 * \fn void init_initEnums()
 * \brief enumerator initialization
 *
 * initializes enumerators, which are then used
 * to define if system is adiabatic or kinetic,
 * electromagnetic or electrostatic,
 * and type of initial conditions
 ***************************************/
void init_initEnums(){
   kinetic = parameters.adiabatic;
   systemType = parameters.electromagnetic;
   initialConditions = parameters.initial;
};

/***************************************
 * \fn void fill_rand(COMPLEX *data)
 * \brief fills the inital conditions randomly
 * \param data: complex 6d array to fill
 * initializes distribution with spectrum defined in
 * #init_energySpec
 * This function is supposed to be used in-module only
 * and should not be used elsewhere outside init.c file.
 ***************************************/
void fill_rand(COMPLEX *ar1) {
    for (size_t i = 0; i < array_local_size.total_comp; i++) {
        ar1[i] = 0;//cexp(2. * M_PI *1.j * (double) rand() / (double) (RAND_MAX)) * (array_global_size.nkx*array_global_size.nky*array_global_size.nz);
                //((0.5 - (double) rand() / (double) (RAND_MAX)) + (0.5 - (double) rand() / (double) (RAND_MAX)) * 1.j) * (array_global_size.nkx*array_global_size.nky*array_global_size.nz);
    }
    for (size_t ix = 0; ix < array_global_size.nkx; ix++){
        for (size_t iy = 0; iy < array_global_size.nky; iy++){
            for (size_t iz = 0; iz < array_global_size.nkz; iz++){
                for (size_t im = 0; im < array_global_size.nm; im++){
                    for (size_t il = 0; il < array_global_size.nl; il++){
                        for (size_t is = 0; is < array_global_size.ns; is++){
                            if (mpi_whereIsX[2*ix] == mpi_my_row_rank && mpi_whereIsM[2*im] == mpi_my_col_rank && (im == 0 || im == 1)){
                                size_t ix_local = mpi_whereIsX[2 * ix + 1];
                                size_t im_local = mpi_whereIsM[2 * im + 1];
                                size_t ind6D = get_flat_c(is,il,im_local,ix_local,iy,iz);
                                size_t ind3D = ix_local * array_local_size.nky * array_local_size.nkz +
                                               iy * array_local_size.nkz +
                                               iz;
                                double theta = 2. * M_PI * (double) rand() / (double) (RAND_MAX);
                                ar1[ind6D] = cexp(1.j * theta) * (array_global_size.nkx*array_global_size.nky*array_global_size.nz);
                                if (space_kSq[ind3D] > 1e-10){
                                    double amplitude = sqrt(init_energySpec(sqrt(space_kSq[ind3D]), 0, 1., .5) / 2.0/ M_PI);
                                    ar1[ind6D] *=amplitude;
                                }
                                if(global_nkx_index[ix_local] == 0 && iy == 0 && iz == 0) ar1[ind6D] = 0;
                            }
                            else{
                                rand();
                            }
                        }
                    }
                }
            }
        }
    }
}


/***************************************
 * \fn void fill_randM0(COMPLEX *data)
 * \brief fill zeroth Hermite moment with random values
 * \param data: complex 6D array to fill
 *
 * fills 0-th Hermite moment of a distribution function ar1 with random values
 * This function is supposed to be used in-module only
 * and should not be used elsewhere outside init.c file.
 ***************************************/
void fill_randM0(COMPLEX *ar1) {
    size_t ind6D;
    for (size_t i = 0; i < array_local_size.total_comp; i++) {
        ar1[i] = 0.;
    }
    for (size_t ix = 0; ix < array_local_size.nkx; ix++) {
        for (size_t iy = 0; iy < array_local_size.nky; iy++) {
            for (size_t iz = 0; iz < array_local_size.nkz; iz++) {
                for (size_t il = 0; il < array_local_size.nl; il++) {
                    for (size_t im = 0; im < array_local_size.nm; im++) {
                        for (size_t is = 0; is < array_local_size.ns; is++) {
                            if (space_globalMIndex[im] == 0 || space_globalMIndex[im] == 1) {
                                //printf("%d\n",im);
                                ind6D = get_flat_c(is, il, im, ix, iy, iz);
                                ar1[ind6D] = ((0.5 - (double) rand() / (double) (RAND_MAX)) +
                                             (0.5 - (double) rand() / (double) (RAND_MAX)) * 1.j) * (array_global_size.nkx*array_global_size.nky*array_global_size.nz);
                            }

                        }
                    }
                }
            }
        }
    }
}

/***************************************
 * \fn void fill_randSingleKM(COMPLEX *ar1)
 * \brief fill single chosen wavevector and Hermite moment
 * \param data: complex 6D array
 *
 * initializes single wavevector and Hermite moment
 * of a distribution function with random variable.
 * This function is only for in-module use
 * and should not be used elsewhere outside init.c file.
 ***************************************/
void fill_randSingleKM(COMPLEX *ar1) {
    size_t ind6D;
    for (size_t i = 0; i < array_local_size.total_comp; i++) {
        ar1[i] = 0.;
    }
    for (size_t ix = 0; ix < array_local_size.nkx; ix++) {
       // printf("[MPI process %d] ix = %zu, %d\n",mpi_my_rank, ix, global_nkx_index[ix]);
        for (size_t iy = 0; iy < array_local_size.nky; iy++) {
            for (size_t iz = 0; iz < array_local_size.nkz; iz++) {
                for (size_t im = 0; im < array_local_size.nm; im++) {
                    //printf("im = %d\n", space_globalMIndex[im]);
                    for (size_t il = 0; il < array_local_size.nl; il++) {
                        for (size_t is = 0; is < array_local_size.ns; is++) {
                            if (space_globalMIndex[im] == 1 || space_globalMIndex[im] == 0 && global_nkx_index[ix] == 0) {

                                //printf("%f\n", space_kz[2]);
                                ind6D = get_flat_c(is, il, im, ix, 0, 2);
                                ar1[ind6D] = ((0.5 - (double) rand() / (double) (RAND_MAX)) +
                                                   (0.5 - (double) rand() / (double) (RAND_MAX)) * 1.j)
                                                           * (array_global_size.nkx*array_global_size.nky*array_global_size.nz);
                            };
                        }
                    }
                }
            }
        }
    }
}

/***************************************
 * \fn void init_conditions(COMPLEX *data)
 * \brief distribution function initialization
 * \param data: complex 6D array
 *
 * initializes distribution function with chosen method
 * (see #fill_rand,
 * #fill_randM0,
 * #fill_randSingleKM)
 ***************************************/
void init_conditions(COMPLEX *data){
    switch(initialConditions)
    {
        case RANDOM:
            fill_rand(data);
            //fill_randM0(data);
            //fill_randSingleKM(data);
            break;

        case FROMFILE:
            hdf_readData(parameters.from_simulationName,data);
            break;

        default:
            printf("[MPI process %d] error with initial conditions! Aborting...",mpi_my_rank);
            exit(1);
    }
    dealiasing23(data);
    distrib_enforceReality(data);
    distrib_setZeroNHalf(data);
};

/***************************************
 * \fn double init_energySpec(double k, double m, double amp, double disp)
 * \brief returns energy spectrum
 * \param k: a wavenumber at which spectrum is computed
 * \param m: Hermite moment at which amplitude is computed
 * \param amp: amplitude of the spectrum
 * \param disp: dispersion of the spectrum
 *
 * computes spectrum of form
 * \f$ A \cdot k^2 exp(-2 k^2/\sigma^2) \f$, where
 * \f$\sigma = disp\f$, and
 * \f$A = amp \f$
 * This function is supposed to be used in-module only
 * and should not be used elsewhere outside init.c file.
 ***************************************/
double init_energySpec(double k, double m, double amp, double disp){
    return amp * k*k * exp( - 2.0 * (k * k / disp / disp  ));
};
