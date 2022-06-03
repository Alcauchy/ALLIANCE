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

#define RANK_IO 0

enum adiabatic kinetic;
enum electromagnetic systemType;
enum initial initialConditions;

/***************************************
 * init_start(char *filename)
 ***************************************/
void init_start(char *filename){
    srand(time(NULL));
    mpi_init();
    srand(mpi_my_rank);
    read_parameters(filename);
    init_initEnums();
    mpi_generateTopology();
    //mpi_getLocalArraySize();
    mpi_initMExchange();
    fftw_init(mpi_row_comm);
    space_init();
    diag_initSpec();
    var_init();
    fields_init();
    hdf_init();
    solver_init();
    init_printParameters();
};

/***************************************
 * init_printParameters()
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
 * init_initEnums()
 ***************************************/
void init_initEnums(){
   kinetic = parameters.adiabatic;
   systemType = parameters.electromagnetic;
   initialConditions = parameters.initial;
};

/***************************************
 * fill_rand(COMPLEX *ar1)
 ***************************************/
void fill_rand(COMPLEX *ar1) {
    for (size_t i = 0; i < array_local_size.total_comp; i++) {
        ar1[i] = 0;//cexp(2. * M_PI *1.j * (double) rand() / (double) (RAND_MAX)) * (array_global_size.nkx*array_global_size.nky*array_global_size.nz);
                //((0.5 - (double) rand() / (double) (RAND_MAX)) + (0.5 - (double) rand() / (double) (RAND_MAX)) * 1.j) * (array_global_size.nkx*array_global_size.nky*array_global_size.nz);
    }
    for (size_t ix = 0; ix < array_local_size.nkx; ix++){
        for (size_t iy = 0; iy < array_local_size.nky; iy++){
            for (size_t iz = 0; iz < array_local_size.nkz; iz++){
                for (size_t im = 0; im < array_local_size.nm; im++){
                    for (size_t il = 0; il < array_local_size.nl; il++){
                        for (size_t is = 0; is < array_local_size.ns; is++){
                            size_t ind6D = get_flat_c(is,il,im,ix,iy,iz);
                            size_t ind3D = ix * array_local_size.nky * array_local_size.nkz +
                                            iy * array_local_size.nkz +
                                            iz;
                            double theta = 2. * M_PI * (double) rand() / (double) (RAND_MAX);
                            ar1[ind6D] = cexp(1.j * theta) * (array_global_size.nkx*array_global_size.nky*array_global_size.nz);
                            if (space_kSq[ind3D] > 1e-10){
                                double amplitude = sqrt(init_energySpec(sqrt(space_kSq[ind3D]), 0, 0.01, 1.0) / 2.0/ M_PI);
                                ar1[ind6D] *=amplitude;
                            }
                            if(global_nkx_index[ix] == 0 && iy == 0 && iz == 0) ar1[ind6D] = 0;
                        }
                    }
                }
            }
        }
    }
}

/***************************************
 * fill_randM0(COMPLEX *ar1)
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
 * fill_randSingleKM(COMPLEX *ar1)
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
 * init_conditions(COMPLEX *data)
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
    distrib_enforceReality(data);
    distrib_setZeroNHalf(data);
    dealiasing23(data);
};

/***************************************
 * init_energySpec(COMPLEX *data)
 ***************************************/
double init_energySpec(double k, double m, double amp, double disp){
    return amp * k*k * exp( - 2.0 * (k * k / disp / disp  ));
};
