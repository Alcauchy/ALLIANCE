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
#define RANK_IO 0

enum adiabatic kinetic;
enum electromagnetic systemType;
enum initial initialConditions;

/***************************************
 * init_start(char *filename)
 ***************************************/
void init_start(char *filename){
    mpi_init();
    read_parameters(filename);

    init_initEnums();
    mpi_generateTopology();
    mpi_getLocalArraySize();
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
 * init_conditions(COMPLEX *data)
 ***************************************/
void init_conditions(COMPLEX *data){
    switch(initialConditions)
    {
        case RANDOM:
            //fill_rand(data);
            fill_randM0(data);
           // fill_randSingleKM(data);
            break;

        case FROMFILE:
            hdf_readData(parameters.from_simulationName,data);
            break;

        default:
            printf("[MPI process %d] error with initial conditions! Aborting...",mpi_my_rank);
            exit(1);
    }
};