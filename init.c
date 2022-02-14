//
// Created by alcauchy on 03/02/2022.
//
#include "init.h"
#define RANK_IO 0

enum adiabatic kinetic;
enum electromagnetic systemType;
enum initial initialConditions;
void init_init(char *filename){
    mpi_init();
    read_parameters(filename);

    init_initEnums();
    mpi_generateTopology();
    mpi_get_local_array_size();
    fftw_init(mpi_row_comm);
    space_init();
    var_init();
    fields_init();
    hdf_init();
    solver_init();
    init_printParameters();
};

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

void init_initEnums(){
   kinetic = parameters.adiabatic;
   systemType = parameters.electromagnetic;
   initialConditions = parameters.initial;
};

void init_conditions(COMPLEX *data){
    switch(initialConditions)
    {
        case RANDOM:
            fill_rand(data);
            break;

        case FROMFILE:
            break;

        default:
            printf("[MPI process %d] error with initial conditions! Aborting...",mpi_my_rank);
            exit(1);
    }
};