#include "mpi_utils.h"
#include "fftw_utils.h"
#include "hdf_utils.h"
#include "parameters_io.h"
#include "diagnostics.h"
#include "space_config.h"
#include "init.h"
#include "fields.h"
#include "tests.h"

int main(int argc, char **argv) {
    //mpi_init();
    char *filename;
    if (argc<2){
        if(mpi_my_rank == 0){
            printf("OI! PROVIDE A FILENAME LAD!\n");
            exit(1);
        }


    }
    else{
        if (strcmp("-f", argv[1]) == 0){
            filename = argv[2];
        }
    }
    init_init(filename);
    //test_mainFunction();
    test_fieldComparison();
    //test_fieldComputation();
    free_wavespace();
    fftw_kill();
    mpi_kill();
    exit(0);
}
