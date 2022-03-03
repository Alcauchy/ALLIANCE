#include "utils_mpi.h"
#include "utils_fftw.h"
#include "utils_hdf.h"
#include "parameters_io.h"
#include "diagnostics.h"
#include "space_config.h"
#include "init.h"
#include "fields.h"
#include "utils_tests.h"

int main(int argc, char **argv) {
    char *filename;
    if (argc<2){
        if(mpi_my_rank == 0){
            printf("provide a filename\n");
            exit(1);
        }
    }
    else{
        if (strcmp("-f", argv[1]) == 0){
            filename = argv[2];
        }
    }

    init_init(filename);
    hdf_createFiles();
    test_mainFunction();
    //test_kSpecComputations();
    free_wavespace();
    fftw_kill();
    mpi_kill();
    exit(0);
}
