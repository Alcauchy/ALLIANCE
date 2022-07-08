#include "utils_mpi.h"
#include "utils_fftw.h"
#include "utils_hdf.h"
#include "parameters_io.h"
#include "diagnostics.h"
#include "space_config.h"
#include "init.h"
#include "fields.h"
#include "utils_tests.h"
#include <unistd.h>

int main(int argc, char **argv) {
    char *filename;
    if (argc < 3) {
        if (mpi_my_rank == 0) {
            printf("provide a filename\n");
            exit(1);
        }
    } else {
        if (strcmp("-f", argv[1]) == 0) {
            filename = argv[2];
        }
    }

    //
    //initializing run
    //
    init_start(filename);
    COMPLEX *h = malloc(array_local_size.total_comp * sizeof(*h));
    COMPLEX *g = malloc(array_local_size.total_comp * sizeof(*g));
    init_conditions(h);
    fields_sendF(h);
    fields_getFieldsFromH(f00, f10, f01);
    fields_getChi();
    distrib_getG(g, h);
    double free_energy0;

    //
    // main loop
    //
    for (int it = 0; it <= solver.Nt; it++) {
        hdf_saveData(h, it);
        if (parameters.save_diagnostics && it % parameters.iter_diagnostics == 0) {
            fields_sendF(g);
            fields_getFields(f00, f10, f01);
            fields_getChi();
            distrib_getH(h, g);
            diag_compute(g, h, it);
            if (it == 0) free_energy0 = diag_freeEnergy;
            if (mpi_my_rank == 0) printf("W = %.16f\n", diag_freeEnergy/free_energy0);
        }
        solver_makeStep(&g, h, it);
    }

    //
    // finalizing run
    //
    free_wavespace();
    fftw_kill();
    mpi_kill();
    exit(0);
}
