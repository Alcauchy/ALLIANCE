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
    if (argc < 2) {
        if (mpi_my_rank == 0) {
            printf("provide a filename\n");
            exit(1);
        }
    } else {
        if (strcmp("-f", argv[1]) == 0) {
            filename = argv[2];
        }
    }

    init_start(filename);
    hdf_createFiles();

    COMPLEX *h = malloc(array_local_size.total_comp * sizeof(*h));
    COMPLEX *g = malloc(array_local_size.total_comp * sizeof(*g));
    init_conditions(h);
    fields_sendG(h);
    fields_getFieldsFromH(g00, g10, g01);
    fields_getChi();
    distrib_getG(g, h);
    for (int it = 0; it < solver.Nt; it++) {
        solver_makeStep();
        //solver_updateDt();
        if (parameters.save_diagnostics && it % parameters.iter_diagnostics == 0) {
            diag_compute(g, h, it);
        }
        hdf_saveData(h, it);
    }
    free_wavespace();
    fftw_kill();
    mpi_kill();
    exit(0);
}
