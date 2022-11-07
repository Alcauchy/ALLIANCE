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
    COMPLEX *h = calloc(array_local_size.total_comp, sizeof(*h));
    COMPLEX *g = calloc(array_local_size.total_comp, sizeof(*g));
    init_conditions(h);
    fields_sendF(h);
    fields_getFieldsFromH(f00, f10, f01);
    hdf_saveFieldPhi("field_h.h5");
    fields_getChi();
    distrib_getG(g, h);
    fields_sendF(g);
    fields_getFields(f00, f10, f01);
    diag_compute(g, h, 0);
    hdf_saveFieldPhi("field_g.h5");
    hdf_saveData(h, 0);
    diag_print(h,0);
    //updating the time step size
    solver_updateDt(g, h, 0);
    //
    // main loop
    //
    char name[64];
    COMPLEX *freeEn = malloc(array_local_size.total_comp * sizeof(*freeEn));
    for (int it = 1; it <= solver.Nt; it++) {
        //integrate over time
        solver_makeStep(&g, h, it);
        //updating the time step size
        solver_updateDt(g, h, it);
        //compute diagnostics
        diag_compute(g, h, it);
        //save data
        hdf_saveData(g, it);
        if (it%10000 == 0){
            char distrib_name[32];
            sprintf(distrib_name, "g_%d.h5", it);
            hdf_create_file_c(distrib_name,g);
        }
        //print data
        diag_print(h,it);
    }
    //
    // finalizing run
    //
    diag_compute(g, h, -1);
    hdf_saveData(h, -1);
    free_wavespace();
    fftw_kill();
    mpi_kill();
    exit(0);
}
