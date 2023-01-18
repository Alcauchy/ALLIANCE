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
    fields_getChi();
    distrib_getG(g, h);
    fields_sendF(g);
    fields_getFields(f00, f10, f01);
    hdf_saveData(h, 0);
    diag_compute(g, h, 0);

    diag_print(h, g, 0);
    //updating the time step size
    solver_updateDt(g, h, 0);

    // if not postprocessing, then we run main solver loop
    if (!parameters.postprocess) {
        //
        // main loop
        //
        char name[64];
        COMPLEX *freeEn = malloc(array_local_size.total_comp * sizeof(*freeEn));
        double start;
        double finish;
        double total;
        start = MPI_Wtime();
        for (int it = 1; it <= solver.Nt; it++) {
            //integrate over time
            solver_makeStep(&g, h, it);
            //updating the time step size
            solver_updateDt(g, h, it);
            //compute diagnostics
            diag_compute(g, h, it);
            //print data
            diag_print(h, g, it);
            //save data
            hdf_saveData(h, it);
            if (it % 10 == 0 && mpi_my_rank == 0) {
                printf("currently at it = %d\n", it);
            }
        }
        finish = MPI_Wtime();
        total = finish - start;
        MPI_Allreduce(MPI_IN_PLACE, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        if (mpi_my_rank == 0) {
            printf("TOTAL TIME OF SIMULATION = %f\n", total / mpi_size);
            printf("PER ONE TIMESTEP = %f\n", total / mpi_size / (solver.Nt - 1));
        }
        //
        // finalizing run
        //
        diag_compute(g, h, -1);
        hdf_saveData(h, -1);
        free(freeEn);
    }

    //
    // finalizing computations
    //
    free_wavespace();
    fftw_kill();
    mpi_kill();
    exit(0);
}
