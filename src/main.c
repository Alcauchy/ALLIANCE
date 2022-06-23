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
    test_transposedFFTW();
    /*COMPLEX *h = malloc(array_local_size.total_comp * sizeof(*h));
    COMPLEX *g = malloc(array_local_size.total_comp * sizeof(*g));
    init_conditions(h);
    fields_sendG(h);
    fields_getFieldsFromH(g00, g10, g01);
    fields_getChi();
    distrib_getG(g, h);
    double free_energy0;
    for (int it = 0; it <= solver.Nt; it++) {
        hdf_saveData(h, it);
        if(it%10 == 0) printf("it = %d\n", it);
        //solver_updateDt();
        if (parameters.save_diagnostics && it % parameters.iter_diagnostics == 0) {

            fields_sendG(g);
            fields_getFields(g00, g10, g01);
            fields_getChi();
            distrib_getH(h, g);
            diag_compute(g, h, it);
            //
            // saving fields and everything
            //
            if (parameters.save_diagnostics && it % 1000 == 0){
                char name[64];
                //save g and h
                sprintf(name, "%s%s%s%d%s", ".","/","g_",it,".h5");
                hdf_create_file_c(name,g);
                sprintf(name, "%s%s%s%d%s", ".","/","h_",it,".h5");
                hdf_create_file_c(name,h);

            }


            // save real g and h
            //fftw_copy_buffer_c(fftw_hBuf,g);
            //fftw_c2r();
            //sprintf(name, "%s%s%s%d%s", ".","/","gr_",it,".h5");
            //hdf_create_file_r(name,fftw_hBuf);

            //fftw_copy_buffer_c(fftw_hBuf,h);
            //fftw_c2r();
            //sprintf(name, "%s%s%s%d%s", ".","/","hr_",it,".h5");
            //hdf_create_file_r(name,fftw_hBuf);

            // save fields
            //fftw_copyFieldBuf_c(fftw_field,fields_fields.phi);
            //fftw_c2r_field();
            //sprintf(name, "%s%s%s%d%s", ".","/","phir_",it,".h5");
            //hdf_saveField_r(fftw_field,name);

            //fftw_copyFieldBuf_c(fftw_field,fields_fields.A);
            //fftw_c2r_field();
            //sprintf(name, "%s%s%s%d%s", ".","/","Ar_",it,".h5");
            //hdf_saveField_r(fftw_field,name);

            //fftw_copyFieldBuf_c(fftw_field,fields_fields.B);
            //fftw_c2r_field();
            //sprintf(name, "%s%s%s%d%s", ".","/","Br_",it,".h5");
            //hdf_saveField_r(fftw_field,name);*/

           // if (it == 0) free_energy0 = diag_freeEnergy;
           // if (mpi_my_rank == 0) printf("W = %.16f\n", diag_freeEnergy/free_energy0);
       // }
        //printf("1)1st = %p, 2nd = %p\n", g, rk4.g_buf);
       // solver_makeStep(&g, h);
        //printf("2)1st = %p, 2nd = %p\n", g, rk4.g_buf);
    //}

    //test_linearRHS();
    //test_inplaceFFTW_chi();

    //test_Poisson();
    //test_Poisson1();
    //test_xGrad();
    //test_enforceRealityConditions();
    //test_enforceZero();
    //test_fieldComparison();
    //test_fieldsFFT();
    //test_Poisson();
    //test_nonlinearTerm();
    //test_everything();
    //test_RHS();
    free_wavespace();
    fftw_kill();
    mpi_kill();
    exit(0);
}
