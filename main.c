#include "mpi_utils.h"
#include "fftw_utils.h"
#include "hdf_utils.h"
#include "parameters_io.h"
#include "diagnostics.h"

int main(int argc, char **argv) {
    mpi_init();
    char *filename;
    if (argc<2){
        if(mpi_my_rank == 0){
            printf("OI! PROVIDE A FILENAME LAD!\n");
        }


    }
    else{
        if (strcmp("-f", argv[1]) == 0){
            filename = argv[2];
            read_parameters(filename);
        }
    }

    mpi_get_local_array_size();
    int start = MPI_Wtime();
    fftw_init(mpi_row_comm);
    printf("[MPI process %d] initialized fftw in t = %.2fs.!\n", mpi_my_rank,MPI_Wtime()-start);
    /*
    COMPLEX *arr_c = malloc(array_local_size.nkx*array_local_size.nky*array_local_size.nkz*array_local_size.nm*array_local_size.nl*array_local_size.ns*sizeof(*arr_c));
    double *arr_r =malloc(array_local_size.nkx*array_local_size.nky*(array_local_size.nz+2)*array_local_size.nm*array_local_size.nl*array_local_size.ns*sizeof(*arr_r));
    double *arr_r_sq =malloc(array_local_size.nkx*array_local_size.nky*(array_local_size.nz+2)*array_local_size.nm*array_local_size.nl*array_local_size.ns*sizeof(*arr_r));
    for (size_t i = 0;i<array_local_size.nkx*array_local_size.nky*(array_local_size.nz+2)*array_local_size.nm*array_local_size.nl*array_local_size.ns;i++){
        arr_r[i] = 0.;
    }
    for (size_t i = 0; i<array_local_size.nkx*array_local_size.nky*array_local_size.nkz*array_local_size.nm*array_local_size.nl*array_local_size.ns; i++){
        arr_c[i] = 0.;
    }

    //fftw_test_fill(arr_r,1.);
    for(int i = 0; i<array_local_size.nkx;i++){
        for(int j = 0; j<array_local_size.nky;j++){
            for(int k = 0; k<array_local_size.nz;k++){
               // arr_r[get_flat_r(0,0,0,i,j,k)] = cosinus(4.,i);
                //r_d[get_flat_r(0,0,0,i,j,k)] = cosinus(1.,i);
            }
        }
    }
    diag_fill_rand(arr_c);
    //printf("[MPI process %d] (%d,%d), m_local = %f\n",mpi_my_rank,mpi_my_coords[0],mpi_my_coords[1],arr_r[10]);

   // fftw_test_fill(arr_r,2.);
     */
    mpi_init_m_exchange();
    COMPLEX* big_array = alloc_complex6D(array_local_size.nkx,array_local_size.nky,array_local_size.nkz,array_local_size.nm,array_local_size.nl,array_local_size.ns);
    COMPLEX* minus_array = alloc_complex6D(array_local_size.nkx,array_local_size.nky,array_local_size.nkz,1,array_local_size.nl,array_local_size.ns);
    COMPLEX* plus_array = alloc_complex6D(array_local_size.nkx,array_local_size.nky,array_local_size.nkz,1,array_local_size.nl,array_local_size.ns);
    mpi_exchange_m_boundaries(big_array,plus_array,minus_array);
  //  arr_c[get_flat_c(0,0,0,1,0,0)] = 1.j;
    //fftw_r2c(arr_r,arr_c);
    COMPLEX *arr_c = malloc(array_local_size.nkx*array_local_size.nky*array_local_size.nkz*array_local_size.nm*array_local_size.nl*array_local_size.ns*sizeof(*arr_c));
    COMPLEX *arr_c1 = malloc(array_local_size.nkx*array_local_size.nky*array_local_size.nkz*array_local_size.nm*array_local_size.nl*array_local_size.ns*sizeof(*arr_c1));
    double *arr_r = malloc(array_local_size.nkx*array_local_size.nky*(array_local_size.nz+2)*array_local_size.nm*array_local_size.nl*array_local_size.ns*sizeof(*arr_r));
    double *arr_r1 = malloc(array_local_size.nkx*array_local_size.nky*(array_local_size.nz+2)*array_local_size.nm*array_local_size.nl*array_local_size.ns*sizeof(*arr_r1));
    double *arr_r_sq =malloc(array_local_size.nkx*array_local_size.nky*(array_local_size.nz+2)*array_local_size.nm*array_local_size.nl*array_local_size.ns*sizeof(*arr_r));

    for (size_t i = 0;i<array_local_size.nkx*array_local_size.nky*(array_local_size.nz+2)*array_local_size.nm*array_local_size.nl*array_local_size.ns;i++){
        arr_r[i] = 0.;
        arr_r1[i] = 0.;
        arr_r_sq[i] = 0.;
    }
    for (size_t i = 0; i<array_local_size.nkx*array_local_size.nky*array_local_size.nkz*array_local_size.nm*array_local_size.nl*array_local_size.ns; i++){
        arr_c[i] = 0.;
        arr_c1[i] = 0.;
    }

    // finding the wavevector closest to aliased modes
    for (size_t ikx = 0; ikx<array_local_size.nkx;ikx++){
        if(global_nkx_index[ikx] == 9){
            arr_c[get_flat_c(0,0,0,ikx,0,0)] = (1.j+1);
        }
        if(global_nkx_index[ikx] == 2){
            arr_c1[get_flat_c(0,0,0,ikx,0,0)] = (1.j+1);
        }
    }
    // computing ifft from the complex array
    fftw_c2r(arr_c,arr_r);
    fftw_c2r(arr_c1,arr_r1);
    //fftw_normalise_data(arr_r);
    //getting square from it
    diag_multiply_ar_r(arr_r1,arr_r,arr_r_sq);
    hdf_init();
    hdf_create_file_r("test.h5",arr_r);
    hdf_create_file_r("test1.h5",arr_r1);

    //now converting it back and dealiasing it
    fftw_r2c(arr_r_sq,arr_c);
    fftw_c2r(arr_c,arr_r_sq);
    fftw_normalise_data(arr_r_sq);
    hdf_create_file_r("test_non_deal.h5",arr_r_sq);
    dealiasing23(arr_c);
    fftw_c2r(arr_c,arr_r_sq);
    fftw_normalise_data(arr_r_sq);
    printf("[process id %d]zero mode equal to %f+%f i\n",mpi_my_rank, creal(arr_r_sq[get_flat_r(0,0,0,0,0,0)]),cimag(arr_r_sq[get_flat_r(0,0,0,0,0,0)]));
    hdf_create_file_r("test_deal.h5",arr_r_sq);

    fftw_normalise_data(arr_r);
    //hdf_create_file_r("test.h5",arr_r);
   // hdf_create_file_r("test_sq.h5",arr_r);
    //hdf_create_file_c("test1.h5",arr_c);
    //


    free(arr_c);
    free(arr_r);
    free(big_array);
    free(minus_array);
    free(plus_array);
    fftw_kill();
    mpi_kill();
    return 0;
}
