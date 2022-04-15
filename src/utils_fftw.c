////////////////////////////////////////////////////////////////////////////////
// 08/12/2021 created by Gene Gorbunov
//                                   FFTW UTILITIES
//
//  fftw_init
//  fftw_kill
//  fftw_r2c
//  fftw_c2r
//  fftw_r2c_chi
//  fftw_c2r_chi
//  fftw_copy_buffer_r
//  fftw_copy_buffer_c
//  fftw_copyChiBuf_r
//  fftw_copyChiBuf_c
//  fftw_normalise_data
//  fftw_test_fill
//  dealiasing23
//  cosinus
//
// VERSION 1.0
////////////////////////////////////////////////////////////////////////////////

#include "utils_fftw.h"


#define FFTW_RANK 3
#define CHI_EL 1
#define CHI_EM 3



ptrdiff_t size_c[3];                            // full size of complex array
ptrdiff_t size_r[3];                            // full size of real array
ptrdiff_t howmany;                              // how many 3D transforms of box of size (nkx,nky,nkz) to perform
ptrdiff_t howmany_chi;
ptrdiff_t local_size, local_n0, local_0_start;  // local size, local start and local kx block size on a given processor
ptrdiff_t local_size_chi, local_n0_chi, local_0_start_chi;  // local size, local start and local kx block size on a given processor
fftw_plan plan_c2r, plan_r2c;                   // plans needed to perform complex to real and real to complex transforms
fftw_plan plan_c2r_chi, plan_r2c_chi;                   // plans needed to perform complex to real and real to complex transforms
COMPLEX* c_d;                                   // complex data array buffer; used for out-of-place
COMPLEX* c_chi_buf;                             // complex data array buffer for chi transformation;
double* r_d;                                    // real data array buffer; used for out-of-place transform;
double* r_chi_buf;                              // real data array buffer for chi transformations;
double fftw_norm;                               //normalization coefficient for backward fft transform
void (*fftw_dealiasing)(COMPLEX*) = NULL;
int *global_nkx_index;                          // array which stores the global position of nkx on the processor. Needed for dealiasing, in order to find Nkx/3 and 2*Nkx/3 and put all modes between those to zeros.

/***************************************
 * fftw_init(MPI_Comm communicator)
 ***************************************/
void fftw_init(MPI_Comm communicator){
    /*initializing fftw and size of the transformation*/
    fftw_mpi_init();
    size_c[0] = array_global_size.nkx;
    size_c[1] = array_global_size.nky;
    size_c[2] = array_global_size.nkz;
    size_r[0] = array_global_size.nkx;
    size_r[1] = array_global_size.nky;
    size_r[2] = array_global_size.nz;

    /*creating plans for full 6D transforms*/
    fftw_norm = 1./(array_global_size.nkx * array_global_size.nky * array_global_size.nz);

    howmany = array_local_size.nm * array_local_size.nl * array_local_size.ns;
    local_size = fftw_mpi_local_size_many(FFTW_RANK,
                                          size_c,
                                          howmany,
                                          array_local_size.nkx,
                                          communicator,
                                          &local_n0,
                                          &local_0_start); // getting local size stored on each processor;
    printf("[MPI process %d] local size is %td, howmany is %d\n", mpi_my_rank,local_size, howmany);

    global_nkx_index = malloc(array_local_size.nkx * sizeof(*global_nkx_index));
    for (size_t i = 0; i < array_local_size.nkx; i++){
        global_nkx_index[i] = array_global_size.nkx / mpi_dims[1] * mpi_my_row_rank + i;
        printf("[MPI process %d] my row rank = %d\t global_nkx_index[%d] = %d\n", mpi_my_rank,mpi_my_row_rank,i, global_nkx_index[i]);
    }

    c_d = fftw_alloc_complex(local_size);
    r_d = fftw_alloc_real(2 * local_size);

    plan_c2r = fftw_mpi_plan_many_dft_c2r(FFTW_RANK,
                                          size_r,
                                          howmany,
                                          local_n0,
                                          FFTW_MPI_DEFAULT_BLOCK,
                                          c_d,
                                          r_d,
                                          communicator,
                                          FFTW_ESTIMATE);
    plan_r2c = fftw_mpi_plan_many_dft_r2c(FFTW_RANK,
                                          size_r,
                                          howmany,
                                          local_n0,
                                          FFTW_MPI_DEFAULT_BLOCK,
                                          r_d,
                                          c_d,
                                          communicator,
                                          FFTW_ESTIMATE);

    /*preparing chi transform*/
    howmany_chi = array_local_size.ns;
    switch(systemType){
        case (ELECTROSTATIC):
            howmany_chi *= CHI_EL;
            break;
        case (ELECTROMAGNETIC):
            howmany_chi *= CHI_EM;
            break;
        default:
            printf("ERROR DEFINING NUMBER OF TRANSFORMS FOR CHI, ABORTING...\n");
            exit(1);
    }
    local_size_chi = fftw_mpi_local_size_many(FFTW_RANK,
                                          size_c,
                                          howmany_chi,
                                          array_local_size.nkx,
                                          communicator,
                                          &local_n0_chi,
                                          &local_0_start_chi); // getting local size stored on each processor;
    printf("[MPI process %d] CHI TRANSFORM: local size is %td, howmany is %d\n", mpi_my_rank,local_size_chi, howmany_chi);
    c_chi_buf = fftw_alloc_complex(local_size_chi);
    r_chi_buf = fftw_alloc_real(2 * local_size_chi);

    plan_c2r_chi = fftw_mpi_plan_many_dft_c2r(FFTW_RANK,
                                          size_r,
                                          howmany_chi,
                                          local_n0_chi,
                                          FFTW_MPI_DEFAULT_BLOCK,
                                          c_chi_buf,
                                          r_chi_buf,
                                          communicator,
                                          FFTW_ESTIMATE);
    plan_r2c_chi = fftw_mpi_plan_many_dft_r2c(FFTW_RANK,
                                          size_r,
                                          howmany_chi,
                                          local_n0_chi,
                                          FFTW_MPI_DEFAULT_BLOCK,
                                          r_chi_buf,
                                          c_chi_buf,
                                          communicator,
                                          FFTW_ESTIMATE);
}

/***************************************
 * fftw_r2c(double *data_r, COMPLEX *data_c)
 ***************************************/
void fftw_r2c(double *data_r, COMPLEX *data_c){
    int start = MPI_Wtime();
    fftw_copy_buffer_r(r_d,data_r);
    fftw_mpi_execute_dft_r2c(plan_r2c,r_d,c_d);
    fftw_copy_buffer_c(data_c, c_d);
    printf("[MPI process %d] r2c transform performed in t = %.2fs.!\n", mpi_my_rank,MPI_Wtime()-start);
};

/***************************************
 * fftw_c2r(COMPLEX *data_c, double *data_r)
 ***************************************/
void fftw_c2r(COMPLEX *data_c, double *data_r){
    int start = MPI_Wtime();
    if (mpi_my_coords[1] == 0){
    //    c_d[get_flat_c(0,0,0,1,0,0)] = array_global_size.nkx*array_global_size.nky*array_global_size.nkz;
    }

    fftw_copy_buffer_c(c_d, data_c);
    fftw_mpi_execute_dft_c2r(plan_c2r,c_d,r_d);
    fftw_copy_buffer_r(data_r,r_d);
   // fftw_normalise_data(data_r);
    printf("[MPI process %d] c2r transform performed in t = %.2fs.!\n", mpi_my_rank,MPI_Wtime()-start);
    if (mpi_my_rank == 0){
        for (int i = 0; i < array_local_size.nkx; i++ ){
            printf("[MPI process %d] data = %f\n", mpi_my_rank,r_d[get_flat_r(0,0,0,i,0,0)]);
        }

    }
}

/***************************************
 * fftw_r2c_chi(double *data_r, COMPLEX *data_c)
 ***************************************/
void fftw_r2c_chi(double *data_r, COMPLEX *data_c){
    int start = MPI_Wtime();
    fftw_copyChiBuf_r(r_d,data_r);
    fftw_mpi_execute_dft_r2c(plan_r2c_chi,r_chi_buf,c_chi_buf);
    fftw_copyChiBuf_c(data_c, c_d);
    printf("[MPI process %d] r2c transform performed in t = %.2fs.!\n", mpi_my_rank,MPI_Wtime()-start);
};

/***************************************
 * fftw_c2r_chi(COMPLEX *data_c, double *data_r)
 ***************************************/
void fftw_c2r_chi(COMPLEX *data_c, double *data_r){
    int start = MPI_Wtime();
    fftw_copyChiBuf_c(c_d, data_c);
    fftw_mpi_execute_dft_c2r(plan_c2r_chi,c_chi_buf,r_chi_buf);
    fftw_copyChiBuf_r(data_r,r_d);
    // fftw_normalise_data(data_r);
    printf("[MPI process %d] c2r transform performed in t = %.2fs.!\n", mpi_my_rank,MPI_Wtime()-start);
    if (mpi_my_rank == 0){
        for (int i = 0; i < array_local_size.nkx; i++ ){
            printf("[MPI process %d] data = %f\n", mpi_my_rank,r_d[get_flat_r(0,0,0,i,0,0)]);
        }

    }

};

/***************************************
 * fftw_kill()
 ***************************************/
void fftw_kill(){
    fftw_destroy_plan(plan_c2r);
    fftw_destroy_plan(plan_r2c);
    fftw_destroy_plan(plan_c2r_chi);
    fftw_destroy_plan(plan_r2c_chi);
    free(c_d);
    free(global_nkx_index);
    fftw_mpi_cleanup();
}

/***************************************
 * fftw_copy_buffer_r(double *ar1, double *ar2)
 ***************************************/
void fftw_copy_buffer_r(double *ar1, double *ar2){
    for(size_t ikx = 0; ikx < array_local_size.nkx; ikx++){
        for(size_t iky = 0; iky < array_local_size.nky; iky++){
            for(size_t iz = 0; iz < array_local_size.nz+2; iz++){
                for(size_t im = 0; im < array_local_size.nm; im++){
                    for(size_t il = 0; il < array_local_size.nl; il++){
                        for(size_t is = 0; is <array_local_size.ns; is++){
                            ar1[get_flat_r(is,il,im,ikx,iky,iz)] = ar2[get_flat_r(is,il,im,ikx,iky,iz)];
                        }
                    }
                }

            }
        }
    }
}

/***************************************
 * fftw_copy_buffer_c(COMPLEX *ar1, COMPLEX *ar2)
 ***************************************/
void fftw_copy_buffer_c(COMPLEX *ar1, COMPLEX *ar2){
    for(size_t ikx = 0; ikx < array_local_size.nkx; ikx++){
        for(size_t iky = 0; iky < array_local_size.nky; iky++){
            for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                for(size_t im = 0; im < array_local_size.nm; im++){
                    for(size_t il = 0; il < array_local_size.nl; il++){
                        for(size_t is = 0; is <array_local_size.ns; is++){
                            ar1[get_flat_c(is,il,im,ikx,iky,iz)] = ar2[get_flat_c(is,il,im,ikx,iky,iz)];
                        }
                    }
                }
            }
        }
    }
}

/***************************************
 * fftw_copyChiBuf_r(double *ar1, double *ar2)
 ***************************************/
void fftw_copyChiBuf_r(double *ar1, double *ar2){
    size_t ind;
    switch (systemType) {
        case ELECTROSTATIC:
            for(size_t ikx = 0; ikx < array_local_size.nkx; ikx++){
                for(size_t iky = 0; iky < array_local_size.nky; iky++){
                    for(size_t iz = 0; iz < array_local_size.nz+2; iz++){
                        for(size_t is = 0; is <array_local_size.ns; is++){
                            ind = getIndChiBufEL_r(ikx,iky,iz,is);
                            ar1[ind] = ar2[ind];
                        }
                    }
                }
            }
            break;
        case ELECTROMAGNETIC:
            for(size_t ikx = 0; ikx < array_local_size.nkx; ikx++){
                for(size_t iky = 0; iky < array_local_size.nky; iky++){
                    for(size_t iz = 0; iz < array_local_size.nz+2; iz++){
                        for(size_t is = 0; is <array_local_size.ns; is++){
                            for(size_t ifield = 0; ifield < CHI_EM; ifield++){
                                ind = getIndChiBufEM_r(ikx,iky,iz,is,ifield);
                                ar1[ind] = ar2[ind];
                            }
                        }
                    }
                }
            }
            break;
        default:
            printf("[MPI process %d] error copying real chi buffer\n", mpi_my_rank);
            exit(1);
    }
}

/***************************************
 * fftw_copyChiBuf_c(COMPLEX *ar1, COMPLEX *ar2)
 ***************************************/
void fftw_copyChiBuf_c(COMPLEX *ar1, COMPLEX *ar2){
    size_t ind;
    switch (systemType) {
        case ELECTROSTATIC:
            for(size_t ikx = 0; ikx < array_local_size.nkx; ikx++){
                for(size_t iky = 0; iky < array_local_size.nky; iky++){
                    for(size_t ikz = 0; ikz < array_local_size.nkz; ikz++){
                        for(size_t is = 0; is <array_local_size.ns; is++){
                            ind = getIndChiBufEL_c(ikx,iky,ikz,is);
                            ar1[ind] = ar2[ind];
                        }
                    }
                }
            }
            break;
        case ELECTROMAGNETIC:
            for(size_t ikx = 0; ikx < array_local_size.nkx; ikx++){
                for(size_t iky = 0; iky < array_local_size.nky; iky++){
                    for(size_t ikz = 0; ikz < array_local_size.nkz; ikz++){
                        for(size_t is = 0; is <array_local_size.ns; is++){
                            for(size_t ifield = 0; ifield < CHI_EM; ifield++){
                                ind = getIndChiBufEM_c(ikx,iky,ikz,is,ifield);
                                ar1[ind] = ar2[ind];
                            }
                        }
                    }
                }
            }
            break;
        default:
            printf("[MPI process %d] error copying complex chi buffer\n", mpi_my_rank);
            exit(1);
    }
}

/***************************************
 * cosinus(double f,int ix)
 ***************************************/
double cosinus(double f,int ix){

    return cos((double)2. * M_PI * f / array_global_size.nkx * (local_0_start + ix));
}

/***************************************
 * fftw_test_fill(double *ar,double f)
 ***************************************/
void fftw_test_fill(double *ar,double f){
    for(int i = 0; i<local_n0;i++){
        for(int j = 0; j<array_local_size.nky;j++){
            for(int k = 0; k<array_local_size.nz;k++){
                ar[get_flat_r(0,0,0,i,j,k)] = cosinus(f,i);
            }
        }
    }
}

/***************************************
 * fftw_normalise_data(double *data)
 ***************************************/
void fftw_normalise_data(double *data){
    for(size_t i = 0; i < array_local_size.total_real; i++) {
        data[i] *= fftw_norm;
    }
}

/***************************************
 * dealiasing23(COMPLEX *data_c)
 ***************************************/
void dealiasing23(COMPLEX *data_c){
    for(size_t ikx = 0; ikx < array_local_size.nkx; ikx++){
        for(size_t iky = 0; iky < array_local_size.nky; iky++){
            for(size_t ikz = 0; ikz < array_local_size.nkz; ikz++){
                for(size_t im = 0; im < array_local_size.nm; im++){
                    for(size_t il = 0; il < array_local_size.nl; il++){
                        for(size_t is = 0; is <array_local_size.ns; is++){
                            data_c[get_flat_c(is,il,im,ikx,iky,ikz)] =  (ikz>array_local_size.nkz*2/3)||
                                                                        ((iky>array_local_size.nky/3)&&(iky<array_local_size.nky*2/3))||
                                                                        ((global_nkx_index[ikx] > array_global_size.nkx / 3) && (global_nkx_index[ikx] < array_global_size.nkx * 2 / 3)) ? 0.j : data_c[get_flat_c(is, il, im, ikx, iky, ikz)];
                        }
                    }
                }

            }
        }
    }
}