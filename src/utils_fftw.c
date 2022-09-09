/**************************************
* @file utils_fftw.c
* \brief FFT module
*
* contains FFT related routines
***************************************/
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
//  fftw_copyFieldBuf_r
//  fftw_copyFieldBuf_c
//  fftw_normalise_data
// fftw_normalise_chi_r
// fftw_normalise_field_r
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
#define VERBOSE 0



ptrdiff_t size_c[3];                            // full size of complex array
ptrdiff_t size_r[3];                            // full size of real array
ptrdiff_t howmany;                              // how many 3D transforms of box of size (nkx,nky,nkz) to perform
ptrdiff_t howmany_chi;
ptrdiff_t local_size, local_nx,local_ny,local_y_start, local_x_start;  // local size, local start and local kx block size on a given processor
ptrdiff_t local_size_chi, local_nx_chi, local_x_start_chi, local_ny_chi, local_y_start_chi;  // local size, local start and local kx block size on a given processor
ptrdiff_t local_size_field, local_nx_field, local_x_start_field, local_ny_field, local_y_start_field;
fftw_plan plan_c2r, plan_r2c;                   // plans needed to perform complex to real and real to complex transforms
fftw_plan plan_transposeToXY, plan_transposeToYX;
fftw_plan plan_c2r_chi, plan_r2c_chi;                   // plans needed to perform complex to real and real to complex transforms
fftw_plan plan_transposeToXY_chi, plan_transposeToYX_chi;
fftw_plan plan_c2r_field, plan_r2c_field;
fftw_plan plan_transposeToXY_field, plan_transposeToYX_field;
COMPLEX* fftw_hBuf;                                   // complex data array buffer; used for in-place
COMPLEX* fftw_chiBuf;                             // complex data array buffer for chi transformation;
COMPLEX *fftw_field;                              //complex data to transform fields
double fftw_norm;                               //normalization coefficient for backward fft transform
int *global_nkx_index;                          // array which stores the global position of nkx on the processor. Needed for dealiasing, in order to find Nkx/3 and 2*Nkx/3 and put all modes between those to zeros.

/***************************************
 * \fn  fftw_init(MPI_Comm communicator)
 * \brief initializes fftw transform.
 ***************************************/
void fftw_init(MPI_Comm communicator){
    /*initializing fftw and size of the transformation*/
    fftw_mpi_init();
    size_c[0] = array_global_size.nkx;
    size_c[1] = array_global_size.nky;
    size_c[2] = array_global_size.nkz;
    size_r[0] = array_global_size.ny;
    size_r[1] = array_global_size.nx;
    size_r[2] = array_global_size.nz;
    //size_r[0] = array_global_size.nkx;
    //size_r[1] = array_global_size.nky;
    //size_r[2] = array_global_size.nz;

    /*creating plans for full 6D transforms*/
    fftw_norm = 1./(array_global_size.nkx * array_global_size.nky * array_global_size.nz);

    howmany = array_local_size.nm * array_local_size.nl * array_local_size.ns;
    local_size = fftw_mpi_local_size_many_transposed(FFTW_RANK,
                                          size_c,
                                          howmany,
                                          array_local_size.nkx,
                                          array_local_size.ny,
                                          communicator,
                                          &local_nx,
                                          &local_x_start,
                                          &local_ny,
                                          &local_y_start); // getting local size stored on each processor;
    if (VERBOSE) printf("[MPI process %d] local size is %td, howmany is %d\n", mpi_my_rank,local_size, howmany);

    global_nkx_index = malloc(array_local_size.nkx * sizeof(*global_nkx_index));
    for (size_t i = 0; i < array_local_size.nkx; i++){
        global_nkx_index[i] = array_global_size.nkx / mpi_dims[1] * mpi_my_row_rank + i;
    }

    fftw_hBuf = fftw_alloc_complex(local_size);
    unsigned int flags_c2r = FFTW_MPI_TRANSPOSED_IN|FFTW_ESTIMATE;
    unsigned int flags_r2c = FFTW_MPI_TRANSPOSED_OUT|FFTW_ESTIMATE;
    plan_c2r = fftw_mpi_plan_many_dft_c2r(FFTW_RANK,
                                          size_r,
                                          howmany,
                                          local_nx,
                                          local_ny,
                                          fftw_hBuf,
                                          fftw_hBuf,
                                          communicator,
                                          flags_c2r);
    plan_r2c = fftw_mpi_plan_many_dft_r2c(FFTW_RANK,
                                          size_r,
                                          howmany,
                                          local_ny,
                                          local_nx,
                                          fftw_hBuf,
                                          fftw_hBuf,
                                          communicator,
                                          flags_r2c);
    plan_transposeToXY = fftw_mpi_plan_many_transpose(size_r[0], size_r[1],
                                                      howmany*(size_r[2]+2),
                                                      local_ny, local_nx,
                                                      fftw_hBuf, fftw_hBuf,
                                                      communicator,
                                                      FFTW_ESTIMATE);

    plan_transposeToYX = fftw_mpi_plan_many_transpose(size_r[1], size_r[0],
                                                      howmany*(size_r[2]+2),
                                                      local_nx, local_ny,
                                                      fftw_hBuf, fftw_hBuf,
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
    local_size_chi = fftw_mpi_local_size_many_transposed(FFTW_RANK,
                                                         size_c,
                                                         howmany_chi,
                                                         array_local_size.nkx,
                                                         array_local_size.ny,
                                                         communicator,
                                                         &local_nx_chi,
                                                         &local_x_start_chi,
                                                         &local_ny_chi,
                                                         &local_y_start_chi); // getting local size stored on each processor;


    if (VERBOSE) printf("[MPI process %d] CHI TRANSFORM: local size is %td, howmany is %d\n", mpi_my_rank,local_size_chi, howmany_chi);
    fftw_chiBuf = fftw_alloc_complex(local_size_chi);
    for (size_t ii = 0; ii < local_size_chi; ii++) fftw_chiBuf[ii] = 0;
    plan_c2r_chi = fftw_mpi_plan_many_dft_c2r(FFTW_RANK,
                                              size_r,
                                              howmany_chi,
                                              local_nx_chi,
                                              local_ny_chi,
                                              fftw_chiBuf,
                                              fftw_chiBuf,
                                              communicator,
                                              flags_c2r);
    plan_r2c_chi = fftw_mpi_plan_many_dft_r2c(FFTW_RANK,
                                              size_r,
                                              howmany_chi,
                                              local_ny_chi,
                                              local_nx_chi,
                                              fftw_chiBuf,
                                              fftw_chiBuf,
                                              communicator,
                                              flags_r2c);

    plan_transposeToXY_chi = fftw_mpi_plan_many_transpose(size_r[0], size_r[1],
                                                      howmany_chi * (size_r[2] + 2),
                                                      local_ny, local_nx,
                                                      fftw_chiBuf, fftw_chiBuf,
                                                      communicator,
                                                      FFTW_ESTIMATE);

    plan_transposeToYX_chi = fftw_mpi_plan_many_transpose(size_r[1], size_r[0],
                                                      howmany_chi * (size_r[2] + 2),
                                                      local_nx, local_ny,
                                                      fftw_chiBuf, fftw_chiBuf,
                                                      communicator,
                                                      FFTW_ESTIMATE);

    /* preparing field transform */
    //getting local size
    local_size_field = fftw_mpi_local_size_many_transposed(FFTW_RANK,
                                                           size_c,
                                                           1,
                                                           array_local_size.nkx,
                                                           array_local_size.ny,
                                                           communicator,
                                                           &local_nx_field,
                                                           &local_x_start_field,
                                                           &local_ny_field,
                                                           &local_y_start_field);

    if (VERBOSE) printf("[MPI process %d] FIELD TRANSFORM: local size is %td, howmany is %d\n", mpi_my_rank,local_size_field, 1);
    fftw_field = fftw_alloc_complex(local_size_field);
    for (size_t ii = 0; ii < local_size_field; ii++) fftw_field[ii] = 0;
    //creating r2c and c2r plans
    plan_c2r_field =  fftw_mpi_plan_dft_c2r_3d(size_r[0],
                                               size_r[1],
                                               size_r[2],
                                               fftw_field,
                                               fftw_field,
                                               communicator,
                                               flags_c2r);

    plan_r2c_field = fftw_mpi_plan_dft_r2c_3d(size_r[0],
                                              size_r[1],
                                              size_r[2],
                                              fftw_field,
                                              fftw_field,
                                              communicator,
                                              flags_r2c);

    plan_transposeToXY_field = fftw_mpi_plan_many_transpose(size_r[0], size_r[1],
                                                            (size_r[2] + 2),
                                                          local_ny, local_nx,
                                                          fftw_field, fftw_field,
                                                          communicator,
                                                          FFTW_ESTIMATE);

    plan_transposeToYX_field = fftw_mpi_plan_many_transpose(size_r[0], size_r[1],
                                                            (size_r[2] + 2),
                                                            local_ny, local_nx,
                                                            fftw_field, fftw_field,
                                                            communicator,
                                                            FFTW_ESTIMATE);
}

/***************************************
 * \fn fftw_r2c()
 * \brief real to complex fft transform.
 *
 * Performs real to complex in-place fft transform of on array <tt>fftw_hBuf</tt>.
 * Used to transform 6D arrays (kx,ky,kz,m,l,s).
 ***************************************/
void fftw_r2c() {
    int start = MPI_Wtime();
    fftw_mpi_execute_dft_r2c(plan_r2c, fftw_hBuf, fftw_hBuf);
    if (VERBOSE) printf("[MPI process %d] r2c transform performed in t = %.2fs.!\n", mpi_my_rank,MPI_Wtime()-start);
};

/***************************************
 * \fn fftw_c2r()
 * \brief complex to real fft transform.
 *
 * Performs complex to real in-place fft transform on array <tt>fftw_hBuf</tt>.
 * Used to transform 6D arrays (x,y,z,m,l,s).
 ***************************************/
void fftw_c2r() {
    int start = MPI_Wtime();
    fftw_mpi_execute_dft_c2r(plan_c2r, fftw_hBuf, fftw_hBuf);
    fftw_normalise_data_r(fftw_hBuf);
    if (VERBOSE) printf("[MPI process %d] c2r transform performed in t = %.2fs.!\n", mpi_my_rank,MPI_Wtime()-start);
}

/***************************************
 * \fn void fftw_r2c_chi()
 * \brief real to complex transform of chi potentials
 *
 * Performs real to complex in-place fft transform on array <tt>fftw_chiBuf</tt>.
 * Used to transform 5D arrays (x,y,z,s,field).
 ***************************************/
void fftw_r2c_chi() {
    int start = MPI_Wtime();
    fftw_mpi_execute_dft_r2c(plan_r2c_chi, fftw_chiBuf, fftw_chiBuf);
    if (VERBOSE) printf("[MPI process %d] r2c transform performed in t = %.2fs.!\n", mpi_my_rank,MPI_Wtime()-start);
};

/***************************************
 * \fn void fftw_c2r_chi()
 * \brief complex to real transform of chi potentials
 *
 * Performs complex to real in-place fft transform on array <tt>fftw_chiBuf</tt>.
 * Used to transform 5D arrays (kx,ky,kz,s,field).
 ***************************************/
void fftw_c2r_chi() {
    int start = MPI_Wtime();
    fftw_mpi_execute_dft_c2r(plan_c2r_chi, fftw_chiBuf, fftw_chiBuf);
    fftw_normalise_chi_r(fftw_chiBuf);
    if (VERBOSE) printf("[MPI process %d] c2r transform performed in t = %.2fs.!\n", mpi_my_rank,MPI_Wtime()-start);
};

/***************************************
 * \fn void fftw_r2c_field()
 * \brief real to complex transform of field potentials
 *
 * Performs real to complex in-place fft transform on array <tt>fftw_field</tt>.
 * Used to transform 3D arrays (x,y,z).
 ***************************************/
void fftw_r2c_field() {
    int start = MPI_Wtime();
    fftw_mpi_execute_dft_r2c(plan_r2c_field, fftw_field, fftw_field);
    if (VERBOSE) printf("[MPI process %d] r2c transform performed in t = %.2fs.!\n", mpi_my_rank,MPI_Wtime()-start);
};


/***************************************
 * \fn void fftw_c2r_field()
 * \brief complex to real transform of field potentials
 *
 * Performs complex to real in-place fft transform on array <tt>fftw_field</tt>.
 * Used to transform 3D arrays (kx,ky,kz).
 ***************************************/
void fftw_c2r_field() {
    int start = MPI_Wtime();
    fftw_mpi_execute_dft_c2r(plan_c2r_field, fftw_field, fftw_field);
    fftw_normalise_field_r(fftw_field);
    if (VERBOSE) printf("[MPI process %d] c2r field transform performed in t = %.2fs.!\n", mpi_my_rank,MPI_Wtime()-start);
};

/***************************************
 * \fn void fftw_kill()
 * \brief kills fftw
 *
 * to be added
 ***************************************/
void fftw_kill(){
    fftw_destroy_plan(plan_c2r);
    fftw_destroy_plan(plan_r2c);
    fftw_destroy_plan(plan_c2r_chi);
    fftw_destroy_plan(plan_r2c_chi);
    fftw_destroy_plan(plan_c2r_field);
    fftw_destroy_plan(plan_r2c_field);
    free(fftw_hBuf);
    free(fftw_chiBuf);
    free(fftw_field);
    free(global_nkx_index);
    fftw_mpi_cleanup();
}

/***************************************
 * \fn fftw_copy_buffer_r(double *to, double *from)
 * \brief copy 6D real array
 * \param to: where to copy array
 * \param from: array which will be copied
 *
 * copies real data <tt>from</tt> array to array <tt>to</tt>
 ***************************************/
void fftw_copy_buffer_r(double *to, double *from){
    for(size_t iky = 0; iky < array_local_size.ny; iky++){
        for(size_t ikx = 0; ikx < array_local_size.nx; ikx++){
            for(size_t iz = 0; iz < array_local_size.nz+2; iz++){
                for(size_t im = 0; im < array_local_size.nm; im++){
                    for(size_t il = 0; il < array_local_size.nl; il++){
                        for(size_t is = 0; is <array_local_size.ns; is++){
                                to[get_flat_r(is, il, im, ikx, iky, iz)] = from[get_flat_r(is, il, im, ikx, iky, iz)];
                        }
                    }
                }
            }
        }
    }
}

/***************************************
 * \fn fftw_copy_buffer_c(COMPLEX *to, COMPLEX *from)
 * \brief copy 6D complex array
 * \param to: where to copy array
 * \param from: array which will be copied
 *
 * copies complex data <tt>from</tt> array to array <tt>to</tt>
 ***************************************/
void fftw_copy_buffer_c(COMPLEX *to, COMPLEX *from){
    for(size_t ikx = 0; ikx < array_local_size.nkx; ikx++){
        for(size_t iky = 0; iky < array_local_size.nky; iky++){
            for(size_t iz = 0; iz < array_local_size.nkz; iz++){
                for(size_t im = 0; im < array_local_size.nm; im++){
                    for(size_t il = 0; il < array_local_size.nl; il++){
                        for(size_t is = 0; is <array_local_size.ns; is++){
                            to[get_flat_c(is, il, im, ikx, iky, iz)] = from[get_flat_c(is, il, im, ikx, iky, iz)];
                        }
                    }
                }
            }
        }
    }
}

/***************************************
 * \fn fftw_copyChiBuf_r(double *ar1, double *ar2)
 * \brief copy 5D real array
 * \param ar1: destination
 * \param ar2: source
 *
 * copies real \f$ \chi \f$ array from <tt>ar1</tt> to <tt>ar2</tt>.
 ***************************************/
void fftw_copyChiBuf_r(double *ar1, double *ar2){
    size_t ind;
    switch (systemType) {
        case ELECTROSTATIC:
            for(size_t iy = 0; iy < array_local_size.ny; iy++){
                for(size_t ix = 0; ix < array_local_size.nx; ix++){
                    for(size_t iz = 0; iz < array_local_size.nz+2; iz++){
                        for(size_t is = 0; is <array_local_size.ns; is++){
                            ind = getIndChiBufEL_r(ix,iy,iz,is);
                            ar1[ind] = ar2[ind];
                        }
                    }
                }
            }
            break;
        case ELECTROMAGNETIC:
            for(size_t iy = 0; iy < array_local_size.ny; iy++){
                for(size_t ix = 0; ix < array_local_size.nx; ix++){
                    for(size_t iz = 0; iz < array_local_size.nz+2; iz++){
                        for(size_t is = 0; is <array_local_size.ns; is++){
                            for(size_t ifield = 0; ifield < CHI_EM; ifield++){
                                ind = getIndChiBufEM_r(ix,iy,iz,is,ifield);
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
 * \fn fftw_copyChiBuf_c(COMPLEX *ar1, COMPLEX *ar2)
 * \brief copy 5D complex array.
 * \param ar1: destination
 * \param ar2: source
 *
 * copies complex \f$ \chi \f$ array from <tt>ar1</tt> to <tt>ar2</tt>.
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
 * \fn fftw_copyFieldBuf_r(double *to, double *from)
 * \brief copy 3D real data array.
 * \param to:
 * \param from:
 *
 * copies 3D data array <tt>from</tt> to <tt></tt>
 ***************************************/
void fftw_copyFieldBuf_r(double *to, double *from){
    for (size_t ii = 0; ii < array_local_size.nx * array_local_size.ny * (array_local_size.nz + 2); ii++){
        to[ii] = from[ii];
    }
}

/***************************************
 * \fn fftw_copyFieldBuf_c(COMPLEX *to, COMPLEX *from)
 * \brief copy 3D complex data array.
 * \param to:
 * \param from:
 *
 * copies 3D complex data array <tt>from</tt> to <tt></tt>
 ***************************************/
void fftw_copyFieldBuf_c(COMPLEX *to, COMPLEX *from){
    for (size_t ii = 0; ii < array_local_size.nkx * array_local_size.nky * array_local_size.nkz; ii++){
        to[ii] = from[ii];
    }
}


/***************************************
 * cosinus(double f,int ix)
 ***************************************/
double cosinus(double f,int ix){

    return cos((double)2. * M_PI * f / array_global_size.nkx * (local_x_start + ix));
}

/***************************************
 *  fftw_test_fill(double *ar,double f)
 ***************************************/
void fftw_test_fill(double *ar,double f){
    for(int i = 0; i < local_nx; i++){
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
void fftw_normalise_data(COMPLEX *data) {
    for(size_t i = 0; i < array_local_size.total_comp; i++) {
        data[i] *= fftw_norm;
    }
}

/***************************************
 * \fn void fftw_normalise_data_r(double *data)
 * \brief normalise data.
 * \param data: 6D data array
 *
 * normalises <tt>data</tt> by #fftw_norm.
 ***************************************/
void fftw_normalise_data_r(double *data) {
    for(size_t i = 0; i < array_local_size.total_real; i++) {
        data[i] *= fftw_norm;
    }
}

/***************************************
 * \fn void fftw_normalise_chi_r(double *data)
 * \brief notmalase chi data
 * \param data: 5D real data array
 *
 * normalises <tt>data</tt> by #fftw_norm.
 ***************************************/
void fftw_normalise_chi_r(double *data) {
    switch(systemType){
        case ELECTROSTATIC:
            for(size_t i = 0; i < array_local_size.nx * array_local_size.ny * (array_local_size.nz + 2) * array_local_size.ns; i++) {
                data[i] *= fftw_norm;
            }
            break;
        case ELECTROMAGNETIC:
            for(size_t i = 0; i < array_local_size.nx * array_local_size.ny * (array_local_size.nz + 2) * array_local_size.ns * 3; i++) {
                data[i] *= fftw_norm;
            }
            break;
        default: exit(1);
    }
}

/***************************************
 * \fn void fftw_normalise_field_r(double *data)
 * \brief normalise real 3D data
 * \param data: 3D real array
 *
 * normalises <tt>data</tt> by #fftw_norm.
 ***************************************/
void fftw_normalise_field_r(double *data) {
    for(size_t i = 0; i < array_local_size.nx * array_local_size.ny * (array_local_size.nz + 2); i++) {
        data[i] *= fftw_norm;
    }
}

/***************************************
 * \fn void fftw_transposeToXY()
 * \brief transposes 6D array
 ***************************************/
void fftw_transposeToXY(){
    fftw_execute(plan_transposeToXY);
}

/***************************************
 * \fn void fftw_transposeToXY()
 * \brief transposes 6D array
 ***************************************/
void fftw_transposeToYX(){
    fftw_execute(plan_transposeToYX);
}

/***************************************
 * \fn void fftw_transposeToXY_chi()
 * \brief transposes chi array
 ***************************************/
void fftw_transposeToXY_chi(){
    fftw_execute(plan_transposeToXY_chi);
}

/***************************************
 * \fn void fftw_transposeToXY_chi()
 * \brief transposes chi array
 ***************************************/
void fftw_transposeToYX_chi(){
    fftw_execute(plan_transposeToYX_chi);
}

/***************************************
 * \fn void fftw_transposeToXY_field()
 * \brief transposes field array
 ***************************************/
void fftw_transposeToXY_field(){
    fftw_execute(plan_transposeToXY_field);
}

/***************************************
 * \fn void fftw_transposeToYX_field()
 * \brief transposes field array
 ***************************************/
void fftw_transposeToYX_field(){
    fftw_execute(plan_transposeToYX_field);
}