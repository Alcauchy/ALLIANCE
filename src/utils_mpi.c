/**************************************
* @file utils_mpi.c
* \brief mpi module
*
* module to generate mpi topology and other mpi related routines
***************************************/
////////////////////////////////////////////////////////////////////////////////
// 06/12/2021 created by Gene Gorbunov
//                                   MPI UTILITIES
//
// mpi_init
// mpi_kill
// mpi_generateTopology
// mpi_createTopology
// mpi_findHermiteNeighbours
// mpi_splitInRows
// mpi_splitInCols
// mpi_initMExchange
// mpi_exchangeMBoundaries
// mpi_getLocalArraySize
// mpi_getLocalArrayOffsets
//
// VERSION 1.0
////////////////////////////////////////////////////////////////////////////////
#include "utils_mpi.h"
#define VERBOSE 0
#define SUBARRAY_COUNT 1
#define SUBARRAY_M_SIZE 1
#define SUBARRAY_DIMS 6
int mpi_my_rank;                        // rank of the process for 2D MPI communicator
int mpi_size;                           // size of the communicator
int mpi_my_row_rank;                    // rank of the process in a kx direction of parallelization
int mpi_my_col_rank;                    // rank of the process in Hermite direction of parallelization
int mpi_my_coords[2];                   // coordinates of the process in 2D topology
int mpi_dims[] = {0, 0};                // size of the dimensions. {0,0} means there is no limits in defining the size. 0th dimension is for Hermite parallelization, and 1st dimension is for kx.
int m_neighbour_ranks[2];               // ranks of the neighbours
int mpi_sub_buf_size;                   // buffer size needed to exchange m+1 and m-1 Hermite moments
int mpi_sub_buf_size_r;                 // buffer size needed to exchange m+1 and m-1 Hermite moments of real type
int *mpi_whereIsX;                       // arrays which returns the rank of a process which has a required kx
int *mpi_whereIsM;
size_t mpi_vectorSliceLength;
MPI_Comm mpi_cube_comm;                 // 2D topology communicator
MPI_Comm mpi_row_comm;                  // row communicator (kx direction)
MPI_Comm mpi_col_comm;                  // column communicator (Hermite direction)
MPI_Datatype mpi_subarray_type_plus;    // subarray type to perform the m+1 boundary exchange
MPI_Datatype mpi_subarray_type_minus;   // subarray type to perform the m-1 boundary exchange
MPI_Datatype mpi_subarray_type_plus_r;    // subarray type to perform the m+1 boundary exchange for real data
MPI_Datatype mpi_subarray_type_minus_r;   // subarray type to perform the m-1 boundary exchange for real data
MPI_Datatype mpi_vector_kxSlice;           // kx slice to enforce reality condition
enum DIRECTIONS {
    MINUS, PLUS
};           // minus and plus neighbours in m direction


/***************************
 *  mpi_init
 * *************************/
void mpi_init() {
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
}

/***************************
 *  mpi_generateTopology
 * *************************/
void mpi_generateTopology(){
    mpi_createTopology();
    mpi_getLocalArraySize();
    mpi_findHermiteNeighbours();
    mpi_splitInRows();
    mpi_splitInCols();
}

/***************************
 *  mpi_kill
 * *************************/
void mpi_kill() {
    MPI_Finalize();
    fftw_mpi_cleanup();
}

/***************************
 *  mpi_createTopology
 * *************************/
void mpi_createTopology() {
    int ndims = 2;
    int periods = {0, 0};
    int reorder = 1;
    MPI_Dims_create(mpi_size, ndims, &mpi_dims);
    MPI_Cart_create(MPI_COMM_WORLD, ndims, &mpi_dims, &periods, reorder, &mpi_cube_comm);
    MPI_Cart_coords(mpi_cube_comm, mpi_my_rank, ndims, mpi_my_coords);
    if (mpi_my_rank == 0) {
        printf("[MPI process %d] created %d x %d communicator\n", mpi_my_rank, mpi_dims[0], mpi_dims[1]);
    }
}

/***************************
 *  mpi_getLocalArraySize
 * *************************/
void mpi_getLocalArraySize() {
    array_local_size.nkx = array_global_size.nkx / mpi_dims[1];
    array_local_size.nky = array_global_size.nky;
    array_local_size.nkz = array_global_size.nkz;
    array_local_size.nm = array_global_size.nm / mpi_dims[0];
    array_local_size.nl = array_global_size.nl;
    array_local_size.ns = array_global_size.ns;
    array_local_size.nz = array_global_size.nz;
    if (mpi_my_coords[0] == mpi_dims[0] - 1) {
        array_local_size.nm = array_global_size.nm - array_global_size.nm / mpi_dims[0] * (mpi_my_coords[0]);
    }

    if (mpi_my_coords[1] == mpi_dims[1] - 1) {
        array_local_size.nkx = array_global_size.nkx - array_global_size.nkx / mpi_dims[1] * (mpi_my_coords[1]);
    }
    array_local_size.total_comp = array_local_size.nkx *
                                  array_local_size.nky *
                                  array_local_size.nkz *
                                  array_local_size.nm *
                                  array_local_size.nl *
                                  array_local_size.ns;

    array_local_size.total_real = array_local_size.nkx *
                                  array_local_size.nky *
                                  (array_local_size.nz + 2) *
                                  array_local_size.nm *
                                  array_local_size.nl *
                                  array_local_size.ns;

    mpi_getLocalArrayOffsets();

    printf("[MPI process %d] my coordinates = (%d,%d), local array size = (%d,%d,%d,%d,%d,%d), total size is %d\n",
           mpi_my_rank,
           mpi_my_coords[0],
           mpi_my_coords[1],
           array_local_size.nkx,
           array_local_size.nky,
           array_local_size.nkz,
           array_local_size.nm,
           array_local_size.nl,
           array_local_size.ns,
           array_local_size.total_comp);
}

/***************************
 *  mpi_getLocalArrayOffsets
 * *************************/
void mpi_getLocalArrayOffsets() {
    array_offset.kx = array_local_size.nky *
                      array_local_size.nkz *
                      array_local_size.nm *
                      array_local_size.nl *
                      array_local_size.ns;

    array_offset.ky = array_local_size.nkz *
                      array_local_size.nm *
                      array_local_size.nl *
                      array_local_size.ns;

    array_offset.kz = array_local_size.nm *
                      array_local_size.nl *
                      array_local_size.ns;

    array_offset.m = array_local_size.nl *
                     array_local_size.ns;

    array_offset.l = array_local_size.ns;

    array_offset.s = 1.0;

    array_offset.x = array_local_size.nky *
                     (array_local_size.nz + 2) *
                     array_local_size.nm *
                     array_local_size.nl *
                     array_local_size.ns;

    array_offset.y = (array_local_size.nz + 2) *
                     array_local_size.nm *
                     array_local_size.nl *
                     array_local_size.ns;

    array_offset.z = array_local_size.nm *
                     array_local_size.nl *
                     array_local_size.ns;

    array_offset3D.kx = array_local_size.nky *
                        array_local_size.nkz;

    array_offset3D.ky = array_local_size.nkz;

    array_offset3D.kz = 1;

};

/***************************
 *  mpi_findHermiteNeighbours
 * *************************/
void mpi_findHermiteNeighbours() {
    char *neighbour_names[2] = {"minus", "plus"};
    //performing shift in cartesian coordinates and finding plus and minus neighbours
    MPI_Cart_shift(mpi_cube_comm, 0, 1, &m_neighbour_ranks[MINUS], &m_neighbour_ranks[PLUS]);
    for (int i = 0; i < 2; i++) {
        if (m_neighbour_ranks[i] == MPI_PROC_NULL)
            printf("[MPI process %d] my cartesian coordinates are (%d,%d), I have no %s neighbour in m.\n", mpi_my_rank,
                   mpi_my_coords[0], mpi_my_coords[1], neighbour_names[i]);
        else
            printf("[MPI process %d] my cartesian coordinates are (%d,%d), I have a %s neighbour in m: process %d.\n",
                   mpi_my_rank, mpi_my_coords[0], mpi_my_coords[1], neighbour_names[i], m_neighbour_ranks[i]);
    }
}

/***************************
 *  mpi_splitInRows
 * *************************/
void mpi_splitInRows(){
    mpi_my_row_rank = mpi_my_coords[1];
    MPI_Comm_split(MPI_COMM_WORLD, mpi_my_coords[0], mpi_my_row_rank, &mpi_row_comm);
    printf("[MPI process %d] I am from row %d with row rank %d\n",mpi_my_rank,mpi_my_coords[0],mpi_my_row_rank);
    mpi_whereIsX = malloc(2 * array_global_size.nkx * sizeof(*mpi_whereIsX));
    mpi_whereIsM = malloc(2 * array_global_size.nm * sizeof(*mpi_whereIsM));
    for(int iproc = 0; iproc < mpi_dims[1]; iproc++){
        int ilocal = 0;
        if(iproc != mpi_dims[1] - 1){
            for (int ix = array_global_size.nkx/mpi_dims[1] * iproc; ix < array_global_size.nkx/mpi_dims[1] * (iproc + 1); ix++){
                mpi_whereIsX[ix * 2 + 0] = iproc;
                mpi_whereIsX[ix * 2 + 1] = ilocal;
                ilocal++;
            }
        }
        else{
            for (int ix = array_global_size.nkx/mpi_dims[1] * iproc; ix < array_global_size.nkx; ix++){
                mpi_whereIsX[ix * 2 + 0] = iproc;
                mpi_whereIsX[ix * 2 + 1] = ilocal;
                ilocal++;
            }
        }
    }
    for(int iproc = 0; iproc < mpi_dims[0]; iproc++){
        int ilocal = 0;
        if(iproc != mpi_dims[0] - 1){
            for (int im = array_global_size.nm/mpi_dims[0] * iproc; im < array_global_size.nm/mpi_dims[0] * (iproc + 1); im++){
                mpi_whereIsM[im * 2 + 0] = iproc;
                mpi_whereIsM[im * 2 + 1] = ilocal;
                ilocal++;
            }
        }
        else{
            for (int im = array_global_size.nm/mpi_dims[0] * iproc; im < array_global_size.nm; im++){
                mpi_whereIsM[im * 2 + 0] = iproc;
                mpi_whereIsM[im * 2 + 1] = ilocal;
                ilocal++;
            }
        }
    }
    if (mpi_my_rank == 0) for (size_t ii = 0; ii < array_global_size.nkx; ii++) printf("[MPI process %d] x = %zu is located at %d proc and local ix = %d\n",mpi_my_rank, ii, mpi_whereIsX[2 * ii], mpi_whereIsX[2 * ii + 1]);
    if (mpi_my_rank == 0) for (size_t ii = 0; ii < array_global_size.nm; ii++) printf("[MPI process %d] m = %zu is located at %d proc and local im = %d\n",mpi_my_rank, ii, mpi_whereIsM[2 * ii], mpi_whereIsM[2 * ii + 1]);
    // create mpi vector type to use for enforcing reality condition
    mpi_vectorSliceLength = array_local_size.nky*array_local_size.nl * array_local_size.nm * array_local_size.ns;
    MPI_Type_vector(array_local_size.nky,
                    array_local_size.nl * array_local_size.nm * array_local_size.ns,
                    array_local_size.nl * array_local_size.nm * array_local_size.ns * (array_local_size.nkz ) ,
                    MPI_C_DOUBLE_COMPLEX,
                    &mpi_vector_kxSlice);
    MPI_Type_commit(&mpi_vector_kxSlice);

}

/***************************
 *  mpi_splitInCols
 * *************************/
void mpi_splitInCols(){
    mpi_my_col_rank = mpi_my_coords[0];
    MPI_Comm_split(MPI_COMM_WORLD, mpi_my_coords[1], mpi_my_col_rank, &mpi_col_comm);
    printf("[MPI process %d] I am from column %d with column rank %d\n", mpi_my_rank, mpi_my_coords[1], mpi_my_col_rank);
};

/***************************
 *  mpi_initMExchange
 * *************************/
void mpi_initMExchange() {
    //Complex data init exchange
    int dimensions_full[6] = {array_local_size.nkx,
                              array_local_size.nky,
                              array_local_size.nkz,
                              array_local_size.nm,
                              array_local_size.nl,
                              array_local_size.ns};
    int dimensions_sub[6] = {array_local_size.nkx,
                             array_local_size.nky,
                             array_local_size.nkz,
                             SUBARRAY_M_SIZE,
                             array_local_size.nl,
                             array_local_size.ns};
    int start_coordinates_plus[6] = {0, 0, 0, array_local_size.nm - SUBARRAY_M_SIZE, 0, 0};
    int start_coordinates_minus[6] = {0, 0, 0, 0, 0, 0};
    mpi_sub_buf_size = array_local_size.nkx *
                       array_local_size.nky *
                       array_local_size.nkz *
                       array_local_size.nl *
                       array_local_size.ns;
    MPI_Type_create_subarray(SUBARRAY_DIMS,
                             dimensions_full,
                             dimensions_sub,
                             start_coordinates_plus,
                             MPI_ORDER_C,
                             MPI_C_DOUBLE_COMPLEX,
                             &mpi_subarray_type_plus);
    MPI_Type_create_subarray(SUBARRAY_DIMS,
                             dimensions_full,
                             dimensions_sub,
                             start_coordinates_minus,
                             MPI_ORDER_C,
                             MPI_C_DOUBLE_COMPLEX,
                             &mpi_subarray_type_minus);
    MPI_Type_commit(&mpi_subarray_type_plus);
    MPI_Type_commit(&mpi_subarray_type_minus);

    //Real data init exchange
    dimensions_full[2] = array_local_size.nz + 2;
    dimensions_sub[2] = array_local_size.nz + 2;
    mpi_sub_buf_size_r = array_local_size.nkx *
                         array_local_size.nky *
                        (array_local_size.nz + 2) *
                         array_local_size.nl *
                         array_local_size.ns;
    MPI_Type_create_subarray(SUBARRAY_DIMS,
                             dimensions_full,
                             dimensions_sub,
                             start_coordinates_plus,
                             MPI_ORDER_C,
                             MPI_DOUBLE,
                             &mpi_subarray_type_plus_r);
    MPI_Type_create_subarray(SUBARRAY_DIMS,
                             dimensions_full,
                             dimensions_sub,
                             start_coordinates_minus,
                             MPI_ORDER_C,
                             MPI_DOUBLE,
                             &mpi_subarray_type_minus_r);
    MPI_Type_commit(&mpi_subarray_type_plus_r);
    MPI_Type_commit(&mpi_subarray_type_minus_r);


}

/***************************
 *  mpi_exchangeMBoundaries
 * *************************/
void mpi_exchangeMBoundaries(COMPLEX *input_array, COMPLEX *plus_boundary, COMPLEX *minus_boundary) {
    double start;
    double end;
    start = MPI_Wtime();
    MPI_Send(input_array,
             SUBARRAY_COUNT,
             mpi_subarray_type_minus,
             m_neighbour_ranks[MINUS],
             0,
             mpi_cube_comm);
    MPI_Recv(plus_boundary,
             mpi_sub_buf_size,
             MPI_C_DOUBLE_COMPLEX,
             m_neighbour_ranks[PLUS],
             0,
             mpi_cube_comm,
             MPI_STATUS_IGNORE);
    if (VERBOSE) printf("[MPI process %d] received plus, t = %.2fs.! buf_size = %d, my_col_rank = %d\n", mpi_my_rank, MPI_Wtime() - start,
           mpi_sub_buf_size,mpi_my_col_rank);
    start = MPI_Wtime();
    MPI_Send(input_array, SUBARRAY_COUNT, mpi_subarray_type_plus, m_neighbour_ranks[PLUS], 1, mpi_cube_comm);
    MPI_Recv(minus_boundary, mpi_sub_buf_size, MPI_C_DOUBLE_COMPLEX, m_neighbour_ranks[MINUS], 1, mpi_cube_comm,
             MPI_STATUS_IGNORE);
    if (VERBOSE)printf("[MPI process %d] received minus, t = %.2fs.! buf_size = %d,my_col_rank = %d\n", mpi_my_rank,
           MPI_Wtime() - start,
           mpi_sub_buf_size,mpi_my_col_rank);
}

/***************************
 *  mpi_exchangeMBoundaries_r
 * *************************/
void mpi_exchangeMBoundaries_r(double *input_array, double *plus_boundary, double *minus_boundary) {
    double start;
    double end;
    start = MPI_Wtime();
    MPI_Send(input_array,
             SUBARRAY_COUNT,
             mpi_subarray_type_minus_r,
             m_neighbour_ranks[MINUS],
             0,
             mpi_cube_comm);
    MPI_Recv(plus_boundary,
             mpi_sub_buf_size_r,
             MPI_DOUBLE,
             m_neighbour_ranks[PLUS],
             0,
             mpi_cube_comm,
             MPI_STATUS_IGNORE);
    if (VERBOSE)  printf("[MPI process %d] received plus, t = %.2fs.! buf_size = %d\n", mpi_my_rank, MPI_Wtime() - start,
                         mpi_sub_buf_size);
    start = MPI_Wtime();
    MPI_Send(input_array,
             SUBARRAY_COUNT,
             mpi_subarray_type_plus_r,
             m_neighbour_ranks[PLUS],
             1,
             mpi_cube_comm);
    MPI_Recv(minus_boundary,
             mpi_sub_buf_size_r,
             MPI_DOUBLE,
             m_neighbour_ranks[MINUS],
             1,
             mpi_cube_comm,
             MPI_STATUS_IGNORE);
    if (VERBOSE) printf("[MPI process %d] received minus, t = %.2fs.! buf_size = %d\n", mpi_my_rank,
                        MPI_Wtime() - start,
                        mpi_sub_buf_size);
}

/***************************
 *  mpi_sendVector(COMPLEX *from_array, COMPLEX *to_buffer, int to_proc)
 * *************************/
void mpi_sendVector(COMPLEX *from_array, COMPLEX *to_buffer, int from_proc, int to_proc) {
    if(mpi_my_row_rank == from_proc) {
        MPI_Send(from_array, 1, mpi_vector_kxSlice, to_proc,1,mpi_row_comm);
        if (VERBOSE)printf("[MPI process %d] sent kx slice to process %d!\n",mpi_my_rank, to_proc);
        if (VERBOSE)printf("[MPI process %d] from =  %f!\n",mpi_my_rank, cabs(from_array[0]));
    }
    if(mpi_my_row_rank == to_proc) {
        MPI_Recv(to_buffer, mpi_vectorSliceLength, MPI_C_DOUBLE_COMPLEX, from_proc, 1, mpi_row_comm, MPI_STATUS_IGNORE);
        if (VERBOSE) printf("[MPI process %d] received kx slice from %d!\n", mpi_my_rank, from_proc);
    }
};