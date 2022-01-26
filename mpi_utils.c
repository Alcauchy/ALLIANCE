//
// Created by alcauchy on 06/12/2021.
//
#include "mpi_utils.h"

int mpi_my_rank;                        // rank of the process for 2D MPI communicator
int mpi_size;                           // size of the communicator
int mpi_my_row_rank;                    // rank of the process in a kx direction of parallelization
int mpi_my_coords[2];                   // coordinates of the process in 2D topology
int mpi_dims[] = {0, 0};                // size of the dimensions. {0,0} means there is no limits in defining the size.
enum DIRECTIONS {MINUS,PLUS};           // minus and plus neighbours in m direction
int m_neighbour_ranks[2];               // ranks of the neighbours
int mpi_sub_buf_size;                   // buffer size needed to exchange m+1 and m-1 Hermite moments
MPI_Comm mpi_cube_comm;                 // 2D topology communicator
MPI_Comm mpi_row_comm;                  // row communicator (kx direction)
MPI_Datatype mpi_subarray_type_plus;    // subarray type to perform the m+1 boundary exchange
MPI_Datatype mpi_subarray_type_minus;   // subarray type to perform the m-1 boundary exchange



void mpi_init(){
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    mpi_create_topology();
    mpi_find_hermite_neighbours();
    mpi_split_in_rows();

}


void mpi_kill(){
    MPI_Finalize();
    fftw_mpi_cleanup();
}


void mpi_create_topology(){
    int ndims = 2;
    int periods = {0,0};
    int reorder = 1;
    MPI_Dims_create(mpi_size, ndims,  &mpi_dims);
    MPI_Cart_create(MPI_COMM_WORLD, ndims, &mpi_dims, &periods, reorder, &mpi_cube_comm);
    MPI_Cart_coords(mpi_cube_comm, mpi_my_rank, ndims, mpi_my_coords);
    if (mpi_my_rank == 0){ printf("[MPI process %d] created %d x %d communicator\n",mpi_my_rank, mpi_dims[0],mpi_dims[1]); }
}


void mpi_get_local_array_size(){
    array_local_size.nkx = parameters.nkx / mpi_dims[1];
    array_local_size.nky = parameters.nky;
    array_local_size.nkz = parameters.nkz;
    array_local_size.nm = parameters.nm / mpi_dims[0];
    array_local_size.nl = parameters.nl;
    array_local_size.ns = parameters.ns;
    array_local_size.nz = parameters.nz;
    if(mpi_my_coords[0] == mpi_dims[0] - 1){
        array_local_size.nm = parameters.nm- parameters.nm / mpi_dims[0] * (mpi_my_coords[0]);
    }

    if(mpi_my_coords[1] == mpi_dims[1] - 1){
        array_local_size.nkx = parameters.nkx- parameters.nkx / mpi_dims[1] * (mpi_my_coords[1]);
    }
    int total_size  = array_local_size.nkx * array_local_size.nky * array_local_size.nkz * array_local_size.nm * array_local_size.nl * array_local_size.ns;

    printf("[MPI process %d] my coordinates = (%d,%d), local array size = (%d,%d,%d,%d,%d,%d), total size is %d\n",mpi_my_rank, mpi_my_coords[0], mpi_my_coords[1],
                                                                                                  array_local_size.nkx, array_local_size.nky, array_local_size.nkz,
                                                                                                  array_local_size.nm, array_local_size.nl, array_local_size.ns, total_size);
}

void mpi_find_hermite_neighbours(){
    char* neighbour_names[2] = { "minus", "plus" };
    //performing shift in cartesian coordinates and finding plus and minus neighbours
    MPI_Cart_shift(mpi_cube_comm, 0, 1, &m_neighbour_ranks[MINUS], &m_neighbour_ranks[PLUS]);
    for (int i = 0; i < 2; i++){
        if (m_neighbour_ranks[i] == MPI_PROC_NULL)
            printf("[MPI process %d] my cartesian coordinates are (%d,%d), I have no %s neighbour in m.\n", mpi_my_rank, mpi_my_coords[0], mpi_my_coords[1], neighbour_names[i]);
        else
            printf("[MPI process %d] my cartesian coordinates are (%d,%d), I have a %s neighbour in m: process %d.\n", mpi_my_rank, mpi_my_coords[0], mpi_my_coords[1], neighbour_names[i],m_neighbour_ranks[i]);
    }
}

void mpi_split_in_rows(){
    mpi_my_row_rank = mpi_my_coords[1];
    MPI_Comm_split(MPI_COMM_WORLD, mpi_my_coords[0], mpi_my_row_rank, &mpi_row_comm);
    printf("[MPI process %d] I am from row %d with row rank %d\n",mpi_my_rank,mpi_my_coords[0],mpi_my_row_rank);
}

void mpi_init_m_exchange() {
    int dimensions_full[6] = {array_local_size.nkx, array_local_size.nky, array_local_size.nkz, array_local_size.nm, array_local_size.nl, array_local_size.ns};
    int dimensions_sub[6] = {array_local_size.nkx, array_local_size.nky, array_local_size.nkz, 1, array_local_size.nl, array_local_size.ns};
    int start_coordinates_plus[6] = {0,0,0,array_local_size.nm-1,0,0};
    int start_coordinates_minus[6] = {0,0,0,0,0,0};
    mpi_sub_buf_size = array_local_size.nkx*array_local_size.nky*array_local_size.nkz*array_local_size.nl*array_local_size.ns;
    MPI_Type_create_subarray(6, dimensions_full, dimensions_sub, start_coordinates_plus, MPI_ORDER_C, MPI_C_DOUBLE_COMPLEX, &mpi_subarray_type_plus);
    MPI_Type_create_subarray(6, dimensions_full, dimensions_sub, start_coordinates_minus, MPI_ORDER_C, MPI_C_DOUBLE_COMPLEX, &mpi_subarray_type_minus);
    MPI_Type_commit(&mpi_subarray_type_plus);
    MPI_Type_commit(&mpi_subarray_type_minus);

}

void mpi_exchange_m_boundaries(COMPLEX *input_array, COMPLEX *plus_boundary, COMPLEX *minus_boundary){
    double start;
    double end;
    start = MPI_Wtime();
    MPI_Send(input_array, 1, mpi_subarray_type_minus, m_neighbour_ranks[MINUS], 0, mpi_cube_comm);
    MPI_Recv(plus_boundary, mpi_sub_buf_size, MPI_C_DOUBLE_COMPLEX, m_neighbour_ranks[PLUS], 0, mpi_cube_comm, MPI_STATUS_IGNORE);
    printf("[MPI process %d] received plus, t = %.2fs.! buf_size = %d\n", mpi_my_rank,MPI_Wtime()-start,mpi_sub_buf_size);
    start = MPI_Wtime();
    MPI_Send(input_array, 1, mpi_subarray_type_plus, m_neighbour_ranks[PLUS], 1, mpi_cube_comm);
    MPI_Recv(minus_boundary, mpi_sub_buf_size, MPI_C_DOUBLE_COMPLEX, m_neighbour_ranks[MINUS], 1, mpi_cube_comm, MPI_STATUS_IGNORE);
    printf("[MPI process %d] received minus, t = %.2fs.! buf_size = %d\n", mpi_my_rank,MPI_Wtime()-start,mpi_sub_buf_size);
}