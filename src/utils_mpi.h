#ifndef ALLIANCE_ALPHA_1_0_UTILS_MPI_H
#define ALLIANCE_ALPHA_1_0_UTILS_MPI_H

#include <stdio.h> // standard C lib
#include <stdlib.h> //Standard C lib
#include <mpi.h> // MPI
#include "array.h"
#include "parameters_io.h"

void mpi_init();                        //initialize mpi environment
void mpi_kill();                        //finalize mpi environment
void mpi_generateTopology();
void mpi_createTopology();             //create 2D topology
void mpi_findHermiteNeighbours();     // finds neighbour processes in hermite direction
void mpi_splitInRows();               // splits communicators in rows. needed for fftw transforms along kx direction of parallelization
void mpi_splitInCols();               // splits communicator in columns. Needed to send data to compute fields for each process.
void mpi_initMExchange();             // used to define all the variables needed to prepare exchange of m+1 and m-1 Hermite moments between processors
void mpi_exchangeMBoundaries(COMPLEX *input_array, COMPLEX *plus_boundary, COMPLEX *minus_boundary);       // required to exchange the m+1 and m-1 Hermite moments between processes
void mpi_exchangeMBoundaries_r(double *input_array, double *plus_boundary, double *minus_boundary);       // required to exchange the m+1 and m-1 Hermite moments between processes for real data
void mpi_getLocalArraySize();        //computes array sizes stored locally by each process
void mpi_getLocalArrayOffsets();     //computes offsets for easy 6d indexing (see function get_flat_r)
void mpi_sendVector(COMPLEX *from_array, COMPLEX *to_buffer, int from_proc, int to_proc);               // send mpi_vector_kxSlice to a chosen processor

//from mpi_init.c file
extern int mpi_my_rank;                // rank of the process for 2D MPI communicator
extern int mpi_size;                   // size of the communicator
extern int mpi_my_row_rank;            // rank of the process in a kx direction of parallelization
extern int mpi_my_col_rank;            // rank of the process in Hermite direction of parallelization
extern int mpi_my_coords[2];           // coordinates of the process in 2D topology
extern int mpi_dims[];                 // size of the dimensions, 0th dimension is for Hermite parallelization, and 1st dimension is for kx.
extern int *mpi_whereIsX;
extern int *mpi_whereIsM;
extern int *mpi_whereIsY;
extern MPI_Comm mpi_cube_comm;         // 2D topology communicator
extern MPI_Comm mpi_row_comm;          // row communicator (kx direction)
extern MPI_Comm mpi_col_comm;          // column communicator (Hermite direction)

#endif //ALLIANCE_ALPHA_1_0_UTILS_MPI_H

