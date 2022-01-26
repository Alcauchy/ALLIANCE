//
// Created by alcauchy on 27/12/2021.
//
#include "parameters_io.h"

void read_parameters(char *filename){
    MPI_Info info;
    MPI_File handle;
    char string[60];
    char *tmp[15];
    FILE *fp;
    fp = fopen(filename, "r");  /* open file for input */
    if (fp)  /* If no error occurred while opening file */
    {           /* input the data from the file. */
        while(fgets(string, 60, fp)){
             /* read the name from the file */
            sscanf(string, "%s : %*s", tmp);
            if(strcmp(tmp,"nkx")==0){
                sscanf(string, "%*s : %zu", &parameters.nkx);
                printf("[MPI process %d] nkx = %zu\n", mpi_my_rank,  parameters.nkx);
            }
            if(strcmp(tmp,"nky")==0){
                sscanf(string, "%*s : %zu", &parameters.nky);
                printf("[MPI process %d] nky = %zu\n", mpi_my_rank,  parameters.nky);
            }
            if(strcmp(tmp,"nz")==0){
                sscanf(string, "%*s : %zu", &parameters.nz);
                printf("[MPI process %d] nz = %zu\n", mpi_my_rank,  parameters.nz);
                parameters.nkz = parameters.nz/2+1;
            }
            if(strcmp(tmp,"nm")==0){
                sscanf(string, "%*s : %zu", &parameters.nm);
                printf("[MPI process %d] nm = %zu\n", mpi_my_rank,  parameters.nm);
            }
            if(strcmp(tmp,"nl")==0){
                sscanf(string, "%*s : %zu", &parameters.nl);
                printf("[MPI process %d] ns = %zu\n", mpi_my_rank, parameters.nl);
            }
            if(strcmp(tmp,"ns")==0){
                sscanf(string, "%*s : %zu", &parameters.ns);
                printf("[MPI process %d] ns = %zu\n", mpi_my_rank, parameters.ns);
            }
        }

    }

}