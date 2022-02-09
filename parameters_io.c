//
// Created by alcauchy on 27/12/2021.
//

#include "parameters_io.h"
#include "fftw_utils.h"


struct system_param parameters;

void init_global_size() {
    array_global_size.nkx = parameters.nkx;
    array_global_size.nky = parameters.nky;
    array_global_size.nkz = parameters.nz / 2 + 1;
    array_global_size.nz = parameters.nz;
    array_global_size.nm = parameters.nm;
    array_global_size.nl = parameters.nl;
    array_global_size.ns = parameters.ns;
    array_global_size.total_comp = array_global_size.nkx *
                                   array_global_size.nky *
                                   array_global_size.nkz *
                                   array_global_size.nm *
                                   array_global_size.nl *
                                   array_global_size.ns;
    array_global_size.total_real = array_global_size.nkx *
                                   array_global_size.nky *
                                   array_global_size.nz *
                                   array_global_size.nm *
                                   array_global_size.nl *
                                   array_global_size.ns;
};

void read_parameters(char *filename) {
    char string[128];
    char tmp[16];
    int particle_index;
    FILE *fp;
    fp = fopen(filename, "r");  /* open file for input */
    if (fp)  /* If no error occurred while opening file */
    {           /* input the data from the file. */
        while (fgets(string, 128, fp))
        {
            /* read the name from the file */
            sscanf(string, "%s : %*s", tmp);
            if (strcmp(tmp, "nkx") == 0)
            {
                sscanf(string, "%*s : %zu", &parameters.nkx);
                printf("[MPI process %d] nkx = %zu\n", mpi_my_rank, parameters.nkx);
            }
            if (strcmp(tmp, "nky") == 0)
            {
                sscanf(string, "%*s : %zu", &parameters.nky);
                printf("[MPI process %d] nky = %zu\n", mpi_my_rank, parameters.nky);
            }
            if (strcmp(tmp, "nz") == 0)
            {
                sscanf(string, "%*s : %zu", &parameters.nz);
                printf("[MPI process %d] nz = %zu\n", mpi_my_rank, parameters.nz);

            }
            if (strcmp(tmp, "nm") == 0)
            {
                sscanf(string, "%*s : %zu", &parameters.nm);
                printf("[MPI process %d] nm = %zu\n", mpi_my_rank, parameters.nm);
            }
            if (strcmp(tmp, "nl") == 0)
            {
                sscanf(string, "%*s : %zu", &parameters.nl);
                printf("[MPI process %d] ns = %zu\n", mpi_my_rank, parameters.nl);
            }
            if (strcmp(tmp, "ns") == 0)
            {
                sscanf(string, "%*s : %zu", &parameters.ns);
                printf("[MPI process %d] ns = %zu\n", mpi_my_rank, parameters.ns);
                parameters.charge = malloc(parameters.ns * sizeof(parameters.charge));
                parameters.mass = malloc(parameters.ns * sizeof(parameters.mass));
                parameters.density = malloc(parameters.ns * sizeof(parameters.density));
                parameters.temperature = malloc(parameters.ns * sizeof(parameters.temperature));
            }
            if (strcmp(tmp, "dealiasing") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.dealiasing);
                if (parameters.dealiasing == TWOTHIRDS)
                {
                    printf("[MPI process %d] dealiasing rule is %d\n", mpi_my_rank, parameters.dealiasing);
                    fftw_dealiasing = dealiasing23;
                }
            }
            if (strcmp(tmp, "particle") == 0)
            {
                sscanf(string, "%*s : %d", &particle_index);
                if (particle_index >= parameters.ns)
                {
                    printf("[MPI process %d] Wrong number of particles (%d)! ns parameter = %d! \n Aborting... \n",
                           mpi_my_rank, particle_index, parameters.ns);
                    exit(1);
                }
                printf("[MPI process %d] particle index is %d\n", mpi_my_rank, particle_index);
            }
            if (strcmp(tmp, "density") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.density[particle_index]);
                printf("[MPI process %d] particle %d, density = %3lf\n", mpi_my_rank, particle_index,
                       parameters.density[particle_index]);
            }
            if (strcmp(tmp, "temperature") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.temperature[particle_index]);
                printf("[MPI process %d] particle %d, temperature = %3lf\n", mpi_my_rank, particle_index,
                       parameters.temperature[particle_index]);
            }
            if (strcmp(tmp, "charge") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.charge[particle_index]);
                printf("[MPI process %d] particle %d, charge = %3lf\n", mpi_my_rank, particle_index,
                       parameters.charge[particle_index]);
            }
            if (strcmp(tmp, "mass") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.mass[particle_index]);
                printf("[MPI process %d] particle %d, mass = %3lf\n", mpi_my_rank, particle_index,
                       parameters.mass[particle_index]);
            }
            if (strcmp(tmp, "electromagnetic") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.electromagnetic);
                printf("[MPI process %d] electromagnetic = %d\n", mpi_my_rank, parameters.electromagnetic);
            }
            if (strcmp(tmp, "adiabatic") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.adiabatic);
                printf("[MPI process %d] adiabatic = %d\n", mpi_my_rank, parameters.adiabatic);
            }
            if (strcmp(tmp, "beta") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.beta);
                printf("[MPI process %d] plasma beta = %f\n", mpi_my_rank, parameters.beta);
            }
        }
    }
    init_global_size();
}

