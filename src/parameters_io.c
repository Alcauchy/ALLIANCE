//
// Created by alcauchy on 27/12/2021.
//

#include "parameters_io.h"
#include "fftw_utils.h"
#define VERBOSE 0

struct system_param parameters;

/***************************************
 * init_global_size():
 ***************************************/
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

/***************************************
 * read_parameters(char *filename):
 ***************************************/
void read_parameters(char *filename) {
    char string[128];
    char tmp[32];
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
                parameters.nkz = parameters.nz / 2 + 1;

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
            if (strcmp(tmp, "T") == 0)
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
            if (strcmp(tmp, "initial") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.initial);
                printf("[MPI process %d] initial conditions set to %d\n", mpi_my_rank, parameters.initial);
            }
            if (strcmp(tmp, "simulation_name") == 0)
            {
                sscanf(string, "%*s : %s", &parameters.from_simulationName);
                printf("[MPI process %d] from simulation %s \n", mpi_my_rank, parameters.from_simulationName);
            }
            if (strcmp(tmp, "beta") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.beta);
                printf("[MPI process %d] plasma beta = %f\n", mpi_my_rank, parameters.beta);
            }
            if (strcmp(tmp, "dt") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.dt);
                printf("[MPI process %d] dt = %f\n", mpi_my_rank, parameters.dt);
            }
            if (strcmp(tmp, "timesteps") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.Nt);
                printf("[MPI process %d] number of timesteps = %d\n", mpi_my_rank, parameters.Nt);
            }
            if (strcmp(tmp, "compute_k") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.compute_k);
                printf("[MPI process %d] compute k spec = %d\n", mpi_my_rank, parameters.compute_k);
            }
            if (strcmp(tmp, "k_shells") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.k_shells);
                printf("[MPI process %d] number of k shells = %d\n", mpi_my_rank, parameters.k_shells);
            }
            if (strcmp(tmp, "compute_m") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.compute_m);
                printf("[MPI process %d] compute m spectra = %d\n", mpi_my_rank, parameters.compute_m);
            }
            if (strcmp(tmp, "compute_k_every") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.compute_k_every);
                printf("[MPI process %d] compute k spectra every %d timesteps\n", mpi_my_rank, parameters.compute_k_every);
            }
            if (strcmp(tmp, "compute_m_every") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.compute_m_every);
                printf("[MPI process %d] compute m spectra every %d timesteps\n", mpi_my_rank, parameters.compute_m_every);
            }
            if (strcmp(tmp, "first_shell") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.firstShell);
                printf("[MPI process %d] first shell = %f\n", mpi_my_rank, parameters.firstShell);
            }
            if (strcmp(tmp, "last_shell") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.lastShell);
                printf("[MPI process %d] last shell =  %f\n", mpi_my_rank, parameters.lastShell);
            }
            if (strcmp(tmp, "save_field_every") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.save_field_step);
                printf("[MPI process %d] save field every %d timesteps\n", mpi_my_rank, parameters.save_field_step);
            }
            if (strcmp(tmp, "checkpoints") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.checkpoints);
                printf("[MPI process %d] number of checkpoints = %d \n", mpi_my_rank, parameters.checkpoints);
            }
            if (strcmp(tmp, "save_checkpoint_every") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.save_checkpoint_step);
                printf("[MPI process %d] save checkpoint every %d timesteps\n", mpi_my_rank, parameters.save_checkpoint_step);
            }
            if (strcmp(tmp, "save_distribution_every") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.save_distrib_step);
                printf("[MPI process %d] save distribution function every %d timesteps\n", mpi_my_rank, parameters.save_distrib_step);
            }
            if (strcmp(tmp, "save_distribution") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.save_distrib);
                printf("[MPI process %d] save distribution = %d\n", mpi_my_rank, parameters.save_distrib);
            }
            if (strcmp(tmp, "save_k_spec_every") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.save_kSpec_step);
                printf("[MPI process %d] save k spectrum every %d timesteps\n", mpi_my_rank, parameters.save_kSpec_step);
            }
            if (strcmp(tmp, "save_m_spec_every") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.save_mSpec_step);
                printf("[MPI process %d] save m spectrum every %d timesteps\n", mpi_my_rank, parameters.save_mSpec_step);
            }
            if (strcmp(tmp, "save_energy_every") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.save_energy_step);
                printf("[MPI process %d] save free energy every %d timesteps\n", mpi_my_rank, parameters.save_energy_step);
            }
            if (strcmp(tmp, "save_k_spec") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.save_kSpec);
                printf("[MPI process %d] save k spectrum = %d \n", mpi_my_rank, parameters.save_kSpec);
            }
            if (strcmp(tmp, "save_m_spec") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.save_mSpec);
                printf("[MPI process %d] save m spectrum = %d \n", mpi_my_rank, parameters.save_mSpec);
            }
            if (strcmp(tmp, "save_energy") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.save_energy);
                printf("[MPI process %d] save free energy = %d \n", mpi_my_rank, parameters.save_energy);
            }
            if (strcmp(tmp, "save_energy") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.save_energy);
                printf("[MPI process %d] save free energy = %d \n", mpi_my_rank, parameters.save_energy);
            }
            if (strcmp(tmp, "save_field") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.save_field);
                printf("[MPI process %d] save free energy = %d \n", mpi_my_rank, parameters.save_field);
            }

            if (strcmp(tmp, "save_directory") == 0)
            {
                sscanf(string, "%*s : %s", &parameters.save_dir);
                printf("[MPI process %d] save to %s directory \n", mpi_my_rank, parameters.save_dir);
            }
        }
    }
    init_global_size();
    if (parameters.initial == 1)
    {
        read_parametersFromFile(parameters.from_simulationName);
    }
}

/***************************************
 * read_parametersFromFile(char *filename):
 ***************************************/
void read_parametersFromFile(char *filename){
    hid_t dset_id, dspace_id, file_id, filespace, memspace;
    hid_t plist_id; //property list id
    MPI_Info info = MPI_INFO_NULL;
    hid_t params_rank = 1;
    hid_t param_dims[1] = {1};

    plist_id = H5Pcreate(H5P_FILE_ACCESS); // access property list
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
    /* open file to read */
    file_id = H5Fopen(filename,H5F_ACC_RDONLY,plist_id);
    H5Pclose(plist_id);
    /* read beta parameter */
    dset_id = H5Dopen2(file_id, "beta", H5P_DEFAULT);
    dspace_id = H5Dget_space(dset_id);
    memspace = H5Screate_simple(params_rank,param_dims,param_dims);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dread(dset_id,H5T_NATIVE_DOUBLE, memspace, dspace_id, plist_id, &parameters.beta);
    H5Sclose(dspace_id);
    H5Sclose(memspace);
    H5Dclose(dset_id);

    /* read nkx parameter */
    dset_id = H5Dopen2(file_id, "nkx", H5P_DEFAULT);
    dspace_id = H5Dget_space(dset_id);
    memspace = H5Screate_simple(params_rank,param_dims,param_dims);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dread(dset_id,H5T_NATIVE_INT, memspace, dspace_id, plist_id, &parameters.nkx);
    H5Sclose(dspace_id);
    H5Sclose(memspace);
    H5Dclose(dset_id);
    if(VERBOSE) printf("[MPI process %d] READ nkx = %d \n", mpi_my_rank, parameters.nkx);

    /* read nky parameter */
    dset_id = H5Dopen2(file_id, "nky", H5P_DEFAULT);
    dspace_id = H5Dget_space(dset_id);
    memspace = H5Screate_simple(params_rank,param_dims,param_dims);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dread(dset_id,H5T_NATIVE_INT, memspace, dspace_id, plist_id, &parameters.nky);
    if(VERBOSE) printf("[MPI process %d] READ nky = %d \n", mpi_my_rank, parameters.nky);
    H5Sclose(dspace_id);
    H5Sclose(memspace);
    H5Dclose(dset_id);

    /* read nkz parameter */
    dset_id = H5Dopen2(file_id, "nkz", H5P_DEFAULT);
    dspace_id = H5Dget_space(dset_id);
    memspace = H5Screate_simple(params_rank,param_dims,param_dims);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dread(dset_id,H5T_NATIVE_INT, memspace, dspace_id, plist_id, &parameters.nkz);
    if(VERBOSE) printf("[MPI process %d] READ nkz = %d \n", mpi_my_rank, parameters.nkz);
    H5Sclose(dspace_id);
    H5Sclose(memspace);
    H5Dclose(dset_id);

    /* read nm parameter */
    dset_id = H5Dopen2(file_id, "nm", H5P_DEFAULT);
    dspace_id = H5Dget_space(dset_id);
    memspace = H5Screate_simple(params_rank,param_dims,param_dims);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dread(dset_id,H5T_NATIVE_INT, memspace, dspace_id, plist_id, &parameters.nm);
    if(VERBOSE) printf("[MPI process %d] READ nm = %d \n", mpi_my_rank, parameters.nm);
    H5Sclose(dspace_id);
    H5Sclose(memspace);
    H5Dclose(dset_id);

    /* read nl parameter */
    dset_id = H5Dopen2(file_id, "nl", H5P_DEFAULT);
    dspace_id = H5Dget_space(dset_id);
    memspace = H5Screate_simple(params_rank,param_dims,param_dims);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dread(dset_id,H5T_NATIVE_INT, memspace, dspace_id, plist_id, &parameters.nl);
    if(VERBOSE) printf("[MPI process %d] READ nl = %d \n", mpi_my_rank, parameters.nl);
    H5Sclose(dspace_id);
    H5Sclose(memspace);
    H5Dclose(dset_id);

    /* read ns parameter */
    dset_id = H5Dopen2(file_id, "ns", H5P_DEFAULT);
    dspace_id = H5Dget_space(dset_id);
    memspace = H5Screate_simple(params_rank,param_dims,param_dims);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dread(dset_id,H5T_NATIVE_INT, memspace, dspace_id, plist_id, &parameters.ns);
    if(VERBOSE) printf("[MPI process %d] READ ns = %d \n", mpi_my_rank, parameters.ns);
    H5Sclose(dspace_id);
    H5Sclose(memspace);
    H5Dclose(dset_id);

    /*read particle parameters*/
    int rank_particle = 1;
    hsize_t dims_particle[1] = {parameters.ns};
    /* read charge information */
    dspace_id = H5Screate_simple(rank_particle, dims_particle, NULL);
    dset_id = H5Dopen2(file_id, "charge", H5P_DEFAULT);
    memspace = H5Screate_simple(rank_particle, dims_particle, NULL);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dread(dset_id,H5T_NATIVE_DOUBLE, memspace, dspace_id, plist_id, parameters.charge);
    if(VERBOSE)
    {
        for(size_t i = 0; i < parameters.ns; i++) printf("[MPI process %d] READ q[%d] = %f \n", mpi_my_rank, i, parameters.charge[i]);
    }
    H5Dclose(dset_id);
    H5Sclose(dspace_id);
    H5Sclose(memspace);
    /* read temperature information */
    dspace_id = H5Screate_simple(rank_particle, dims_particle, NULL);
    dset_id = H5Dopen2(file_id, "temperature", H5P_DEFAULT);
    memspace = H5Screate_simple(rank_particle, dims_particle, NULL);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dread(dset_id,H5T_NATIVE_DOUBLE, memspace, dspace_id, plist_id, parameters.temperature);
    if(VERBOSE)
    {
        for(size_t i = 0; i < parameters.ns; i++) printf("[MPI process %d] READ T[%d] = %f \n", mpi_my_rank, i, parameters.temperature[i]);
    }
    H5Dclose(dset_id);
    H5Sclose(dspace_id);
    H5Sclose(memspace);

    /* read mass information */
    dspace_id = H5Screate_simple(rank_particle, dims_particle, NULL);
    dset_id = H5Dopen2(file_id, "mass", H5P_DEFAULT);
    memspace = H5Screate_simple(rank_particle, dims_particle, NULL);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dread(dset_id,H5T_NATIVE_DOUBLE, memspace, dspace_id, plist_id, parameters.mass);
    if(VERBOSE)
    {
        for(size_t i = 0; i < parameters.ns; i++) printf("[MPI process %d] READ m[%d] = %f \n", mpi_my_rank, i, parameters.mass[i]);
    }
    H5Dclose(dset_id);
    H5Sclose(dspace_id);
    H5Sclose(memspace);

    /* read density information */
    dspace_id = H5Screate_simple(rank_particle, dims_particle, NULL);
    dset_id = H5Dopen2(file_id, "density", H5P_DEFAULT);
    memspace = H5Screate_simple(rank_particle, dims_particle, NULL);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dread(dset_id,H5T_NATIVE_DOUBLE, memspace, dspace_id, plist_id, parameters.density);
    if(VERBOSE)
    {
        for(size_t i = 0; i < parameters.ns; i++) printf("[MPI process %d] READ n[%d] = %f \n", mpi_my_rank, i, parameters.density[i]);
    }
    H5Dclose(dset_id);
    H5Sclose(dspace_id);
    H5Sclose(memspace);

    H5Fclose(file_id);


};
