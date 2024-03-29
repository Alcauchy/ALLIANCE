/**************************************
* @file parameters_io.c
* \brief reads inpuit parameters from parameter file provided by user
*
***************************************/
////////////////////////////////////////////////////////////////////////////////
// 27/12/2021 created by Gene Gorbunov
//                                   PARAMETERS
//
// read_parameters
// read_parametersFromFile
//
// VERSION 1.0
////////////////////////////////////////////////////////////////////////////////

#include "parameters_io.h"
#include "utils_fftw.h"
#define VERBOSE 0
#define IO_RANK 0

struct system_param parameters;

/***************************************
 * \fn void init_global_size():
 * \brief initializes global size of the 6D array
 *
 * initializes <tt>array_local_size</tt> structure with global simulation size.
 ***************************************/
void init_global_size() {
    array_global_size.nkx = parameters.nkx;
    array_global_size.nky = parameters.nky;
    array_global_size.nkz = parameters.nz / 2 + 1;
    array_global_size.nx = parameters.nkx;
    array_global_size.ny = parameters.nky;
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
 * \fn read_parameters(char *filename)
 * \brief reads parameters from user parameter file.
 *
 * Reads parameters from user parameter file. All the parameters
 * are stored in the <tt>parameters</tt> structure
 ***************************************/
void read_parameters(char *filename) {
    char string[128];
    char tmp[32];
    int particle_index;
    FILE *fp;
    fp = fopen(filename, "r");  /* open file for input */
    if (mpi_my_rank == 0) printf("==============PARAMETERS==============\n");
    if (fp)  /* If no error occurred while opening file */
    {           /* input the data from the file. */
        while (fgets(string, 128, fp))
        {
            /* read the name from the file */
            sscanf(string, "%s : %*s", tmp);
            if (strcmp(tmp, "nkx") == 0)
            {
                sscanf(string, "%*s : %zu", &parameters.nkx);
                if (mpi_my_rank == IO_RANK) printf("nkx = %zu\n", parameters.nkx);
            }
            if (strcmp(tmp, "nky") == 0)
            {
                sscanf(string, "%*s : %zu", &parameters.nky);
                if (mpi_my_rank == IO_RANK) printf("nky = %zu\n", parameters.nky);
            }
            if (strcmp(tmp, "nz") == 0)
            {
                sscanf(string, "%*s : %zu", &parameters.nz);
                if (mpi_my_rank == IO_RANK) printf("nz = %zu\n", parameters.nz);
                parameters.nkz = parameters.nz / 2 + 1;

            }
            if (strcmp(tmp, "nm") == 0)
            {
                sscanf(string, "%*s : %zu", &parameters.nm);
                if (mpi_my_rank == IO_RANK)  printf("nm = %zu\n", parameters.nm);
            }
            if (strcmp(tmp, "nl") == 0)
            {
                sscanf(string, "%*s : %zu", &parameters.nl);
                if (mpi_my_rank == IO_RANK) printf("ns = %zu\n", parameters.nl);
            }
            if (strcmp(tmp, "ns") == 0)
            {
                sscanf(string, "%*s : %zu", &parameters.ns);
                if (mpi_my_rank == IO_RANK) printf("ns = %zu\n", parameters.ns);
                parameters.charge = malloc(parameters.ns * sizeof(parameters.charge));
                parameters.mass = malloc(parameters.ns * sizeof(parameters.mass));
                parameters.density = malloc(parameters.ns * sizeof(parameters.density));
                parameters.temperature = malloc(parameters.ns * sizeof(parameters.temperature));
            }
            if (strcmp(tmp, "dealiasing") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.dealiasing);
                if (mpi_my_rank == IO_RANK) printf("dealiasing rule is %d\n", parameters.dealiasing);

            }
            if (strcmp(tmp, "particle") == 0)
            {
                sscanf(string, "%*s : %d", &particle_index);
                if (particle_index >= parameters.ns)
                {
                    if (mpi_my_rank == IO_RANK) printf("[MPI process %d] Wrong number of particles (%d)! ns parameter = %d! \n Aborting... \n",
                           mpi_my_rank, particle_index, parameters.ns);
                    exit(1);
                }
                if (mpi_my_rank == IO_RANK) printf("particle index is %d\n", particle_index);
            }
            if (strcmp(tmp, "density") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.density[particle_index]);
                if (mpi_my_rank == IO_RANK) printf("particle %d, density = %3lf\n", particle_index,
                       parameters.density[particle_index]);
            }
            if (strcmp(tmp, "T") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.temperature[particle_index]);
                if (mpi_my_rank == IO_RANK) printf("particle %d, temperature = %3lf\n", particle_index,
                       parameters.temperature[particle_index]);
            }
            if (strcmp(tmp, "charge") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.charge[particle_index]);
                if (mpi_my_rank == IO_RANK) printf("particle %d, charge = %3lf\n", particle_index,
                       parameters.charge[particle_index]);
            }
            if (strcmp(tmp, "mass") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.mass[particle_index]);
                if (mpi_my_rank == IO_RANK) printf("particle %d, mass = %3lf\n", particle_index,
                       parameters.mass[particle_index]);
            }
            if (strcmp(tmp, "electromagnetic") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.electromagnetic);
                if (mpi_my_rank == IO_RANK) printf("electromagnetic = %d\n", parameters.electromagnetic);
            }
            if (strcmp(tmp, "adiabatic") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.adiabatic);
                if (mpi_my_rank == IO_RANK) printf("adiabatic = %d\n", parameters.adiabatic);
            }
            if (strcmp(tmp, "initial") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.initial);
                if (mpi_my_rank == IO_RANK) printf("initial conditions set to %d\n", parameters.initial);
            }
            if (strcmp(tmp, "simulation_name") == 0)
            {
                sscanf(string, "%*s : %s", &parameters.from_simulationName);
                if (mpi_my_rank == IO_RANK) printf("from simulation %s \n", parameters.from_simulationName);
            }
            if (strcmp(tmp, "beta") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.beta);
                if (mpi_my_rank == IO_RANK) printf("plasma beta = %f\n", parameters.beta);
            }
            if (strcmp(tmp, "dt") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.dt);
                if (mpi_my_rank == IO_RANK) printf("dt = %f\n",parameters.dt);
            }
            if (strcmp(tmp, "timesteps") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.Nt);
                if (mpi_my_rank == IO_RANK) printf("number of timesteps = %d\n", parameters.Nt);
            }
            if (strcmp(tmp, "compute_k") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.compute_k);
                if (mpi_my_rank == IO_RANK) printf("compute k spec = %d\n",parameters.compute_k);
            }
            if (strcmp(tmp, "compute_m") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.compute_m);
                if (mpi_my_rank == IO_RANK) printf("compute m spectra = %d\n",parameters.compute_m);
            }
            if (strcmp(tmp, "first_shell") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.firstShell);
                if (mpi_my_rank == IO_RANK) printf("first shell = %f\n",parameters.firstShell);
            }
            if (strcmp(tmp, "last_shell") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.lastShell);
                if (mpi_my_rank == IO_RANK) printf("last shell =  %f\n",parameters.lastShell);
            }
            if (strcmp(tmp, "iter_EMfield") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.iter_EMfield);
                if (mpi_my_rank == IO_RANK) printf("save field every %d timesteps\n",parameters.iter_EMfield);
            }
            if (strcmp(tmp, "checkpoints") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.checkpoints);
                if (mpi_my_rank == IO_RANK) printf("number of checkpoints = %d \n",parameters.checkpoints);
            }
            if (strcmp(tmp, "iter_checkpoint") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.iter_checkpoint);
                if (mpi_my_rank == IO_RANK) printf("save checkpoint every %d timesteps\n", parameters.iter_checkpoint);
            }
            if (strcmp(tmp, "iter_distribution") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.iter_distribution);
                if (mpi_my_rank == IO_RANK) printf("save distribution function every %d timesteps\n", parameters.iter_distribution);
            }
            if (strcmp(tmp, "save_distribution") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.save_distrib);
                if (mpi_my_rank == IO_RANK) printf("save distribution = %d\n", parameters.save_distrib);
            }
            if (strcmp(tmp, "iter_diagnostics") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.iter_diagnostics);
                if (mpi_my_rank == IO_RANK) printf("save free energy every %d timesteps\n", parameters.iter_diagnostics);
            }
            if (strcmp(tmp, "save_diagnostics") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.save_diagnostics);
                if (mpi_my_rank == IO_RANK) printf("save free energy = %d \n", parameters.save_diagnostics);
            }
            if (strcmp(tmp, "save_energy") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.save_diagnostics);
                if (mpi_my_rank == IO_RANK) printf("save free energy = %d \n", parameters.save_diagnostics);
            }
            if (strcmp(tmp, "save_EMfield") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.save_EMfield);
                if (mpi_my_rank == IO_RANK) printf("save electromagnetic fields = %d \n", parameters.save_EMfield);
            }
            if (strcmp(tmp, "Lx") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.Lx);
                if (mpi_my_rank == IO_RANK) printf("Lx = %f \n",parameters.Lx);
            }
            if (strcmp(tmp, "Ly") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.Ly);
                if (mpi_my_rank == IO_RANK) printf("Ly = %f \n",parameters.Ly);
            }
            if (strcmp(tmp, "Lz") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.Lz);
                if (mpi_my_rank == IO_RANK) printf("Lz = %f \n", parameters.Lz);
            }
            if (strcmp(tmp, "nproc_m") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.nproc_m);
                if (mpi_my_rank == IO_RANK) printf("num of processes in m used is %d \n", parameters.nproc_m);
            }
            if (strcmp(tmp, "nproc_k") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.nproc_k);
                if (mpi_my_rank == IO_RANK) printf("num of processes in k used is %d \n",parameters.nproc_k);
            }
            if (strcmp(tmp, "iter_dt") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.iter_dt);
                if (mpi_my_rank == IO_RANK) printf("update dt every %d steps \n", parameters.iter_dt);
            }
            if (strcmp(tmp, "dissip_k") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.mu_k);
                if (mpi_my_rank == IO_RANK) printf("dissipation in k is %f \n", parameters.mu_k);
            }
            if (strcmp(tmp, "dissip_m") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.mu_m);
                if (mpi_my_rank == IO_RANK) printf("dissipation in m is %f \n", parameters.mu_m);
            }
            if (strcmp(tmp, "dt_dissip") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.dissipDt);
                if (mpi_my_rank == IO_RANK) printf("dissipation time step is %f \n", parameters.dissipDt);
            }
            if (strcmp(tmp, "dt_linear") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.linDt);
                if (mpi_my_rank == IO_RANK) printf("linear time step is %f \n", parameters.linDt);
            }
            if (strcmp(tmp, "k_max") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.forceKmax);
                if (mpi_my_rank == IO_RANK) printf("forcing kmax = %f \n",parameters.forceKmax);
            }
            if (strcmp(tmp, "k_min") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.forceKmin);
                if (mpi_my_rank == IO_RANK) printf("forcing k_min = %f \n",parameters.forceKmin);
            }
            if (strcmp(tmp, "power") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.forcePower);
                if (mpi_my_rank == IO_RANK) printf("forcing power = %f \n",parameters.forcePower);
            }
            if (strcmp(tmp, "spectrum_type") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.spectrum);
                if (mpi_my_rank == IO_RANK) printf("spectrum type = %d \n",parameters.spectrum);
            }
            if (strcmp(tmp, "k_unit") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.unitK);
                if (mpi_my_rank == IO_RANK) printf("unit k = %lf \n", parameters.unitK);
            }
            if (strcmp(tmp, "allow_rescale") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.allow_rescale);
                if (mpi_my_rank == IO_RANK) printf("allow_rescale = %d \n", parameters.allow_rescale);
            }
            if (strcmp(tmp, "dissip_kz") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.mu_kz);
                if (mpi_my_rank == IO_RANK) printf("dissip_kz = %f \n", parameters.mu_kz);
            }
            if (strcmp(tmp, "hyper_k") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.lap_k);
                if (mpi_my_rank == IO_RANK) printf("hyperlaplacian power for k_perp = %f \n", parameters.lap_k);
            }
            if (strcmp(tmp, "hyper_kz") == 0)
            {
                sscanf(string, "%*s : %lf", &parameters.lap_kz);
                if (mpi_my_rank == IO_RANK) printf("hyperlaplacian power for k_perp = %f \n", parameters.lap_kz);
            }
            if (strcmp(tmp, "compute_nonlinear") == 0)
            {
                sscanf(string, "%*s : %d", &parameters.compute_nonlinear);
                if (mpi_my_rank == IO_RANK) printf("compute nonlinear flux every %d steps\n",&parameters.compute_nonlinear);
            }
        }
    }
    if (mpi_my_rank == 0) printf("=====================================\n");
    init_global_size();
   // if (parameters.initial == 1)
   // {
    //    read_parametersFromFile(parameters.from_simulationName);
   // }
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
