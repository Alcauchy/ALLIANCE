//
// Created by alcauchy on 16/12/2021.
//

#include "utils_hdf.h"
#include <unistd.h>
#include <sys/stat.h>

#define BASE_DIR "."
#define WORK_DIR "."
#define CHCK_DIR "checkpoint"
#define CHCK_NAME "chk"
#if defined(WIN32) || defined(_WIN32)
#define PATH_SEPARATOR "\\"
#else
#define PATH_SEPARATOR "/"
#endif
#define VERBOSE 1
#define FILENAME_ID_LEN 128
#define CHECKPOINT_ROOT 0
#define PATH_LEN 128

int hdf_rank = 6;
int hdf_rankFields = 3;
int hdf_freeEnergyCalls = 0;
char **hdf_checkpointNames;
char hdf_newCheckpointName[FILENAME_ID_LEN];
char SIMULATION_PATH[PATH_LEN];
char CHECKPOINT_PATH[PATH_LEN];
char PARAMETER_FILENAME[FILENAME_ID_LEN];
char DISTRIBUTION_FILENAME[FILENAME_ID_LEN];
char FIELD_FILENAME[FILENAME_ID_LEN];
size_t hdf_checkpointCount = 0;
hid_t complex_id;
hsize_t dataspace_dims_r[6];
hsize_t dataspace_dims_c[6];
hsize_t dataspace_dimsFields[3];
hsize_t chunk_dims_r[6];
hsize_t chunk_dims_c[6];
hsize_t chunk_dimsFields[3];
hsize_t offset[6];
hsize_t offsetFields[3];
hsize_t count[6] = {1,1,1,1,1,1};
hsize_t stride[6] = {1,1,1,1,1,1};
hsize_t countFields[3] = {1,1,1};
hsize_t strideFields[3] = {1,1,1};
herr_t	status;
MPI_Info info = MPI_INFO_NULL;
complex_t tmp;  //used only to compute offsets

void complex_t_init(){
    complex_id = H5Tcreate(H5T_COMPOUND, sizeof(tmp));
    H5Tinsert(complex_id, "r", HOFFSET(complex_t,re), H5T_NATIVE_DOUBLE);
    H5Tinsert(complex_id, "i", HOFFSET(complex_t,im), H5T_NATIVE_DOUBLE);
}

/***************************
 *
 *  INITIALIZE HDF5
 *
 * *************************/

/***************************
 *  hdf_init
 * *************************/
void hdf_init(){
    complex_t_init();
    hdf_initField();
    hdf_createSaveDirs();
    if (parameters.checkpoints > 0)
    {
        hdf_initCheckpoints();
    }
    dataspace_dims_r[0] = array_global_size.nkx;
    dataspace_dims_r[1] = array_global_size.nky;
    dataspace_dims_r[2] = array_global_size.nz + 2;
    dataspace_dims_r[3] = array_global_size.nm;
    dataspace_dims_r[4] = array_global_size.nl;
    dataspace_dims_r[5] = array_global_size.ns;

    dataspace_dims_c[0] = array_global_size.nkx;
    dataspace_dims_c[1] = array_global_size.nky;
    dataspace_dims_c[2] = array_global_size.nkz;
    dataspace_dims_c[3] = array_global_size.nm;
    dataspace_dims_c[4] = array_global_size.nl;
    dataspace_dims_c[5] = array_global_size.ns;

    chunk_dims_r[0] = array_local_size.nkx;
    chunk_dims_r[1] = array_local_size.nky;
    chunk_dims_r[2] = array_local_size.nz+2;
    chunk_dims_r[3] = array_local_size.nm;
    chunk_dims_r[4] = array_local_size.nl;
    chunk_dims_r[5] = array_local_size.ns;

    chunk_dims_c[0] = array_local_size.nkx;
    chunk_dims_c[1] = array_local_size.nky;
    chunk_dims_c[2] = array_local_size.nkz;
    chunk_dims_c[3] = array_local_size.nm;
    chunk_dims_c[4] = array_local_size.nl;
    chunk_dims_c[5] = array_local_size.ns;
    offset[0] = array_local_size.nkx * mpi_my_coords[1];
    offset[1] = 0;
    offset[2] = 0;
    offset[3] = array_local_size.nm * mpi_my_coords[0];
    offset[4] = 0;
    offset[5] = 0;

    printf("[MPI process %d] my coords are (%d, %d), offsets are (%d,%d,%d,%d,%d,%d)\n",mpi_my_rank,
           mpi_my_coords[0],
           mpi_my_coords[1],
           offset[0],
           offset[1],
           offset[2],
           offset[3],
           offset[4],
           offset[5]);
    printf("[MPI process %d] my coords are (%d, %d), chunk dims for real are (%d,%d,%d,%d,%d,%d)\n",mpi_my_rank,
           mpi_my_coords[0],
           mpi_my_coords[1],
           chunk_dims_r[0],
           chunk_dims_r[1],
           chunk_dims_r[2],
           chunk_dims_r[3],
           chunk_dims_r[4],
           chunk_dims_r[5]);
    printf("[MPI process %d] my coords are (%d, %d), chunk dims for complex are (%d,%d,%d,%d,%d,%d)\n",mpi_my_rank,
           mpi_my_coords[0],
           mpi_my_coords[1],
           chunk_dims_c[0],
           chunk_dims_c[1],
           chunk_dims_c[2],
           chunk_dims_c[3],
           chunk_dims_c[4],
           chunk_dims_c[5]);
};

/***************************
 *  hdf_createSaveDirs
 * *************************/
void hdf_createSaveDirs() {
    sprintf(SIMULATION_PATH, "%s%s", BASE_DIR,PATH_SEPARATOR);
    sprintf(CHECKPOINT_PATH, "%s%s%s%s", BASE_DIR,PATH_SEPARATOR,CHCK_DIR,PATH_SEPARATOR);
    sprintf(PARAMETER_FILENAME, "%s%s", SIMULATION_PATH,"parameters.h5");
    if (VERBOSE) printf("SIMULATION PATH = %s\n",SIMULATION_PATH);
    if (VERBOSE) printf("CHECKPOINT PATH = %s\n",CHECKPOINT_PATH);
    if (VERBOSE) printf("PARAMETER FILENAME = %s\n",PARAMETER_FILENAME);
    if(mpi_my_rank == CHECKPOINT_ROOT)
    {
        if (access(CHECKPOINT_PATH, F_OK) == 0)
        {
            //printf("directory %s already exists, please ensure you don't need it and delete it manually. Aborting...\n",SIMULATION_PATH);
            //exit(1);
        }
        else
        {
           // mkdir(SIMULATION_PATH,0777);
            mkdir(CHECKPOINT_PATH,0777);
        }
    }
    if(parameters.save_EMfield)
    {
        sprintf(FIELD_FILENAME, "%s%s", SIMULATION_PATH,"fields.h5");
        if (VERBOSE) printf("FIELD FILENAME = %s\n",FIELD_FILENAME);
    }
    if(parameters.save_distrib)
    {
        sprintf(DISTRIBUTION_FILENAME, "%s%s", SIMULATION_PATH,"h_");
        if (VERBOSE) printf("H FILENAME = %s\n",DISTRIBUTION_FILENAME);
    }
};

/***************************
 *  hdf_create_file_c
 * *************************/
void hdf_create_file_c(char *filename, COMPLEX *data){
    // now we will create a file, and write a dataset into it.
    hid_t file_id, dset_id;
    hid_t file_space, memory_space;
    hid_t plist_id; //property list id
    plist_id = H5Pcreate(H5P_FILE_ACCESS); // access property list
    H5Pset_fapl_mpio(plist_id, mpi_cube_comm, info);
    printf("[process id %d] trying to create a file\n",mpi_my_rank);
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id); //creating new file collectively
    printf("[process id %d] file created\n",mpi_my_rank);
    H5Pclose(plist_id);
    file_space = H5Screate_simple(hdf_rank, dataspace_dims_c, NULL);
    memory_space = H5Screate_simple(hdf_rank, chunk_dims_c, NULL);

    plist_id = H5Pcreate(H5P_DATASET_CREATE); //creating chunked dataset
    H5Pset_chunk(plist_id, hdf_rank, chunk_dims_c);
    dset_id = H5Dcreate(file_id,"g", complex_id, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    H5Sclose(file_space );

    file_space = H5Dget_space(dset_id);
    status = H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset, stride, count, chunk_dims_c);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    status = H5Dwrite(dset_id, complex_id, memory_space, file_space,plist_id, data);

    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Sclose(memory_space);
    H5Pclose(plist_id);
    H5Fclose(file_id);
}

/***************************
 *  hdf_create_file_r
 * *************************/
void hdf_create_file_r(char *filename, double *data){
    // now we will create a file, and write a dataset into it.
    hid_t file_id, dset_id;
    hid_t file_space, memory_space;
    hid_t plist_id; //property list id
    plist_id = H5Pcreate(H5P_FILE_ACCESS); // access property list
    H5Pset_fapl_mpio(plist_id, mpi_cube_comm, info);
    printf("[process id %d] trying to create a file\n",mpi_my_rank);
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id); //creating new file collectively
    printf("[process id %d] file created\n",mpi_my_rank);
    H5Pclose(plist_id);
    file_space = H5Screate_simple(hdf_rank, dataspace_dims_r, NULL);
    memory_space = H5Screate_simple(hdf_rank, chunk_dims_r, NULL);

    plist_id = H5Pcreate(H5P_DATASET_CREATE); //creating chunked dataset
    H5Pset_chunk(plist_id, hdf_rank, chunk_dims_r);
    dset_id = H5Dcreate(file_id,"g", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    H5Sclose(file_space);

    file_space = H5Dget_space(dset_id);
    status = H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset, stride, count, chunk_dims_r);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memory_space, file_space,plist_id, data);

    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Sclose(memory_space);
    H5Pclose(plist_id);
    H5Fclose(file_id);
}

/***************************
 *  hdf_initField
 * *************************/
void hdf_initField(){
    dataspace_dimsFields[0] = array_global_size.nkx;
    dataspace_dimsFields[1] = array_global_size.nky;
    dataspace_dimsFields[2] = array_global_size.nkz;

    chunk_dimsFields[0] = array_local_size.nkx;
    chunk_dimsFields[1] = array_local_size.nky;
    chunk_dimsFields[2] = array_local_size.nkz;
    printf("[MPI process %d] %zu %zu %zu\n",mpi_my_rank, chunk_dimsFields[0],chunk_dimsFields[1],chunk_dimsFields[2]);

    offsetFields[0] = array_local_size.nkx * mpi_my_coords[1];
    offsetFields[1] = 0;
    offsetFields[2] = 0;

};

/***************************
 *  hdf_saveFieldA
 * *************************/
void hdf_saveFieldA(char *filename){
    // now we will create a file, and write a dataset into it.
    hid_t file_id, dset_id;
    hid_t file_space, memory_space;
    hid_t plist_id; //property list id
    plist_id = H5Pcreate(H5P_FILE_ACCESS); // access property list
    H5Pset_fapl_mpio(plist_id, mpi_row_comm, info);
    printf("[process id %d] trying to create a file to write fields\n",mpi_my_rank);
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id); //creating new file collectively
    printf("[process id %d] file created\n",mpi_my_rank);
    H5Pclose(plist_id);
    file_space = H5Screate_simple(hdf_rankFields, dataspace_dimsFields, NULL);
    memory_space = H5Screate_simple(hdf_rankFields, chunk_dimsFields, NULL);

    plist_id = H5Pcreate(H5P_DATASET_CREATE); //creating chunked dataset
    H5Pset_chunk(plist_id, hdf_rankFields, chunk_dimsFields);
    dset_id = H5Dcreate(file_id,"A", complex_id, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    H5Sclose(file_space);

    file_space = H5Dget_space(dset_id);
    status = H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offsetFields, strideFields, countFields, chunk_dimsFields);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    status = H5Dwrite(dset_id, complex_id, memory_space, file_space,plist_id, fields_fields.A);

    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Sclose(memory_space);
    H5Pclose(plist_id);
    H5Fclose(file_id);

};

/***************************
 *  hdf_saveFieldB
 * *************************/
void hdf_saveFieldB(char *filename){
    // now we will create a file, and write a dataset into it.
    hid_t file_id, dset_id;
    hid_t file_space, memory_space;
    hid_t plist_id; //property list id
    plist_id = H5Pcreate(H5P_FILE_ACCESS); // access property list
    H5Pset_fapl_mpio(plist_id, mpi_row_comm, info);
    printf("[process id %d] trying to create a file to write fields\n",mpi_my_rank);
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id); //creating new file collectively
    printf("[process id %d] file created\n",mpi_my_rank);
    H5Pclose(plist_id);
    file_space = H5Screate_simple(hdf_rankFields, dataspace_dimsFields, NULL);
    memory_space = H5Screate_simple(hdf_rankFields, chunk_dimsFields, NULL);

    plist_id = H5Pcreate(H5P_DATASET_CREATE); //creating chunked dataset
    H5Pset_chunk(plist_id, hdf_rankFields, chunk_dimsFields);
    dset_id = H5Dcreate(file_id,"B", complex_id, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    H5Sclose(file_space);

    file_space = H5Dget_space(dset_id);
    status = H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offsetFields, strideFields, countFields, chunk_dimsFields);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    status = H5Dwrite(dset_id, complex_id, memory_space, file_space,plist_id, fields_fields.B);

    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Sclose(memory_space);
    H5Pclose(plist_id);
    H5Fclose(file_id);

};

/***************************
 *  hdf_saveFieldPhi
 * *************************/
void hdf_saveFieldPhi(char *filename){
    // now we will create a file, and write a dataset into it.
    hid_t file_id, dset_id;
    hid_t file_space, memory_space;
    hid_t plist_id; //property list id
    plist_id = H5Pcreate(H5P_FILE_ACCESS); // access property list
    H5Pset_fapl_mpio(plist_id, mpi_row_comm, info);
    printf("[process id %d] trying to create a file to write fields\n",mpi_my_rank);
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id); //creating new file collectively
    printf("[process id %d] file created\n",mpi_my_rank);
    H5Pclose(plist_id);
    file_space = H5Screate_simple(hdf_rankFields, dataspace_dimsFields, NULL);
    memory_space = H5Screate_simple(hdf_rankFields, chunk_dimsFields, NULL);

    plist_id = H5Pcreate(H5P_DATASET_CREATE); //creating chunked dataset
    H5Pset_chunk(plist_id, hdf_rankFields, chunk_dimsFields);
    dset_id = H5Dcreate(file_id,"Phi", complex_id, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    H5Sclose(file_space);

    file_space = H5Dget_space(dset_id);
    status = H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offsetFields, strideFields, countFields, chunk_dimsFields);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    status = H5Dwrite(dset_id, complex_id, memory_space, file_space,plist_id, fields_fields.phi);

    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Sclose(memory_space);
    H5Pclose(plist_id);
    H5Fclose(file_id);

};

/***************************
 *  hdf_saveEnergy
 * *************************/
void hdf_saveEnergy(int timestep)
{
    hid_t file_id, dset_id,dspace_id,group_id,filespace,memspace;
    hid_t plist_id; //property list id
    hid_t dims_ext[1] = {1};
    hid_t size[1];
    hid_t offset[1];
    plist_id = H5Pcreate(H5P_FILE_ACCESS); // access property list
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
    /* open file to read and write */
    file_id = H5Fopen(PARAMETER_FILENAME,H5F_ACC_RDWR,plist_id);
    H5Pclose(plist_id);
    /*opening a group and a dataset*/
    /*opening a group*/
    dset_id = H5Dopen2(file_id, "/freeEnergy/freeEnergy", H5P_DEFAULT);
    /*open a dataset*/
    dspace_id = H5Dget_space(dset_id);
    /*get dataset's dimensions */
    int ndims = H5Sget_simple_extent_ndims(dspace_id);
    hsize_t *dims = malloc(ndims * sizeof(*dims));
    H5Sget_simple_extent_dims(dspace_id, dims, NULL);
    H5Sclose(dspace_id);
    /*extend dataset size*/
    size[0] = dims[0] + dims_ext[0];
    offset[0] = dims[0];
    H5Dset_extent(dset_id, size);
    /*write free energy to the file*/
    dspace_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, NULL, dims_ext, NULL);
    memspace = H5Screate_simple(1,dims_ext,NULL);
    if(mpi_my_rank == 0)
    {
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace, dspace_id,plist_id,&diag_freeEnergy);
    }
    H5Sclose(dspace_id);
    H5Dclose(dset_id);
    /*write a timestep dataset*/
    /*opening a group*/
    dset_id = H5Dopen2(file_id, "/freeEnergy/timestep", H5P_DEFAULT);
    /*open a dataset...*/
    dspace_id = H5Dget_space(dset_id);
    /*... and extend it*/
    H5Dset_extent(dset_id, size);
    dspace_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, NULL, dims_ext, NULL);
    memspace = H5Screate_simple(1,dims_ext,NULL);
    if(mpi_my_rank == 0)
    {
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Dwrite(dset_id, H5T_NATIVE_INT,memspace, dspace_id,plist_id,&timestep);
    }
    H5Sclose(dspace_id);
    H5Dclose(dset_id);
    H5Fclose(file_id);
}

/***************************
 *  hdf_saveData
 * *************************/
void hdf_saveData(COMPLEX *h, int timestep) {
    if (parameters.save_EMfield && timestep % parameters.iter_EMfield == 0) {
        hdf_saveFields(timestep);
    }
    if (parameters.checkpoints && timestep % parameters.iter_checkpoint == 0)
    {
        hdf_createCheckpoint(h, timestep);
    }
    if (parameters.save_distrib && timestep % parameters.iter_distribution == 0)
    {
        hdf_saveDistrib(h,timestep);
    }
}

/***************************
 *
 *
 * PARAMETER FILE
 *
 *
 * *************************/

/***************************
 *  hdf_createParamFile
 * *************************/
void hdf_createParamFile()
{
    // now we will create a file, and write a dataset into it.
    hid_t file_id, dset_id,dspace_id,group_id,filespace;
    hid_t file_space, memory_space;
    hid_t plist_id; //property list id
    hid_t dims[1];
    hid_t dimsf[1];
    hid_t offsetKx[1];
    plist_id = H5Pcreate(H5P_FILE_ACCESS); // access property list
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);

    /* create a new file */
    file_id = H5Fcreate(PARAMETER_FILENAME,
                        H5F_ACC_TRUNC,
                        H5P_DEFAULT,
                        plist_id);
    H5Pclose(plist_id);
    /*create a group for storing kx,ky and kz arrays*/
    group_id = H5Gcreate2(file_id,
                          "/kSpace",
                          H5P_DEFAULT,
                          H5P_DEFAULT,
                          H5P_DEFAULT);
    /*create a datasets for storing kx,ky and kz*/
    /*create ky dataset*/
    dims[0] = array_local_size.nky;
    dspace_id = H5Screate_simple(1, dims, NULL);
    dset_id = H5Dcreate2(file_id, "/kSpace/ky", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT,
                            H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, space_ky);
    H5Sclose(dspace_id);
    H5Dclose(dset_id);
    /*create kz dataset*/
    dims[0] = array_local_size.nkz;
    dspace_id = H5Screate_simple(1, dims, NULL);
    dset_id = H5Dcreate2(file_id, "/kSpace/kz", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT,
                         H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, space_kz);
    H5Sclose(dspace_id);
    H5Dclose(dset_id);
    /*create kx dataset*/
    dims[0] = array_local_size.nkx;
    dimsf[0] = array_global_size.nkx;
    offsetKx[0] = offset[0];
    dspace_id = H5Screate_simple(1, dimsf, NULL);
    dset_id = H5Dcreate2(file_id, "/kSpace/kx", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT,
                         H5P_DEFAULT);
    H5Sclose(dspace_id);
    dspace_id = H5Screate_simple(1, dims, NULL);
    filespace = H5Dget_space(dset_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetKx, NULL, dims, NULL);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, dspace_id, filespace, plist_id, space_kx);


    H5Sclose(dspace_id);
    H5Sclose(filespace);
    H5Dclose(dset_id);
    /* Close the group. */
    H5Gclose(group_id);

    /*creating a group to save spectra*/
    if(parameters.save_diagnostics || parameters.compute_k)
    {
        group_id = H5Gcreate2(file_id, "/spectra", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Gclose(group_id);

        hid_t dims_spec_k[2] = {0,parameters.k_shells};
        hid_t dims_spec_m[2] = {0,parameters.nm};
        hid_t chunk_spec_k[2] = {1,parameters.k_shells};
        hid_t chunk_spec_m[2] = {1,parameters.nm};
        hid_t max_dims_k[2] = {H5S_UNLIMITED,parameters.k_shells};
        hid_t max_dims_m[2] = {H5S_UNLIMITED,parameters.nm};
        /*creating a dt and timestep datasets*/
        /* creating timestep dataset*/
        plist_id   = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 1, &chunk_spec_k[0]);
        dspace_id = H5Screate_simple(1,&dims_spec_k[0],&max_dims_k[0]);
        dset_id = H5Dcreate2(file_id, "/spectra/timestep", H5T_NATIVE_INT, dspace_id, H5P_DEFAULT, plist_id,
                             H5P_DEFAULT);
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);

        /* creating dt dataset*/
        plist_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 1, &chunk_spec_k[0]);
        dspace_id = H5Screate_simple(1,&dims_spec_k[0],&max_dims_k[0]);
        dset_id = H5Dcreate2(file_id, "/spectra/dt", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, plist_id,
                             H5P_DEFAULT);
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);
        if(parameters.compute_k)
        {
            /* creating shells dataset */
            plist_id = H5Pcreate(H5P_DATASET_CREATE);
            H5Pset_chunk(plist_id, 1, &chunk_spec_k[1]);
            dspace_id = H5Screate_simple(1,&chunk_spec_k[1],&max_dims_k[1]);
            dset_id = H5Dcreate2(file_id, "/spectra/shells", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, plist_id,
                                 H5P_DEFAULT);
            //printf("%f",diag_shells[2]);
            H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, diag_shells);
            H5Sclose(dspace_id);
            H5Dclose(dset_id);
            H5Pclose(plist_id);

            /*creating a kSpec dataset*/
            plist_id   = H5Pcreate(H5P_DATASET_CREATE);
            H5Pset_chunk(plist_id, 2, chunk_spec_k);
            dspace_id = H5Screate_simple(2,dims_spec_k,max_dims_k);
            dset_id = H5Dcreate2(file_id,
                                 "/spectra/kSpec",
                                 H5T_NATIVE_DOUBLE,
                                 dspace_id,
                                 H5P_DEFAULT,
                                 plist_id,
                                 H5P_DEFAULT);

            H5Sclose(dspace_id);
            H5Dclose(dset_id);
            H5Pclose(plist_id);
        }
        if(parameters.save_diagnostics)
        {
            /*creating an mSpec dataset*/
            plist_id   = H5Pcreate(H5P_DATASET_CREATE);
            H5Pset_chunk(plist_id, 2, chunk_spec_m);
            dspace_id = H5Screate_simple(2,dims_spec_m,max_dims_m);
            dset_id = H5Dcreate2(file_id, "/spectra/mSpec",
                                 H5T_NATIVE_DOUBLE,
                                 dspace_id,
                                 H5P_DEFAULT,
                                 plist_id,
                                 H5P_DEFAULT);
            H5Sclose(dspace_id);
            H5Dclose(dset_id);
            H5Pclose(plist_id);
        }

    }

    /*creating a group to save free energy*/
    if(parameters.save_diagnostics)
    {
        group_id = H5Gcreate2(file_id, "/freeEnergy", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Gclose(group_id);
        /* creating free energy dataset, dt dataset and timestep dataset*/
        /* creating free energy dataset */
        hid_t dims_energy[1] = {0};
        hid_t chunk_energy[1] = {1};
        plist_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 1, chunk_energy);
        hid_t maxdims[1] = {H5S_UNLIMITED};
        dspace_id = H5Screate_simple(1,dims_energy,maxdims);
        dset_id = H5Dcreate2(file_id, "/freeEnergy/freeEnergy", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, plist_id,
                             H5P_DEFAULT);
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);

        /*creating timestep dataset*/
        plist_id   = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 1, chunk_energy);
        dspace_id = H5Screate_simple(1,dims_energy,maxdims);
        dset_id = H5Dcreate2(file_id, "/freeEnergy/timestep", H5T_NATIVE_INT, dspace_id, H5P_DEFAULT, plist_id,
                             H5P_DEFAULT);
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);

        /*creating dt dataset*/
        plist_id  = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 1, chunk_energy);
        dspace_id = H5Screate_simple(1,dims_energy,maxdims);
        dset_id = H5Dcreate2(file_id, "/freeEnergy/dt", H5T_NATIVE_INT, dspace_id, H5P_DEFAULT, plist_id,
                             H5P_DEFAULT);
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);

    }
    /* Terminate access to the file. */
    H5Fclose(file_id);
}

/***************************
 *  hdf_createFiles
 * *************************/
void hdf_createFiles(){
    hdf_createParamFile();
    if (parameters.save_EMfield)
    {
        hdf_createFieldFile();
    }
}

/***************************
 *  hdf_saveKSpec
 * *************************/
void hdf_saveKSpec(int timestep) {
    hid_t file_id, dset_id,dspace_id,group_id,filespace,memspace;
    hid_t plist_id; //property list id
    hid_t dims_ext[2] = {1,parameters.k_shells};
    hid_t size[2];
    hid_t offset[2];
    plist_id = H5Pcreate(H5P_FILE_ACCESS); // access property list
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
    /* open file to read and write */
    file_id = H5Fopen(PARAMETER_FILENAME,H5F_ACC_RDWR,plist_id);
    H5Pclose(plist_id);
    /*opening a group and a dataset*/
    /*opening a group*/
    dset_id = H5Dopen2(file_id, "/spectra/kSpec", H5P_DEFAULT);
    /*open a dataset*/
    dspace_id = H5Dget_space(dset_id);
    /*get dataset's dimensions */
    int ndims = H5Sget_simple_extent_ndims(dspace_id);
    hsize_t *dims = malloc(ndims * sizeof(*dims));
    H5Sget_simple_extent_dims(dspace_id, dims, NULL);
    H5Sclose(dspace_id);
    /*extend dataset size*/
    size[0] = dims[0] + dims_ext[0];
    size[1] = dims[1];
    offset[0] = dims[0];
    offset[1] = 0;
    H5Dset_extent(dset_id, size);
    /*write free energy to the file*/
    dspace_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, NULL, dims_ext, NULL);
    memspace = H5Screate_simple(2,dims_ext,NULL);

    if(mpi_my_rank == 0)
    {
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace, dspace_id,plist_id,diag_kSpec);
    }
    H5Sclose(dspace_id);
    H5Dclose(dset_id);
    /*write a timestep dataset*/
    /*opening a group*/
    dset_id = H5Dopen2(file_id, "/spectra/timestep", H5P_DEFAULT);
    /*open a dataset...*/
    dspace_id = H5Dget_space(dset_id);
    /*... and extend it*/
    H5Dset_extent(dset_id, &size[0]);
    dspace_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, &offset[0], NULL, &dims_ext[0], NULL);
    memspace = H5Screate_simple(1,&dims_ext[0],NULL);
    if(mpi_my_rank == 0)
    {
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Dwrite(dset_id, H5T_NATIVE_INT,memspace, dspace_id,plist_id,&timestep);
    }
    H5Sclose(dspace_id);
    H5Dclose(dset_id);
    H5Fclose(file_id);
};

/***************************
 *  hdf_saveMSpec
 * *************************/
void hdf_saveMSpec(int timestep){
    hid_t file_id, dset_id,dspace_id,group_id,filespace,memspace;
    hid_t plist_id; //property list id
    hid_t dims_ext_full[2] = {1,parameters.nm};
    hid_t dims_ext_local[2] = {1,array_local_size.nm};
    hid_t max_dims_local[2] = {H5S_UNLIMITED,array_local_size.nm};
    hid_t max_dims_full[2] = {H5S_UNLIMITED,parameters.nm};
    hid_t size[2];
    hid_t offset[2];
    hid_t stride[2] = {1,1};
    hid_t count[2] = {1,1};
    plist_id = H5Pcreate(H5P_FILE_ACCESS); // access property list
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
    /* open file to read and write */
    file_id = H5Fopen(PARAMETER_FILENAME,H5F_ACC_RDWR,plist_id);
    H5Pclose(plist_id);
    /*opening a group and a dataset*/
    /*opening a group*/
    dset_id = H5Dopen2(file_id, "/spectra/mSpec", H5P_DEFAULT);
    /*open a dataset*/
    dspace_id = H5Dget_space(dset_id);
    /*get dataset's dimensions */
    int ndims = H5Sget_simple_extent_ndims(dspace_id);
    hsize_t *dims = malloc(ndims * sizeof(*dims));
    H5Sget_simple_extent_dims(dspace_id, dims, NULL);
    H5Sclose(dspace_id);
    /*extend dataset size*/
    size[0] = dims[0] + dims_ext_full[0];
    size[1] = dims[1];
    offset[0] = dims[0];
    offset[1] =  mpi_my_col_rank * parameters.nm / mpi_dims[0];
    H5Dset_extent(dset_id, size);
    /*write m spec to the file*/
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    dspace_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, stride,count, dims_ext_local);
    memspace = H5Screate_simple(1,&dims_ext_local[1],&max_dims_local[1]);
    if (mpi_my_row_rank != 0)
    {
        H5Sselect_none(dspace_id);
        H5Sselect_none(memspace);
    }
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, dspace_id, plist_id, diag_mSpec);
    H5Sclose(memspace);
    H5Sclose(dspace_id);
    H5Pclose(plist_id);
    H5Dclose(dset_id);
    H5Fclose(file_id);
}

/***************************
 *
 *
 * CHECKPOINTS
 *
 *
 * *************************/

/***************************
 *  hdf_initCheckpoints
 * *************************/
void hdf_initCheckpoints(){
    if (VERBOSE) printf("initialising checkpoints\n");
    hdf_checkpointNames = malloc(parameters.checkpoints * sizeof(hdf_checkpointNames));
    for (size_t i = 0; i < parameters.checkpoints; i++){
        hdf_checkpointNames[i] = malloc(FILENAME_ID_LEN * sizeof(*hdf_checkpointNames[i]));
    }
};

/***************************
 *  hdf_createCheckpoint
 * *************************/
void hdf_createCheckpoint(COMPLEX *h, int timestep) {
    if (VERBOSE) printf("create checkpoint\n");
    sprintf(hdf_newCheckpointName, "%s%s_%d.h5", CHECKPOINT_PATH,CHCK_NAME,timestep);
    if (VERBOSE) printf("checkpoint name is %s\n",hdf_newCheckpointName);
    if (hdf_checkpointCount > parameters.checkpoints - 1)
    {
        if (VERBOSE) printf("condition checkpoint count > checkpoint reached!\n");
        if (mpi_my_rank == CHECKPOINT_ROOT)
        {
            if (access(hdf_checkpointNames[0],F_OK) == 0)
            {
                if(VERBOSE) printf("file %s exists", hdf_checkpointNames[0]);
                remove(hdf_checkpointNames[0]);
            }
        }
        for(size_t i = 0; i < parameters.checkpoints - 1; i++){
            strcpy(hdf_checkpointNames[i],hdf_checkpointNames[i + 1]);
        }
        strcpy(hdf_checkpointNames[parameters.checkpoints - 1], hdf_newCheckpointName);
    }
    else
    {
        strcpy(hdf_checkpointNames[hdf_checkpointCount], hdf_newCheckpointName);

    }
    hdf_dumpCheckpoint(h, timestep, hdf_newCheckpointName);
    hdf_checkpointCount++;
};

/***************************
 *  hdf_dumpCheckpoint
 * *************************/
void hdf_dumpCheckpoint(COMPLEX *h, int timestep, char *filename){
    // now we will create a file, and write a dataset into it.
    hid_t file_id, dset_id;
    hid_t file_space, memory_space;
    hid_t plist_id; //property list id
    hid_t dims_timestep[1] = {1};
    hid_t rank_timestep = 1;
    plist_id = H5Pcreate(H5P_FILE_ACCESS); // access property list
    H5Pset_fapl_mpio(plist_id, mpi_cube_comm, info);
    /* creating the file */
    if (VERBOSE) printf("[process id %d] trying to create a file\n",mpi_my_rank);
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id); //creating new file collectively
    if (VERBOSE) printf("[process id %d] file created\n",mpi_my_rank);
    H5Pclose(plist_id);
    /* create file space and memory space. Memory space holds each process' data, file space going to hold whole data */
    file_space = H5Screate_simple(hdf_rank, dataspace_dims_c, NULL);
    memory_space = H5Screate_simple(hdf_rank, chunk_dims_c, NULL);

    /*
     * creating h dataset
     */
    plist_id = H5Pcreate(H5P_DATASET_CREATE); //creating chunked dataset
    H5Pset_chunk(plist_id, hdf_rank, chunk_dims_c);
    dset_id = H5Dcreate(file_id,"h", complex_id, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    H5Sclose(file_space);
    /* writing dataset collectively */
    file_space = H5Dget_space(dset_id);
    status = H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset, stride, count, chunk_dims_c);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    status = H5Dwrite(dset_id, complex_id, memory_space, file_space,plist_id, h);
    /*closing dataset identifier, memory space and file space*/
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Sclose(memory_space);
    H5Pclose(plist_id);

    /* writing timestep information */
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_timestep, dims_timestep, NULL);
    dset_id = H5Dcreate(file_id,"timestep", H5T_NATIVE_INT, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_timestep, dims_timestep, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id,H5T_NATIVE_INT, memory_space, file_space, plist_id, &timestep);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);

    /* writing beta information */
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_timestep, dims_timestep, NULL);
    dset_id = H5Dcreate(file_id,"beta", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_timestep, dims_timestep, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id,H5T_NATIVE_DOUBLE, memory_space, file_space, plist_id, &parameters.beta);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);

    /* writing nkx information */
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_timestep, dims_timestep, NULL);
    dset_id = H5Dcreate(file_id,"nkx", H5T_NATIVE_INT, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_timestep, dims_timestep, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id,H5T_NATIVE_INT, memory_space, file_space, plist_id, &parameters.nkx);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);

    /* writing nky information */
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_timestep, dims_timestep, NULL);
    dset_id = H5Dcreate(file_id,"nky", H5T_NATIVE_INT, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_timestep, dims_timestep, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id,H5T_NATIVE_INT, memory_space, file_space, plist_id, &parameters.nky);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);

    /* writing nkz information */
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_timestep, dims_timestep, NULL);
    dset_id = H5Dcreate(file_id,"nkz", H5T_NATIVE_INT, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_timestep, dims_timestep, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id,H5T_NATIVE_INT, memory_space, file_space, plist_id, &parameters.nkz);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);

    /* writing nm information */
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_timestep, dims_timestep, NULL);
    dset_id = H5Dcreate(file_id,"nm", H5T_NATIVE_INT, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_timestep, dims_timestep, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id,H5T_NATIVE_INT, memory_space, file_space, plist_id, &parameters.nm);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);

    /* writing nl information */
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_timestep, dims_timestep, NULL);
    dset_id = H5Dcreate(file_id,"nl", H5T_NATIVE_INT, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_timestep, dims_timestep, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id,H5T_NATIVE_INT, memory_space, file_space, plist_id, &parameters.nl);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);

    /* writing ns information */
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_timestep, dims_timestep, NULL);
    dset_id = H5Dcreate(file_id,"ns", H5T_NATIVE_INT, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_timestep, dims_timestep, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id,H5T_NATIVE_INT, memory_space, file_space, plist_id, &parameters.ns);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);

    /*particle data write*/
    int rank_particle = 1;
    hsize_t *dims_particle = malloc( sizeof(*dims_particle));
    dims_particle[0] = parameters.ns;
    /* writing charge information */
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_particle, dims_particle, NULL);
    dset_id = H5Dcreate(file_id,"charge", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_particle, dims_particle, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id,H5T_NATIVE_DOUBLE, memory_space, file_space, plist_id, parameters.charge);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);
    /* writing temperature information */
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_particle, dims_particle, NULL);
    dset_id = H5Dcreate(file_id,"temperature", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_particle, dims_particle, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id,H5T_NATIVE_DOUBLE, memory_space, file_space, plist_id, parameters.temperature);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);
    /* writing density information */
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_particle, dims_particle, NULL);
    dset_id = H5Dcreate(file_id,"density", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_particle, dims_particle, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id,H5T_NATIVE_DOUBLE, memory_space, file_space, plist_id, parameters.density);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);
    /* writing mass information */
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_particle, dims_particle, NULL);
    dset_id = H5Dcreate(file_id,"mass", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_particle, dims_particle, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id,H5T_NATIVE_DOUBLE, memory_space, file_space, plist_id, parameters.mass);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);

    H5Fclose(file_id);
}

void hdf_saveDistrib(COMPLEX* h, int timestep){
    char distrib_name[FILENAME_ID_LEN];
    sprintf(distrib_name, "%s%d.h5", DISTRIBUTION_FILENAME, timestep);
    hdf_dumpCheckpoint(h,  timestep, distrib_name);
};

/***************************
 *
 *
 * FIELD FILE
 *
 *
 * *************************/

/***************************
 *  hdf_createFieldFile
 * *************************/
void hdf_createFieldFile(){
    // now we will create a file, and write a dataset into it.
    hid_t file_id, dset_id,dspace_id,filespace;
    hid_t file_space, memory_space;
    hid_t plist_id; //property list id
    hid_t dims[1];
    hid_t dimsf[1];
    hid_t offsetKx[1];

    /* variables needed to write fields to the file  */
    hid_t dims_field[4] = {0,array_global_size.nkx, array_global_size.nky, array_global_size.nkz};
    hid_t chunk_field[4] = {1,array_global_size.nkx, array_global_size.nky, array_global_size.nkz};
    hid_t max_dims_field[4] = {H5S_UNLIMITED,array_global_size.nkx, array_global_size.nky, array_global_size.nkz};
    plist_id = H5Pcreate(H5P_FILE_ACCESS); // access property list
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);

    /* create a new file */
    file_id = H5Fcreate(FIELD_FILENAME,
                        H5F_ACC_TRUNC,
                        H5P_DEFAULT,
                        plist_id);
    H5Pclose(plist_id);

    /*create a datasets for storing kx,ky and kz*/
    /*create ky dataset*/
    dims[0] = array_local_size.nky;
    dspace_id = H5Screate_simple(1, dims, NULL);
    dset_id = H5Dcreate2(file_id, "ky", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT,
                         H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, space_ky);
    H5Sclose(dspace_id);
    H5Dclose(dset_id);

    /*create kz dataset*/
    dims[0] = array_local_size.nkz;
    dspace_id = H5Screate_simple(1, dims, NULL);
    dset_id = H5Dcreate2(file_id, "kz", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT,
                         H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, space_kz);
    H5Sclose(dspace_id);
    H5Dclose(dset_id);

    /*create kx dataset*/
    dims[0] = array_local_size.nkx;
    dimsf[0] = array_global_size.nkx;
    offsetKx[0] = offset[0];
    dspace_id = H5Screate_simple(1, dimsf, NULL);
    dset_id = H5Dcreate2(file_id, "kx", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT,
                         H5P_DEFAULT);
    H5Sclose(dspace_id);
    dspace_id = H5Screate_simple(1, dims, NULL);
    filespace = H5Dget_space(dset_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offsetKx, NULL, dims, NULL);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, dspace_id, filespace, plist_id, space_kx);
    H5Sclose(dspace_id);
    H5Dclose(dset_id);

    /*creating timestep dataset*/
    plist_id   = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id, 1, &chunk_field[0]);
    dspace_id = H5Screate_simple(1,&dims_field[0],&max_dims_field[0]);
    dset_id = H5Dcreate2(file_id, "timestep", H5T_NATIVE_INT, dspace_id, H5P_DEFAULT, plist_id,
                         H5P_DEFAULT);
    H5Sclose(dspace_id);
    H5Dclose(dset_id);
    H5Pclose(plist_id);

    /* creating fields dataset */
    if (VERBOSE) printf("systemType = %d\n",systemType);
    if (systemType == ELECTROSTATIC)
    {

        /* create phi dataset */
        plist_id   = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 4, chunk_field);
        dspace_id = H5Screate_simple(4,dims_field,max_dims_field);
        dset_id = H5Dcreate2(file_id,
                             "phi",
                             complex_id,
                             dspace_id,
                             H5P_DEFAULT,
                             plist_id,
                             H5P_DEFAULT);

        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);

    }
    if (systemType == ELECTROMAGNETIC)
    {
        /* create phi dataset */
        plist_id   = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 4, chunk_field);
        dspace_id = H5Screate_simple(4,dims_field,max_dims_field);
        dset_id = H5Dcreate2(file_id,
                             "phi",
                             complex_id,
                             dspace_id,
                             H5P_DEFAULT,
                             plist_id,
                             H5P_DEFAULT);

        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);

        /* create A dataset */
        plist_id   = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 4, chunk_field);
        dspace_id = H5Screate_simple(4,dims_field,max_dims_field);
        dset_id = H5Dcreate2(file_id,
                             "A",
                             complex_id,
                             dspace_id,
                             H5P_DEFAULT,
                             plist_id,
                             H5P_DEFAULT);

        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);

        /* create B dataset */
        plist_id   = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 4, chunk_field);
        dspace_id = H5Screate_simple(4,dims_field,max_dims_field);
        dset_id = H5Dcreate2(file_id,
                             "B",
                             complex_id,
                             dspace_id,
                             H5P_DEFAULT,
                             plist_id,
                             H5P_DEFAULT);

        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);
    }
    /*create phi dataset*/
   // H5Sclose(dspace_id);
    //H5Sclose(filespace);
    //H5Dclose(dset_id);
    H5Fclose(file_id);
};

/***************************
 *  hdf_saveFields
 * *************************/
void hdf_saveFields(int timestep){
    hid_t file_id, dset_id,dspace_id,group_id,filespace,memspace;
    hid_t plist_id; //property list id
    hid_t dims_ext_full[4] = {1,array_global_size.nkx, array_global_size.nky, array_global_size.nkz};
    hid_t dims_ext_local[4] = {1,array_local_size.nkx, array_local_size.nky, array_local_size.nkz};
    hid_t max_dims_local[4] = {H5S_UNLIMITED,array_local_size.nkx, array_local_size.nky, array_local_size.nkz};
    hid_t max_dims_full[4] = {H5S_UNLIMITED,array_local_size.nkx, array_local_size.nky, array_local_size.nkz};
    hid_t size[4],size_timestep[1];
    hid_t offset[4],offset_timestep[1];
    hid_t stride[4] = {1,1,1,1};
    hid_t count[4] = {1,1,1,1};

    plist_id = H5Pcreate(H5P_FILE_ACCESS); // access property list
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
    /* open file to read and write */
    file_id = H5Fopen(FIELD_FILENAME,H5F_ACC_RDWR,plist_id);
    H5Pclose(plist_id);

    if (systemType == ELECTROSTATIC)
    {
        /*writing phi dataset*/
        /*opening a group and a dataset*/
        /*opening a group*/
        dset_id = H5Dopen2(file_id, "phi", H5P_DEFAULT);
        /*open a dataset*/
        dspace_id = H5Dget_space(dset_id);
        /*get dataset's dimensions */
        int ndims = H5Sget_simple_extent_ndims(dspace_id);
        hsize_t *dims = malloc(ndims * sizeof(*dims));
        H5Sget_simple_extent_dims(dspace_id, dims, NULL);
        H5Sclose(dspace_id);
        /*extend dataset size*/
        size[0] = dims[0] + dims_ext_full[0];
        size[1] = dims[1];
        size[2] = dims[2];
        size[3] = dims[3];
        offset[0] = dims[0];
        offset[1] =  mpi_my_row_rank * parameters.nkx / mpi_dims[1];
        offset[2] =  0;
        offset[3] =  0;
        H5Dset_extent(dset_id, size);
        /*write m spec to the file*/
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        dspace_id = H5Dget_space(dset_id);
        H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, stride,count, dims_ext_local);
        memspace = H5Screate_simple(3,&dims_ext_local[1],&max_dims_local[1]);
        if (mpi_my_col_rank != 0)
        {
            H5Sselect_none(dspace_id);
            H5Sselect_none(memspace);
        }
        H5Dwrite(dset_id, complex_id, memspace, dspace_id, plist_id, fields_fields.phi);
        H5Sclose(memspace);
        H5Sclose(dspace_id);
        H5Pclose(plist_id);
        H5Dclose(dset_id);
        //H5Fclose(file_id);
    }
    if (systemType == ELECTROMAGNETIC)
    {
        /* writing phi dataset */
        /*open a dataset*/
        dset_id = H5Dopen2(file_id, "phi", H5P_DEFAULT);
        dspace_id = H5Dget_space(dset_id);
        /*get dataset's dimensions */
        int ndims = H5Sget_simple_extent_ndims(dspace_id);
        hsize_t *dims = malloc(ndims * sizeof(*dims));
        H5Sget_simple_extent_dims(dspace_id, dims, NULL);
        H5Sclose(dspace_id);
        /*extend dataset size*/
        size[0] = dims[0] + dims_ext_full[0];
        size[1] = dims[1];
        size[2] = dims[2];
        size[3] = dims[3];
        offset[0] = dims[0];
        offset[1] =  mpi_my_row_rank * parameters.nkx / mpi_dims[1];
        offset[2] =  0;
        offset[3] =  0;
        H5Dset_extent(dset_id, size);
        /*write phi to the file*/
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        dspace_id = H5Dget_space(dset_id);
        H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, stride,count, dims_ext_local);
        memspace = H5Screate_simple(3,&dims_ext_local[1],&max_dims_local[1]);
        if (mpi_my_col_rank != 0)
        {
            H5Sselect_none(dspace_id);
            H5Sselect_none(memspace);
        }
        H5Dwrite(dset_id, complex_id, memspace, dspace_id, plist_id, fields_fields.phi);
        H5Sclose(memspace);
        H5Sclose(dspace_id);
        H5Pclose(plist_id);
        H5Dclose(dset_id);
        /* writing B dataset */
        dset_id = H5Dopen2(file_id, "B", H5P_DEFAULT);
        dspace_id = H5Dget_space(dset_id);
        H5Dset_extent(dset_id, size);
        /*write B to the file*/
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        dspace_id = H5Dget_space(dset_id);
        H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, stride,count, dims_ext_local);
        memspace = H5Screate_simple(3,&dims_ext_local[1],&max_dims_local[1]);
        if (mpi_my_col_rank != 0)
        {
            H5Sselect_none(dspace_id);
            H5Sselect_none(memspace);
        }
        H5Dwrite(dset_id, complex_id, memspace, dspace_id, plist_id, fields_fields.B);
        H5Sclose(memspace);
        H5Sclose(dspace_id);
        H5Pclose(plist_id);
        H5Dclose(dset_id);
        /* writing A dataset */
        dset_id = H5Dopen2(file_id, "A", H5P_DEFAULT);
        dspace_id = H5Dget_space(dset_id);
        H5Dset_extent(dset_id, size);
        /*write A to the file*/
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        dspace_id = H5Dget_space(dset_id);
        H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, stride,count, dims_ext_local);
        memspace = H5Screate_simple(3,&dims_ext_local[1],&max_dims_local[1]);
        if (mpi_my_col_rank != 0)
        {
            H5Sselect_none(dspace_id);
            H5Sselect_none(memspace);
        }
        H5Dwrite(dset_id, complex_id, memspace, dspace_id, plist_id, fields_fields.A);
        H5Sclose(memspace);
        H5Sclose(dspace_id);
        H5Pclose(plist_id);
        H5Dclose(dset_id);
        //H5Fclose(file_id);
    }
    /*write a timestep dataset*/
    /*opening a group*/
    dset_id = H5Dopen2(file_id, "timestep", H5P_DEFAULT);
    /*open a dataset...*/
    dspace_id = H5Dget_space(dset_id);
    /*get dataset's dimensions */
    int ndims = H5Sget_simple_extent_ndims(dspace_id);
    hsize_t *dims = malloc(ndims * sizeof(*dims));
    H5Sget_simple_extent_dims(dspace_id, dims, NULL);
    H5Sclose(dspace_id);
    /*... and extend it*/
    size_timestep[0] = dims[0] + 1;
    H5Dset_extent(dset_id, size_timestep);
    dspace_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, &offset[0], NULL, &dims_ext_full[0], NULL);
    memspace = H5Screate_simple(1,&dims_ext_full[0],NULL);
    if(mpi_my_rank == 0)
    {
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Dwrite(dset_id, H5T_NATIVE_INT,memspace, dspace_id,plist_id,&timestep);
    }
    H5Sclose(dspace_id);
    H5Dclose(dset_id);
    H5Fclose(file_id);
};

/***************************
 *
 *
 * READ FILE
 *
 *
 * *************************/

/***************************
 *  hdf_readData
 * *************************/
void hdf_readData(char *filename, COMPLEX *h) {
    // now we will create a file, and write a dataset into it.
    hid_t file_id, dset_id;
    hid_t file_space, memory_space;
    hid_t plist_id; //property list id
    hid_t dims_timestep[1] = {1};
    hid_t rank_timestep = 1;
    plist_id = H5Pcreate(H5P_FILE_ACCESS); // access property list
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
    /* open file to read */
    file_id = H5Fopen(filename,H5F_ACC_RDONLY,plist_id);
    H5Pclose(plist_id);
    /*open dataset*/
    dset_id = H5Dopen2(file_id, "h", H5P_DEFAULT);
    file_space = H5Dget_space(dset_id);
    H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset, stride, count, chunk_dims_c);
    memory_space = H5Screate_simple(hdf_rank, chunk_dims_c, NULL);
    /*
     * reading h dataset
     */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dread(dset_id,complex_id, memory_space, file_space, plist_id, h);
    /*closing dataset identifier, memory space and file space*/
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Sclose(memory_space);
    H5Pclose(plist_id);

};