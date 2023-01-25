/**************************************
* @file utils_hdf.c
* \brief hdf module
*
* contains HDF related routines to save and read hdf files
***************************************/
////////////////////////////////////////////////////////////////////////////////
// 16/12/2021 created by Gene Gorbunov
//                                   HDF5 UTILITIES
//
// hdf_init
// hdf_create_file_c
// hdf_create_file_r
// hdf_initField
// hdf_saveFieldA
// hdf_saveField_r
// hdf_saveFieldB
// hdf_saveFieldPhi
// hdf_saveData
// hdf_createFiles
// hdf_saveEnergy
// hdf_saveKSpec
// hdf_saveMSpec
// hdf_createFieldFile
// hdf_saveFields
// hdf_createCheckpoint
// hdf_initCheckpoints
// hdf_dumpCheckpoint
// hdf_dumpCheckpointReal
// hdf_readData
// hdf_saveDistrib
// hdf_createSaveDirs
//
// VERSION 1.0
////////////////////////////////////////////////////////////////////////////////

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
#define VERBOSE 0
#define FILENAME_ID_LEN 128
#define CHECKPOINT_ROOT 0
#define PATH_LEN 128

int hdf_rank = 6;
int hdf_rankFields = 3;
int hdf_rankChi = 5;
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
hsize_t dataspace_dimsFields_r[3];
hsize_t dataspace_dimsChi[5];
hsize_t dataspace_dimsChi_r[5];
hsize_t chunk_dims_r[6];
hsize_t chunk_dims_c[6];
hsize_t chunk_dimsFields[3];
hsize_t chunk_dimsFields_r[3];
hsize_t chunk_dimsChi[5];
hsize_t chunk_dimsChi_r[5];
hsize_t offset[6];
hsize_t offsetFields[3];
hsize_t offsetFields_r[3];
hsize_t offsetChi[5];
hsize_t offsetChi_r[5];
hsize_t count[6] = {1,1,1,1,1,1};
hsize_t stride[6] = {1,1,1,1,1,1};
hsize_t countFields[3] = {1,1,1};
hsize_t strideFields[3] = {1,1,1};
hsize_t countChi[5] = {1,1,1,1,1};
hsize_t strideChi[5] = {1,1,1,1,1};
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
    hdf_initChi();
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
    chunk_dims_r[2] = array_local_size.nz + 2;
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

    if(VERBOSE) printf("[MPI process %d] my coords are (%d, %d), offsets are (%lld,%lld,%lld,%lld,%lld,%lld)\n",mpi_my_rank,
           mpi_my_coords[0],
           mpi_my_coords[1],
           offset[0],
           offset[1],
           offset[2],
           offset[3],
           offset[4],
           offset[5]);
    if(VERBOSE) printf("[MPI process %d] my coords are (%d, %d), chunk dims for real are (%lld,%lld,%lld,%lld,%lld,%lld)\n",mpi_my_rank,
           mpi_my_coords[0],
           mpi_my_coords[1],
           chunk_dims_r[0],
           chunk_dims_r[1],
           chunk_dims_r[2],
           chunk_dims_r[3],
           chunk_dims_r[4],
           chunk_dims_r[5]);
    if(VERBOSE) printf("[MPI process %d] my coords are (%d, %d), chunk dims for complex are (%lld,%lld,%lld,%lld,%lld,%lld)\n",mpi_my_rank,
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
    if(VERBOSE) printf("[process id %d] trying to create a file\n",mpi_my_rank);
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id); //creating new file collectively
    if(VERBOSE) printf("[process id %d] file created\n",mpi_my_rank);
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
    //transpose the data
    fftw_copy_buffer_r(fftw_hBuf,data);
    fftw_transposeToXY();
    // now we will create a file, and write a dataset into it.
    hid_t file_id, dset_id;
    hid_t file_space, memory_space;
    hid_t plist_id; //property list id
    plist_id = H5Pcreate(H5P_FILE_ACCESS); // access property list
    H5Pset_fapl_mpio(plist_id, mpi_cube_comm, info);
    if(VERBOSE) printf("[process id %d] trying to create a file\n",mpi_my_rank);
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id); //creating new file collectively
    if(VERBOSE) printf("[process id %d] file created\n",mpi_my_rank);
    H5Pclose(plist_id);
    file_space = H5Screate_simple(hdf_rank, dataspace_dims_r, NULL);
    memory_space = H5Screate_simple(hdf_rank, chunk_dims_r, NULL);

    plist_id = H5Pcreate(H5P_DATASET_CREATE); //creating chunked dataset
    H5Pset_chunk(plist_id, hdf_rank, chunk_dims_r);
    dset_id = H5Dcreate(file_id,"hr", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    H5Sclose(file_space);

    file_space = H5Dget_space(dset_id);
    status = H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset, stride, count, chunk_dims_r);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memory_space, file_space,plist_id, (double*)fftw_hBuf);

    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Sclose(memory_space);
    H5Pclose(plist_id);
    H5Fclose(file_id);

    //transpose data back (in general, not required, added to make function foolproof)
    fftw_transposeToYX();
}

/***************************
 *  hdf_initChi
 * *************************/
void hdf_initChi(){
    if (systemType == ELECTROSTATIC){
        dataspace_dimsChi[0] = array_global_size.nkx;
        dataspace_dimsChi[1] = array_global_size.nky;
        dataspace_dimsChi[2] = array_global_size.nkz;
        dataspace_dimsChi[3] = array_global_size.ns;
        dataspace_dimsChi[4] = 1;

        dataspace_dimsChi_r[0] = array_global_size.nkx;
        dataspace_dimsChi_r[1] = array_global_size.nky;
        dataspace_dimsChi_r[2] = array_global_size.nz+2;
        dataspace_dimsChi_r[3] = array_global_size.ns;
        dataspace_dimsChi_r[4] = 1;

        chunk_dimsChi[0] = array_local_size.nkx;
        chunk_dimsChi[1] = array_local_size.nky;
        chunk_dimsChi[2] = array_local_size.nkz;
        chunk_dimsChi[3] = array_local_size.ns;
        chunk_dimsChi[4] = 1;
        if(VERBOSE) printf("[MPI process %d] %llu %llu %llu\n",mpi_my_rank, chunk_dimsFields[0],chunk_dimsFields[1],chunk_dimsFields[2]);

        chunk_dimsChi_r[0] = array_local_size.nkx;
        chunk_dimsChi_r[1] = array_local_size.nky;
        chunk_dimsChi_r[2] = array_local_size.nz + 2;
        chunk_dimsChi_r[3] = array_local_size.ns;
        chunk_dimsChi_r[4] = 1;

        offsetChi[0] = array_local_size.nkx * mpi_my_coords[1];
        offsetChi[1] = 0;
        offsetChi[2] = 0;
        offsetChi[3] = 0;
        offsetChi[4] = 0;

        offsetChi_r[0] = array_local_size.nkx * mpi_my_coords[1];
        offsetChi_r[1] = 0;
        offsetChi_r[2] = 0;
        offsetChi_r[3] = 0;
        offsetChi_r[4] = 0;

    }
    else{
        dataspace_dimsChi[0] = array_global_size.nkx;
        dataspace_dimsChi[1] = array_global_size.nky;
        dataspace_dimsChi[2] = array_global_size.nkz;
        dataspace_dimsChi[3] = array_global_size.ns;
        dataspace_dimsChi[4] = 3;

        dataspace_dimsChi_r[0] = array_global_size.nkx;
        dataspace_dimsChi_r[1] = array_global_size.nky;
        dataspace_dimsChi_r[2] = array_global_size.nz+2;
        dataspace_dimsChi_r[3] = array_global_size.ns;
        dataspace_dimsChi_r[4] = 3;

        chunk_dimsChi[0] = array_local_size.nkx;
        chunk_dimsChi[1] = array_local_size.nky;
        chunk_dimsChi[2] = array_local_size.nkz;
        chunk_dimsChi[3] = array_local_size.ns;
        chunk_dimsChi[4] = 3;
        if(VERBOSE) printf("[MPI process %d] %llu %llu %llu\n",mpi_my_rank, chunk_dimsFields[0],chunk_dimsFields[1],chunk_dimsFields[2]);

        chunk_dimsChi_r[0] = array_local_size.nkx;
        chunk_dimsChi_r[1] = array_local_size.nky;
        chunk_dimsChi_r[2] = array_local_size.nz + 2;
        chunk_dimsChi_r[3] = array_local_size.ns;
        chunk_dimsChi_r[4] = 3;

        offsetChi[0] = array_local_size.nkx * mpi_my_coords[1];
        offsetChi[1] = 0;
        offsetChi[2] = 0;
        offsetChi[3] = 0;
        offsetChi[4] = 0;

        offsetChi_r[0] = array_local_size.nkx * mpi_my_coords[1];
        offsetChi_r[1] = 0;
        offsetChi_r[2] = 0;
        offsetChi_r[3] = 0;
        offsetChi_r[4] = 0;

    }
};

/***************************
 *  hdf_createChiFile_r
 * *************************/
void hdf_createChiFile_r(char *filename, double *data){
    //transpose data before saving
    fftw_copyChiBuf_r(fftw_chiBuf,data);
    fftw_transposeToXY_chi();
    // now we will create a file, and write a dataset into it.
    hid_t file_id, dset_id;
    hid_t file_space, memory_space;
    hid_t plist_id; //property list id
    plist_id = H5Pcreate(H5P_FILE_ACCESS); // access property list
    H5Pset_fapl_mpio(plist_id, mpi_cube_comm, info);
    if(VERBOSE) printf("[process id %d] trying to create a file\n",mpi_my_rank);
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id); //creating new file collectively
    if(VERBOSE) printf("[process id %d] file created\n",mpi_my_rank);
    H5Pclose(plist_id);
    file_space = H5Screate_simple(hdf_rankChi, dataspace_dimsChi_r, NULL);
    memory_space = H5Screate_simple(hdf_rankChi, chunk_dimsChi_r, NULL);

    plist_id = H5Pcreate(H5P_DATASET_CREATE); //creating chunked dataset
    H5Pset_chunk(plist_id, hdf_rankChi, chunk_dimsChi_r);
    dset_id = H5Dcreate(file_id,"chi", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    H5Sclose(file_space);

    file_space = H5Dget_space(dset_id);
    status = H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offsetChi_r, strideChi, countChi, chunk_dimsChi_r);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memory_space, file_space,plist_id, fftw_chiBuf);

    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Sclose(memory_space);
    H5Pclose(plist_id);
    H5Fclose(file_id);

    // transposing data to original shape
    fftw_transposeToYX_chi();

}

/***************************
 *  hdf_createChiFile_c
 * *************************/
void hdf_createChiFile_c(char *filename, COMPLEX *data){
    // now we will create a file, and write a dataset into it.
    hid_t file_id, dset_id;
    hid_t file_space, memory_space;
    hid_t plist_id; //property list id
    plist_id = H5Pcreate(H5P_FILE_ACCESS); // access property list
    H5Pset_fapl_mpio(plist_id, mpi_cube_comm, info);
    if(VERBOSE) printf("[process id %d] trying to create a file\n",mpi_my_rank);
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id); //creating new file collectively
    if(VERBOSE) printf("[process id %d] file created\n",mpi_my_rank);
    H5Pclose(plist_id);
    file_space = H5Screate_simple(hdf_rank, dataspace_dimsChi, NULL);
    memory_space = H5Screate_simple(hdf_rank, chunk_dimsChi, NULL);

    plist_id = H5Pcreate(H5P_DATASET_CREATE); //creating chunked dataset
    H5Pset_chunk(plist_id, hdf_rankChi, chunk_dimsChi);
    dset_id = H5Dcreate(file_id,"chi", complex_id, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    H5Sclose(file_space );

    file_space = H5Dget_space(dset_id);
    status = H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offsetChi, strideChi, countChi, chunk_dimsChi);
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
 *  hdf_initField
 * *************************/
void hdf_initField(){
    dataspace_dimsFields[0] = array_global_size.nkx;
    dataspace_dimsFields[1] = array_global_size.nky;
    dataspace_dimsFields[2] = array_global_size.nkz;

    dataspace_dimsFields_r[0] = array_global_size.nkx;
    dataspace_dimsFields_r[1] = array_global_size.nky;
    dataspace_dimsFields_r[2] = array_global_size.nz+2;



    chunk_dimsFields[0] = array_local_size.nkx;
    chunk_dimsFields[1] = array_local_size.nky;
    chunk_dimsFields[2] = array_local_size.nkz;
    if(VERBOSE) printf("[MPI process %d] %llu %llu %llu\n",mpi_my_rank, chunk_dimsFields[0],chunk_dimsFields[1],chunk_dimsFields[2]);

    chunk_dimsFields_r[0] = array_local_size.nkx;
    chunk_dimsFields_r[1] = array_local_size.nky;
    chunk_dimsFields_r[2] = array_local_size.nz + 2;

    offsetFields[0] = array_local_size.nkx * mpi_my_coords[1];
    offsetFields[1] = 0;
    offsetFields[2] = 0;

    offsetFields_r[0] = array_local_size.nkx * mpi_my_coords[1];
    offsetFields_r[1] = 0;
    offsetFields_r[2] = 0;

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
    if(VERBOSE) printf("[process id %d] trying to create a file to write fields\n",mpi_my_rank);
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id); //creating new file collectively
    if(VERBOSE) printf("[process id %d] file created\n",mpi_my_rank);
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
 *  hdf_saveField_r
 * *************************/
void hdf_saveField_r(double *f, char *filename){
    //transpose data before saving
    fftw_copyFieldBuf_r(fftw_field,f);
    fftw_transposeToXY_field();
    // now we will create a file, and write a dataset into it.
    hid_t file_id, dset_id;
    hid_t file_space, memory_space;
    hid_t plist_id; //property list id
    plist_id = H5Pcreate(H5P_FILE_ACCESS); // access property list
    H5Pset_fapl_mpio(plist_id, mpi_row_comm, info);
    if(VERBOSE) printf("[process id %d] trying to create a file to write fields\n",mpi_my_rank);
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id); //creating new file collectively
    if(VERBOSE) printf("[process id %d] file created\n",mpi_my_rank);
    H5Pclose(plist_id);
    file_space = H5Screate_simple(hdf_rankFields, dataspace_dimsFields_r, NULL);
    memory_space = H5Screate_simple(hdf_rankFields, chunk_dimsFields_r, NULL);

    plist_id = H5Pcreate(H5P_DATASET_CREATE); //creating chunked dataset
    H5Pset_chunk(plist_id, hdf_rankFields, chunk_dimsFields_r);
    dset_id = H5Dcreate(file_id,"f", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    H5Sclose(file_space);

    file_space = H5Dget_space(dset_id);
    status = H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offsetFields_r, strideFields, countFields, chunk_dimsFields_r);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memory_space, file_space,plist_id, fftw_field);

    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Sclose(memory_space);
    H5Pclose(plist_id);
    H5Fclose(file_id);

    // transposing data to original shape
    fftw_transposeToYX_field();

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
    if(VERBOSE) printf("[process id %d] trying to create a file to write fields\n",mpi_my_rank);
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id); //creating new file collectively
    if(VERBOSE) printf("[process id %d] file created\n",mpi_my_rank);
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
    if(VERBOSE) printf("[process id %d] trying to create a file to write fields\n",mpi_my_rank);
    file_id = H5Fcreate(filename, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id); //creating new file collectively
    if(VERBOSE) printf("[process id %d] file created\n",mpi_my_rank);
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
    /*write injected energy to the file*/
    dset_id = H5Dopen2(file_id, "/freeEnergy/injected", H5P_DEFAULT);
    /*open a dataset*/
    dspace_id = H5Dget_space(dset_id);
    H5Dset_extent(dset_id, size);
    dspace_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, NULL, dims_ext, NULL);
    memspace = H5Screate_simple(1,dims_ext,NULL);
    if(mpi_my_rank == 0)
    {
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace, dspace_id,plist_id,&diag_injected);
    }
    H5Sclose(dspace_id);
    H5Dclose(dset_id);
    /*write dissipated energy to the file*/
    dset_id = H5Dopen2(file_id, "/freeEnergy/dissipated", H5P_DEFAULT);
    /*open a dataset*/
    dspace_id = H5Dget_space(dset_id);
    H5Dset_extent(dset_id, size);
    dspace_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, NULL, dims_ext, NULL);
    memspace = H5Screate_simple(1,dims_ext,NULL);
    if(mpi_my_rank == 0)
    {
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace, dspace_id,plist_id,&diag_dissipated);
    }
    H5Sclose(dspace_id);
    H5Dclose(dset_id);
    /*write dissipated kPerp energy to the file*/
    dset_id = H5Dopen2(file_id, "/freeEnergy/dissipated_kPerp", H5P_DEFAULT);
    /*open a dataset*/
    dspace_id = H5Dget_space(dset_id);
    H5Dset_extent(dset_id, size);
    dspace_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, NULL, dims_ext, NULL);
    memspace = H5Screate_simple(1,dims_ext,NULL);
    if(mpi_my_rank == 0)
    {
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace, dspace_id,plist_id,&diag_dissipated_kPerp);
    }
    H5Sclose(dspace_id);
    H5Dclose(dset_id);
    /*write dissipated kZ energy to the file*/
    dset_id = H5Dopen2(file_id, "/freeEnergy/dissipated_kZ", H5P_DEFAULT);
    /*open a dataset*/
    dspace_id = H5Dget_space(dset_id);
    H5Dset_extent(dset_id, size);
    dspace_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, NULL, dims_ext, NULL);
    memspace = H5Screate_simple(1,dims_ext,NULL);
    if(mpi_my_rank == 0)
    {
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace, dspace_id,plist_id,&diag_dissipated_kZ);
    }
    H5Sclose(dspace_id);
    H5Dclose(dset_id);
    /*write dissipated m energy to the file*/
    dset_id = H5Dopen2(file_id, "/freeEnergy/dissipated_m", H5P_DEFAULT);
    /*open a dataset*/
    dspace_id = H5Dget_space(dset_id);
    H5Dset_extent(dset_id, size);
    dspace_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, NULL, dims_ext, NULL);
    memspace = H5Screate_simple(1,dims_ext,NULL);
    if(mpi_my_rank == 0)
    {
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace, dspace_id,plist_id,&diag_dissipated_m);
    }
    H5Sclose(dspace_id);
    H5Dclose(dset_id);
    /*write h energy to the file*/
    dset_id = H5Dopen2(file_id, "/freeEnergy/hEnergy", H5P_DEFAULT);
    /*open a dataset*/
    dspace_id = H5Dget_space(dset_id);
    H5Dset_extent(dset_id, size);
    dspace_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, NULL, dims_ext, NULL);
    memspace = H5Screate_simple(1,dims_ext,NULL);
    if(mpi_my_rank == 0)
    {
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace, dspace_id,plist_id,&diag_energyH);
    }
    H5Sclose(dspace_id);
    H5Dclose(dset_id);
    /*write phi energy to the file*/
    dset_id = H5Dopen2(file_id, "/freeEnergy/phiEnergy", H5P_DEFAULT);
    /*open a dataset*/
    dspace_id = H5Dget_space(dset_id);
    H5Dset_extent(dset_id, size);
    dspace_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, NULL, dims_ext, NULL);
    memspace = H5Screate_simple(1,dims_ext,NULL);
    if(mpi_my_rank == 0)
    {
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace, dspace_id,plist_id,&diag_energyPhi);
    }
    H5Sclose(dspace_id);
    H5Dclose(dset_id);
    if (systemType == ELECTROMAGNETIC){
        /*write B_perp energy to the file*/
        dset_id = H5Dopen2(file_id, "/freeEnergy/BperpEnergy", H5P_DEFAULT);
        /*open a dataset*/
        dspace_id = H5Dget_space(dset_id);
        H5Dset_extent(dset_id, size);
        dspace_id = H5Dget_space(dset_id);
        H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, NULL, dims_ext, NULL);
        memspace = H5Screate_simple(1,dims_ext,NULL);
        if(mpi_my_rank == 0)
        {
            plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace, dspace_id,plist_id,&diag_energyBperp);
        }
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        /*write B_par energy to the file*/
        dset_id = H5Dopen2(file_id, "/freeEnergy/BparEnergy", H5P_DEFAULT);
        /*open a dataset*/
        dspace_id = H5Dget_space(dset_id);
        H5Dset_extent(dset_id, size);
        dspace_id = H5Dget_space(dset_id);
        H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, NULL, dims_ext, NULL);
        memspace = H5Screate_simple(1,dims_ext,NULL);
        if(mpi_my_rank == 0)
        {
            plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace, dspace_id,plist_id,&diag_energyBpar);
        }
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
    }

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
    /*write a timestep dataset*/
    /*opening a group*/
    dset_id = H5Dopen2(file_id, "/freeEnergy/time", H5P_DEFAULT);
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
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace, dspace_id,plist_id,&solver.curTime);
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
        if (mpi_my_rank == 0) printf("saving fields\n");
        hdf_saveFields(timestep);
        if (mpi_my_rank == 0) printf("fields are saved\n");
    }
    if (parameters.checkpoints && timestep % parameters.iter_checkpoint == 0)
    {
        if (mpi_my_rank == 0) printf("saving checkpoint\n");
        hdf_createCheckpoint(h, timestep);
        if (mpi_my_rank == 0) printf("checkpoint is saved\n");
    }
    if (parameters.save_distrib && timestep % parameters.iter_distribution == 0)
    {
        if (mpi_my_rank == 0) printf("saving distribution function\n");
        hdf_saveDistrib(h,timestep);
        if (mpi_my_rank == 0) printf("distribution function is saved\n");
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
                          "/geometry",
                          H5P_DEFAULT,
                          H5P_DEFAULT,
                          H5P_DEFAULT);
    /*create a datasets for storing kx,ky and kz*/
    /*create ky dataset*/
    dims[0] = array_local_size.nky;
    dspace_id = H5Screate_simple(1, dims, NULL);
    dset_id = H5Dcreate2(file_id, "/geometry/ky", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT,
                            H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, space_ky);
    H5Sclose(dspace_id);
    H5Dclose(dset_id);
    /*create kz dataset*/
    dims[0] = array_local_size.nkz;
    dspace_id = H5Screate_simple(1, dims, NULL);
    dset_id = H5Dcreate2(file_id, "/geometry/kz", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT,
                         H5P_DEFAULT);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, space_kz);
    H5Sclose(dspace_id);
    H5Dclose(dset_id);
    /*create kx dataset*/
    dims[0] = array_local_size.nkx;
    dimsf[0] = array_global_size.nkx;
    offsetKx[0] = offset[0];
    dspace_id = H5Screate_simple(1, dimsf, NULL);
    dset_id = H5Dcreate2(file_id, "/geometry/kx", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT,
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

    //writing Lx
    int rank_L = 1;
    hid_t dims_L[1] = {1};

    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_L, dims_L, NULL);
    dset_id = H5Dcreate(file_id,"geometry/Lx", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_L, dims_L, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dwrite(dset_id,H5T_NATIVE_DOUBLE, memory_space, file_space, plist_id, &parameters.Lx);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);
    //writing Ly
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_L, dims_L, NULL);
    dset_id = H5Dcreate(file_id,"geometry/Ly", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_L, dims_L, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dwrite(dset_id,H5T_NATIVE_DOUBLE, memory_space, file_space, plist_id, &parameters.Ly);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);
    //writing Lz
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_L, dims_L, NULL);
    dset_id = H5Dcreate(file_id,"geometry/Lz", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_L, dims_L, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dwrite(dset_id,H5T_NATIVE_DOUBLE, memory_space, file_space, plist_id, &parameters.Lz);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);
    //writing dx
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_L, dims_L, NULL);
    dset_id = H5Dcreate(file_id,"geometry/dx", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_L, dims_L, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dwrite(dset_id,H5T_NATIVE_DOUBLE, memory_space, file_space, plist_id, &space_dx);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);
    //writing dx
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_L, dims_L, NULL);
    dset_id = H5Dcreate(file_id,"geometry/dy", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_L, dims_L, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dwrite(dset_id,H5T_NATIVE_DOUBLE, memory_space, file_space, plist_id, &space_dy);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);
    //writing dz
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_L, dims_L, NULL);
    dset_id = H5Dcreate(file_id,"geometry/dz", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_L, dims_L, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dwrite(dset_id,H5T_NATIVE_DOUBLE, memory_space, file_space, plist_id, &space_dz);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);
    /* Close the group. */
    H5Gclose(group_id);

    /* Create a group to save physical parameters: charge, temperature, mass, density, beta, J00 and J10*/
    group_id = H5Gcreate2(file_id,
                          "/parameters",
                          H5P_DEFAULT,
                          H5P_DEFAULT,
                          H5P_DEFAULT);
    /*writing J00*/
    hid_t local_dims[3] = {array_local_size.nkx,array_local_size.nky,array_local_size.ns};
    hid_t glob_dims[3] = {array_global_size.nkx,array_global_size.nky,array_global_size.ns};
    hid_t offset_dims[3] = {offset[0],0,0};
    dspace_id = H5Screate_simple(3, glob_dims, NULL);
    dset_id = H5Dcreate2(file_id, "parameters/J00", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT,
                         H5P_DEFAULT);
    H5Sclose(dspace_id);
    dspace_id = H5Screate_simple(3, local_dims, NULL);
    filespace = H5Dget_space(dset_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset_dims, NULL, local_dims, NULL);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, dspace_id, filespace, plist_id, var_J0);


    H5Sclose(dspace_id);
    H5Sclose(filespace);
    H5Dclose(dset_id);

    /*writing J10*/
    dspace_id = H5Screate_simple(3, glob_dims, NULL);
    dset_id = H5Dcreate2(file_id, "parameters/J10", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, H5P_DEFAULT,
                         H5P_DEFAULT);
    H5Sclose(dspace_id);
    dspace_id = H5Screate_simple(3, local_dims, NULL);
    filespace = H5Dget_space(dset_id);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset_dims, NULL, local_dims, NULL);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, dspace_id, filespace, plist_id, var_J1);


    H5Sclose(dspace_id);
    H5Sclose(filespace);
    H5Dclose(dset_id);


    /* writing beta information */
    int rank_beta = 1;
    hid_t dims_beta[1] = {1};

    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_beta, dims_beta, NULL);
    dset_id = H5Dcreate(file_id,"parameters/beta", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_beta, dims_beta, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id,H5T_NATIVE_DOUBLE, memory_space, file_space, plist_id, &parameters.beta);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);

    int rank_particle = 1;
    hsize_t *dims_particle = malloc( sizeof(*dims_particle));
    dims_particle[0] = parameters.ns;
    /* writing charge */
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_particle, dims_particle, NULL);
    dset_id = H5Dcreate(file_id,"parameters/charge", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
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
    /* writing temperature */
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_particle, dims_particle, NULL);
    dset_id = H5Dcreate(file_id,"parameters/temperature", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_particle, dims_particle, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dwrite(dset_id,H5T_NATIVE_DOUBLE, memory_space, file_space, plist_id, parameters.temperature);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);
    /* writing density information */
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_particle, dims_particle, NULL);
    dset_id = H5Dcreate(file_id,"parameters/density", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
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
    dset_id = H5Dcreate(file_id,"parameters/mass", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
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

    /* Close the group. */
    H5Gclose(group_id);

    /* Create a group to save forcing */
    group_id = H5Gcreate2(file_id,
                          "/force",
                          H5P_DEFAULT,
                          H5P_DEFAULT,
                          H5P_DEFAULT);
    /* writing forcing power information */
    int rank_power = 1;
    hid_t dims_power[1] = {1};

    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_power, dims_power, NULL);
    dset_id = H5Dcreate(file_id,"force/power", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_power, dims_power, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    //H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id,H5T_NATIVE_DOUBLE, memory_space, file_space, plist_id, &parameters.forcePower);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);

    /* writing forcing shell boundaries */
    int rank_forceShell = 1;
    hid_t dims_forceShell[1] = {2};
    double forceShell[2] = {parameters.forceKmin,parameters.forceKmax};

    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_forceShell, dims_forceShell, NULL);
    dset_id = H5Dcreate(file_id,"force/forceShell", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_forceShell, dims_forceShell, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dwrite(dset_id,H5T_NATIVE_DOUBLE, memory_space, file_space, plist_id, &forceShell);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);

    /* writing kx indices which are forced */
    int rank_forcedIndices = 1;
    hid_t dims_forcedIndices[1] = {equation_forceNorm};

    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_forcedIndices, dims_forcedIndices, NULL);
    dset_id = H5Dcreate(file_id,"force/kx", H5T_NATIVE_INT, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_forcedIndices, dims_forcedIndices, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dwrite(dset_id,H5T_NATIVE_INT, memory_space, file_space, plist_id, equation_forceKxIndGathered);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);

    /* writing ky indices which are forced */
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_forcedIndices, dims_forcedIndices, NULL);
    dset_id = H5Dcreate(file_id,"force/ky", H5T_NATIVE_INT, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_forcedIndices, dims_forcedIndices, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dwrite(dset_id,H5T_NATIVE_INT, memory_space, file_space, plist_id, equation_forceKyIndGathered);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);

    /* writing kz indices which are forced */
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_forcedIndices, dims_forcedIndices, NULL);
    dset_id = H5Dcreate(file_id,"force/kz", H5T_NATIVE_INT, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_forcedIndices, dims_forcedIndices, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Dwrite(dset_id,H5T_NATIVE_INT, memory_space, file_space, plist_id, equation_forceKzIndGathered);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);

    /*creating dataset in which forcing will be stored*/
    hid_t dims_forceArray[3] = {0, equation_forceNorm, array_local_size.ns};
    hid_t chunk_forceArray[3] = {1, equation_forceNorm, array_local_size.ns};
    hid_t max_forceArray[3] = {H5S_UNLIMITED, equation_forceNorm, array_local_size.ns};

    plist_id   = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id, 3, chunk_forceArray);
    dspace_id = H5Screate_simple(3, dims_forceArray, max_forceArray);
    dset_id = H5Dcreate2(file_id,
                         "/force/forcing",
                         complex_id,
                         dspace_id,
                         H5P_DEFAULT,
                         plist_id,
                         H5P_DEFAULT);
    H5Sclose(dspace_id);
    H5Dclose(dset_id);
    H5Pclose(plist_id);

    /* Close the group. */
    H5Gclose(group_id);


    /*creating a group to save spectra*/
    if(parameters.save_diagnostics || parameters.compute_k)
    {
        group_id = H5Gcreate2(file_id, "/spectra", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Gclose(group_id);

        hid_t dims_spec_k[2] = {0,diag_numOfShells};
        hid_t dims_spec_k2D[3] = {0,diag_numOfShells,parameters.nm};
        hid_t dims_spec_m[2] = {0,parameters.nm};
        hid_t chunk_spec_k[2] = {1, diag_numOfShells};
        hid_t chunk_spec_k2D[3] = {1, diag_numOfShells, parameters.nm};
        hid_t chunk_spec_k_borders[2] = {1,diag_numOfShellBounds};
        hid_t chunk_spec_m[2] = {1,parameters.nm};
        hid_t max_dims_k[2] = {H5S_UNLIMITED,diag_numOfShells};
        hid_t max_dims_k2D[3] = {H5S_UNLIMITED,diag_numOfShells,parameters.nm};
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
        dset_id = H5Dcreate2(file_id, "/spectra/time", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, plist_id,
                             H5P_DEFAULT);
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);

        /* creating L_dissipation dataset*/
        plist_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 1, &chunk_spec_k[0]);
        dspace_id = H5Screate_simple(1,&dims_spec_k[0],&max_dims_k[0]);
        dset_id = H5Dcreate2(file_id, "/spectra/L_dissipation", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, plist_id,
                             H5P_DEFAULT);
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);

        /* creating L_integral dataset*/
        plist_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 1, &chunk_spec_k[0]);
        dspace_id = H5Screate_simple(1,&dims_spec_k[0],&max_dims_k[0]);
        dset_id = H5Dcreate2(file_id, "/spectra/L_integral", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, plist_id,
                             H5P_DEFAULT);
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);

        if(parameters.compute_k)
        {
            /* creating shells borders dataset */
            plist_id = H5Pcreate(H5P_DATASET_CREATE);
            H5Pset_chunk(plist_id, 1, &chunk_spec_k_borders[1]);
            dspace_id = H5Screate_simple(1,&chunk_spec_k_borders[1],&chunk_spec_k_borders[1]);
            dset_id = H5Dcreate2(file_id, "/spectra/shellBorders", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, plist_id,
                                 H5P_DEFAULT);
            //printf("%f",diag_shells[2]);
            H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, diag_shells);
            H5Sclose(dspace_id);
            H5Dclose(dset_id);
            H5Pclose(plist_id);

            /* creating shells centres dataset */
            plist_id = H5Pcreate(H5P_DATASET_CREATE);
            H5Pset_chunk(plist_id, 1, &chunk_spec_k[1]);
            dspace_id = H5Screate_simple(1,&chunk_spec_k[1],&chunk_spec_k[1]);
            dset_id = H5Dcreate2(file_id, "/spectra/shellCentres", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, plist_id,
                                 H5P_DEFAULT);
            H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, diag_shellCentres);
            H5Sclose(dspace_id);
            H5Dclose(dset_id);
            H5Pclose(plist_id);

            /* creating normalizations dataset dataset */
            plist_id = H5Pcreate(H5P_DATASET_CREATE);
            H5Pset_chunk(plist_id, 1, &chunk_spec_k[1]);
            dspace_id = H5Screate_simple(1,&chunk_spec_k[1],&chunk_spec_k[1]);
            dset_id = H5Dcreate2(file_id, "/spectra/norm", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, plist_id,
                                 H5P_DEFAULT);
            H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, diag_shellNorm);
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

            /*creating a kSpecH dataset*/
            plist_id   = H5Pcreate(H5P_DATASET_CREATE);
            H5Pset_chunk(plist_id, 3, chunk_spec_k2D);
            dspace_id = H5Screate_simple(3,dims_spec_k2D,max_dims_k2D);
            dset_id = H5Dcreate2(file_id,
                                 "/spectra/kSpecH",
                                 H5T_NATIVE_DOUBLE,
                                 dspace_id,
                                 H5P_DEFAULT,
                                 plist_id,
                                 H5P_DEFAULT);
            H5Sclose(dspace_id);
            H5Dclose(dset_id);
            H5Pclose(plist_id);

            /*creating a kSpecPhi dataset*/
            plist_id   = H5Pcreate(H5P_DATASET_CREATE);
            H5Pset_chunk(plist_id, 2, chunk_spec_k);
            dspace_id = H5Screate_simple(2,dims_spec_k,max_dims_k);
            dset_id = H5Dcreate2(file_id,
                                 "/spectra/kSpecPhi",
                                 H5T_NATIVE_DOUBLE,
                                 dspace_id,
                                 H5P_DEFAULT,
                                 plist_id,
                                 H5P_DEFAULT);
            H5Sclose(dspace_id);
            H5Dclose(dset_id);
            H5Pclose(plist_id);
            /*creating a kSpecBperp kSpecBperp dataset*/
            if (systemType == ELECTROMAGNETIC){
                /*creating a kSpecBperp dataset*/
                plist_id   = H5Pcreate(H5P_DATASET_CREATE);
                H5Pset_chunk(plist_id, 2, chunk_spec_k);
                dspace_id = H5Screate_simple(2,dims_spec_k,max_dims_k);
                dset_id = H5Dcreate2(file_id,
                                     "/spectra/kSpecBperp",
                                     H5T_NATIVE_DOUBLE,
                                     dspace_id,
                                     H5P_DEFAULT,
                                     plist_id,
                                     H5P_DEFAULT);
                H5Sclose(dspace_id);
                H5Dclose(dset_id);
                H5Pclose(plist_id);
                /*creating a kSpecBpar dataset*/
                plist_id   = H5Pcreate(H5P_DATASET_CREATE);
                H5Pset_chunk(plist_id, 2, chunk_spec_k);
                dspace_id = H5Screate_simple(2,dims_spec_k,max_dims_k);
                dset_id = H5Dcreate2(file_id,
                                     "/spectra/kSpecBpar",
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
    /*creating a group to save nonlinear flux*/
    if(parameters.compute_nonlinear){
        hid_t chunk_spec_nonl_flux_shells[2] = {1,diag_numOfShells};
        hid_t max_dims_flux[2] = {H5S_UNLIMITED,diag_numOfShells};
        hid_t dims_spec_flux[2] = {0, diag_numOfShells};
        hid_t chunk_flux_timestep[1] = {1};
        hid_t maxdims_flux_timestep[1] = {H5S_UNLIMITED};
        hid_t dims_flux_timestep[1] = {0};
        group_id = H5Gcreate2(file_id, "/nonlinearFlux", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        H5Gclose(group_id);

        /*create nonlinear flux dataset, shells dataset, timestep dataset and dt dataset*/
        /*creating nonlinear flux dataset*/
        plist_id   = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 2, chunk_spec_nonl_flux_shells);
        dspace_id = H5Screate_simple(2,dims_spec_flux,max_dims_flux);
        dset_id = H5Dcreate2(file_id,
                             "/nonlinearFlux/flux",
                             H5T_NATIVE_DOUBLE,
                             dspace_id,
                             H5P_DEFAULT,
                             plist_id,
                             H5P_DEFAULT);
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);
        /*creating nonlinear inverse flux dataset*/
        plist_id   = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 2, chunk_spec_nonl_flux_shells);
        dspace_id = H5Screate_simple(2,dims_spec_flux,max_dims_flux);
        dset_id = H5Dcreate2(file_id,
                             "/nonlinearFlux/fluxInverse",
                             H5T_NATIVE_DOUBLE,
                             dspace_id,
                             H5P_DEFAULT,
                             plist_id,
                             H5P_DEFAULT);
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);

        /*creating nonlinear forward flux dataset*/
        plist_id   = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 2, chunk_spec_nonl_flux_shells);
        dspace_id = H5Screate_simple(2,dims_spec_flux,max_dims_flux);
        dset_id = H5Dcreate2(file_id,
                             "/nonlinearFlux/fluxForward",
                             H5T_NATIVE_DOUBLE,
                             dspace_id,
                             H5P_DEFAULT,
                             plist_id,
                             H5P_DEFAULT);
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);

        /* creating shells borders dataset */
        plist_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 1, &chunk_spec_nonl_flux_shells[1]);
        dspace_id = H5Screate_simple(1,&chunk_spec_nonl_flux_shells[1],&chunk_spec_nonl_flux_shells[1]);
        dset_id = H5Dcreate2(file_id, "/nonlinearFlux/shells", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, plist_id,
                             H5P_DEFAULT);
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &diag_shells[1]);
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);

        /* creating shells norms dataset */
        plist_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 1, &chunk_spec_nonl_flux_shells[1]);
        dspace_id = H5Screate_simple(1,&chunk_spec_nonl_flux_shells[1],&chunk_spec_nonl_flux_shells[1]);
        dset_id = H5Dcreate2(file_id, "/nonlinearFlux/norm", H5T_NATIVE_DOUBLE, dspace_id, H5P_DEFAULT, plist_id,
                             H5P_DEFAULT);
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, diag_nonlinearNorm);
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);

        /* creating timestep dataset */
        plist_id   = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 1, chunk_flux_timestep);
        dspace_id = H5Screate_simple(1,dims_flux_timestep,maxdims_flux_timestep);
        dset_id = H5Dcreate2(file_id, "/nonlinearFlux/timestep",
                             H5T_NATIVE_INT,
                             dspace_id,
                             H5P_DEFAULT,
                             plist_id,
                             H5P_DEFAULT);
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);
        /* creating dt dataset */
        plist_id   = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 1, chunk_flux_timestep);
        dspace_id = H5Screate_simple(1,dims_flux_timestep,maxdims_flux_timestep);
        dset_id = H5Dcreate2(file_id, "/nonlinearFlux/time",
                             H5T_NATIVE_DOUBLE,
                             dspace_id,
                             H5P_DEFAULT,
                             plist_id,
                             H5P_DEFAULT);
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);
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
        dset_id = H5Dcreate2(file_id, "/freeEnergy/freeEnergy",
                             H5T_NATIVE_DOUBLE,
                             dspace_id,
                             H5P_DEFAULT,
                             plist_id,
                             H5P_DEFAULT);
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);

        /*creating timestep dataset*/
        plist_id   = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 1, chunk_energy);
        dspace_id = H5Screate_simple(1,dims_energy,maxdims);
        dset_id = H5Dcreate2(file_id, "/freeEnergy/timestep",
                             H5T_NATIVE_INT,
                             dspace_id,
                             H5P_DEFAULT,
                             plist_id,
                             H5P_DEFAULT);
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);

        /*creating dt dataset*/
        plist_id  = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 1, chunk_energy);
        dspace_id = H5Screate_simple(1,dims_energy,maxdims);
        dset_id = H5Dcreate2(file_id, "/freeEnergy/time",
                             H5T_NATIVE_DOUBLE,
                             dspace_id,
                             H5P_DEFAULT,
                             plist_id,
                             H5P_DEFAULT);
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);

        /*creating injected dataset*/
        plist_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 1, chunk_energy);
        dspace_id = H5Screate_simple(1,dims_energy,maxdims);
        dset_id = H5Dcreate2(file_id, "/freeEnergy/injected",
                             H5T_NATIVE_DOUBLE,
                             dspace_id,
                             H5P_DEFAULT,
                             plist_id,
                             H5P_DEFAULT);
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);

        /*creating dissipated dataset*/
        plist_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 1, chunk_energy);
        dspace_id = H5Screate_simple(1,dims_energy,maxdims);
        dset_id = H5Dcreate2(file_id, "/freeEnergy/dissipated",
                             H5T_NATIVE_DOUBLE,
                             dspace_id,
                             H5P_DEFAULT,
                             plist_id,
                             H5P_DEFAULT);
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);

        /*creating dissipated by perp dissipation dataset*/
        plist_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 1, chunk_energy);
        dspace_id = H5Screate_simple(1,dims_energy,maxdims);
        dset_id = H5Dcreate2(file_id, "/freeEnergy/dissipated_kPerp",
                             H5T_NATIVE_DOUBLE,
                             dspace_id,
                             H5P_DEFAULT,
                             plist_id,
                             H5P_DEFAULT);
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);

        /*creating dissipated by kz dissipation dataset*/
        plist_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 1, chunk_energy);
        dspace_id = H5Screate_simple(1,dims_energy,maxdims);
        dset_id = H5Dcreate2(file_id, "/freeEnergy/dissipated_kZ",
                             H5T_NATIVE_DOUBLE,
                             dspace_id,
                             H5P_DEFAULT,
                             plist_id,
                             H5P_DEFAULT);
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);

        /*creating dissipated by m dissipation dataset*/
        plist_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 1, chunk_energy);
        dspace_id = H5Screate_simple(1,dims_energy,maxdims);
        dset_id = H5Dcreate2(file_id, "/freeEnergy/dissipated_m",
                             H5T_NATIVE_DOUBLE,
                             dspace_id,
                             H5P_DEFAULT,
                             plist_id,
                             H5P_DEFAULT);
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);

        /*creating H energy dataset*/
        plist_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 1, chunk_energy);
        dspace_id = H5Screate_simple(1,dims_energy,maxdims);
        dset_id = H5Dcreate2(file_id, "/freeEnergy/hEnergy",
                             H5T_NATIVE_DOUBLE,
                             dspace_id,
                             H5P_DEFAULT,
                             plist_id,
                             H5P_DEFAULT);
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);

        /*creating phi energy dataset*/
        plist_id = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 1, chunk_energy);
        dspace_id = H5Screate_simple(1,dims_energy,maxdims);
        dset_id = H5Dcreate2(file_id, "/freeEnergy/phiEnergy",
                             H5T_NATIVE_DOUBLE,
                             dspace_id,
                             H5P_DEFAULT,
                             plist_id,
                             H5P_DEFAULT);
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);
        if (systemType = ELECTROMAGNETIC){
            /*creating B_perp energy dataset*/
            plist_id = H5Pcreate(H5P_DATASET_CREATE);
            H5Pset_chunk(plist_id, 1, chunk_energy);
            dspace_id = H5Screate_simple(1,dims_energy,maxdims);
            dset_id = H5Dcreate2(file_id, "/freeEnergy/BperpEnergy",
                                 H5T_NATIVE_DOUBLE,
                                 dspace_id,
                                 H5P_DEFAULT,
                                 plist_id,
                                 H5P_DEFAULT);
            H5Sclose(dspace_id);
            H5Dclose(dset_id);
            H5Pclose(plist_id);

            /*creating B_par energy dataset*/
            plist_id = H5Pcreate(H5P_DATASET_CREATE);
            H5Pset_chunk(plist_id, 1, chunk_energy);
            dspace_id = H5Screate_simple(1,dims_energy,maxdims);
            dset_id = H5Dcreate2(file_id, "/freeEnergy/BparEnergy",
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
    hid_t dims_ext[2] = {1,diag_numOfShells};
    hid_t dims_ext2D[3] = {1,diag_numOfShells,array_local_size.nm};
    hid_t size[2];
    hid_t size2D[3];
    hid_t offset[2];
    hid_t offset2D[3];
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

    /* writing h spec to the file*/
    /*opening a group and a dataset*/
    /*opening a group*/
    dset_id = H5Dopen2(file_id, "/spectra/kSpecH", H5P_DEFAULT);
    dspace_id = H5Dget_space(dset_id);
    /*get dataset's dimensions */
    int ndims2D = H5Sget_simple_extent_ndims(dspace_id);
    hsize_t *dims2D = malloc(ndims * sizeof(*dims2D));
    H5Sget_simple_extent_dims(dspace_id, dims2D, NULL);
    H5Sclose(dspace_id);
    /*extend dataset size*/
    size2D[0] = dims2D[0] + dims_ext2D[0];
    size2D[1] = dims2D[1];
    size2D[2] = dims2D[2];
    offset2D[0] = dims2D[0];
    offset2D[1] = 0;
    offset2D[2] = mpi_my_col_rank * parameters.nm / mpi_dims[0];
    /*extent dataset dims*/
    H5Dset_extent(dset_id, size2D);
    /*write free energy to the file*/
    dspace_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset2D, NULL, dims_ext2D, NULL);
    memspace = H5Screate_simple(3,dims_ext2D,NULL);

    if(mpi_my_row_rank == 0)
    {
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace, dspace_id,plist_id,diag_kSpecH);
    }
    H5Sclose(dspace_id);
    H5Dclose(dset_id);

    /* writing phi spec to the file*/
    /*opening a group and a dataset*/
    /*opening a group*/
    dset_id = H5Dopen2(file_id, "/spectra/kSpecPhi", H5P_DEFAULT);
    /*extent dataset dims*/
    H5Dset_extent(dset_id, size);
    /*write free energy to the file*/
    dspace_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, NULL, dims_ext, NULL);
    memspace = H5Screate_simple(2,dims_ext,NULL);

    if(mpi_my_rank == 0)
    {
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace, dspace_id,plist_id,diag_kSpecPhi);
    }
    H5Sclose(dspace_id);
    H5Dclose(dset_id);
    if (systemType == ELECTROMAGNETIC){
        /* writing Bpar spec to the file*/
        /*opening a group and a dataset*/
        /*opening a group*/
        dset_id = H5Dopen2(file_id, "/spectra/kSpecBpar", H5P_DEFAULT);
        /*extent dataset dims*/
        H5Dset_extent(dset_id, size);
        /*write free energy to the file*/
        dspace_id = H5Dget_space(dset_id);
        H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, NULL, dims_ext, NULL);
        memspace = H5Screate_simple(2,dims_ext,NULL);

        if(mpi_my_rank == 0)
        {
            plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace, dspace_id,plist_id,diag_kSpecBpar);
        }
        H5Sclose(dspace_id);
        H5Dclose(dset_id);

        /* writing Bperp spec to the file*/
        /*opening a group and a dataset*/
        /*opening a group*/
        dset_id = H5Dopen2(file_id, "/spectra/kSpecBperp", H5P_DEFAULT);
        /*extent dataset dims*/
        H5Dset_extent(dset_id, size);
        /*write free energy to the file*/
        dspace_id = H5Dget_space(dset_id);
        H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, NULL, dims_ext, NULL);
        memspace = H5Screate_simple(2,dims_ext,NULL);

        if(mpi_my_rank == 0)
        {
            plist_id = H5Pcreate(H5P_DATASET_XFER);
            H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace, dspace_id,plist_id,diag_kSpecBperp);
        }
        H5Sclose(dspace_id);
        H5Dclose(dset_id);
    }
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
    /*write a time dataset*/
    /*opening a group*/
    dset_id = H5Dopen2(file_id, "/spectra/time", H5P_DEFAULT);
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
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace, dspace_id,plist_id,&solver.curTime);
    }
    H5Sclose(dspace_id);
    H5Dclose(dset_id);

    /*write a L_dissipation dataset*/
    /*opening a group*/
    dset_id = H5Dopen2(file_id, "/spectra/L_dissipation", H5P_DEFAULT);
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
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace, dspace_id,plist_id,&diag_LDis);
    }
    H5Sclose(dspace_id);
    H5Dclose(dset_id);

    /*write a L_dissipation dataset*/
    /*opening a group*/
    dset_id = H5Dopen2(file_id, "/spectra/L_integral", H5P_DEFAULT);
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
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace, dspace_id,plist_id,&diag_LInt);
    }
    H5Sclose(dspace_id);
    H5Dclose(dset_id);
    H5Fclose(file_id);
};

/***************************
 *  hdf_saveNonlinearFlux
 * *************************/
void hdf_saveNonlinearFlux(int timestep) {
    hid_t file_id, dset_id,dspace_id,group_id,filespace,memspace;
    hid_t plist_id; //property list id
    hid_t dims_ext[2] = {1,diag_numOfShells};
    hid_t dims_ext2D[3] = {1,diag_numOfShells,array_local_size.nm};
    hid_t size[2];
    hid_t size2D[3];
    hid_t offset[2];
    hid_t offset2D[3];
    plist_id = H5Pcreate(H5P_FILE_ACCESS); // access property list
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
    /* open file to read and write */
    file_id = H5Fopen(PARAMETER_FILENAME,H5F_ACC_RDWR,plist_id);
    H5Pclose(plist_id);
    /*opening a group and a dataset*/
    /*opening a group*/
    dset_id = H5Dopen2(file_id, "/nonlinearFlux/flux", H5P_DEFAULT);
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
    /*write nonlinear flux to the file*/
    dspace_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, NULL, dims_ext, NULL);
    memspace = H5Screate_simple(2,dims_ext,NULL);
    if(mpi_my_rank == 0)
    {
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace, dspace_id,plist_id,diag_nonlinearFlux);
    }
    H5Sclose(dspace_id);
    H5Dclose(dset_id);

    //saving inverse flux now
    /*opening a group*/
    dset_id = H5Dopen2(file_id, "/nonlinearFlux/fluxInverse", H5P_DEFAULT);
    /*open a dataset*/
    dspace_id = H5Dget_space(dset_id);
    /*get dataset's dimensions */
    H5Sget_simple_extent_dims(dspace_id, dims, NULL);
    H5Sclose(dspace_id);
    /*extend dataset size*/
    size[0] = dims[0] + dims_ext[0];
    size[1] = dims[1];
    offset[0] = dims[0];
    offset[1] = 0;
    H5Dset_extent(dset_id, size);
    /*write nonlinear flux to the file*/
    dspace_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, NULL, dims_ext, NULL);
    memspace = H5Screate_simple(2,dims_ext,NULL);
    if(mpi_my_rank == 0)
    {
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace, dspace_id,plist_id,diag_nonlinearFluxInverse);
    }
    H5Sclose(dspace_id);
    H5Dclose(dset_id);

    //saving forward flux now
    /*opening a group*/
    dset_id = H5Dopen2(file_id, "/nonlinearFlux/fluxForward", H5P_DEFAULT);
    /*open a dataset*/
    dspace_id = H5Dget_space(dset_id);
    /*get dataset's dimensions */
    H5Sget_simple_extent_dims(dspace_id, dims, NULL);
    H5Sclose(dspace_id);
    /*extend dataset size*/
    size[0] = dims[0] + dims_ext[0];
    size[1] = dims[1];
    offset[0] = dims[0];
    offset[1] = 0;
    H5Dset_extent(dset_id, size);
    /*write nonlinear flux to the file*/
    dspace_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, NULL, dims_ext, NULL);
    memspace = H5Screate_simple(2,dims_ext,NULL);
    if(mpi_my_rank == 0)
    {
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace, dspace_id,plist_id,diag_nonlinearFluxForward);
    }
    H5Sclose(dspace_id);
    H5Dclose(dset_id);
    H5Fclose(file_id);


}

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
 *  hdf_saveForcing
 * *************************/
void hdf_saveForcing(){
    hid_t file_id, dset_id,dspace_id,group_id,filespace,memspace;
    hid_t plist_id; //property list id
    hid_t dims_ext_full[3] = {1, equation_forceNorm, array_local_size.ns};
    hid_t dims_ext_local[3] = {1, equation_forceKn, array_local_size.ns};
    hid_t max_dims_local[3] = {H5S_UNLIMITED, equation_forceKn, array_local_size.ns};
    hid_t max_dims_full[3] = {H5S_UNLIMITED, equation_forceNorm, array_local_size.ns};
    hid_t size[3];
    hid_t offset[3];
    hid_t stride[3] = {1,1,1};
    hid_t count[3] = {1,1,1};
    plist_id = H5Pcreate(H5P_FILE_ACCESS); // access property list
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
    /* open file to read and write */
    file_id = H5Fopen(PARAMETER_FILENAME,H5F_ACC_RDWR,plist_id);
    H5Pclose(plist_id);
    /*opening a group and a dataset*/
    /*opening a group*/
    dset_id = H5Dopen2(file_id, "/force/forcing", H5P_DEFAULT);
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
    offset[0] = dims[0];
    offset[1] = equation_displacements[mpi_my_row_rank];
    offset[2] = 0;
    H5Dset_extent(dset_id, size);
    /*write forcing to the file*/
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    dspace_id = H5Dget_space(dset_id);
    H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, stride,count, dims_ext_local);
    memspace = H5Screate_simple(2,&dims_ext_local[1],&max_dims_local[1]);
    int forced_proc = mpi_whereIsM[equation_forcedM];
    if (mpi_my_col_rank != forced_proc)
    {
        H5Sselect_none(dspace_id);
        H5Sselect_none(memspace);
    }
    H5Dwrite(dset_id, complex_id, memspace, dspace_id, plist_id, equation_forcingAr);
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
    hdf_dumpCheckpointReal(h, timestep, hdf_newCheckpointName);
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

    /* writing current time information */
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_timestep, dims_timestep, NULL);
    dset_id = H5Dcreate(file_id,"time", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_timestep, dims_timestep, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id,H5T_NATIVE_DOUBLE, memory_space, file_space, plist_id, &solver.curTime);
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

    /* writing Lx information */
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_timestep, dims_timestep, NULL);
    dset_id = H5Dcreate(file_id,"Lx", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_timestep, dims_timestep, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id,H5T_NATIVE_DOUBLE, memory_space, file_space, plist_id, &parameters.Lx);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);

    /* writing Ly information */
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_timestep, dims_timestep, NULL);
    dset_id = H5Dcreate(file_id,"Ly", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_timestep, dims_timestep, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id,H5T_NATIVE_DOUBLE, memory_space, file_space, plist_id, &parameters.Ly);
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Pclose(plist_id);

    /* writing Lz information */
    plist_id = H5Pcreate(H5P_DATASET_CREATE);
    file_space = H5Screate_simple(rank_timestep, dims_timestep, NULL);
    dset_id = H5Dcreate(file_id,"Lz", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    memory_space = H5Screate_simple(rank_timestep, dims_timestep, NULL);
    if(mpi_my_rank != CHECKPOINT_ROOT)
    {
        H5Sselect_none(file_space);
        H5Sselect_none(memory_space);
    }
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset_id,H5T_NATIVE_DOUBLE, memory_space, file_space, plist_id, &parameters.Lz);
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

/***************************
 *  hdf_dumpCheckpointReal
 * *************************/
void hdf_dumpCheckpointReal(COMPLEX *h, int timestep, char *filename){
    //fftw of h distribution function
    fftw_copy_buffer_c(fftw_hBuf, h);
    fftw_c2r();
    fftw_transposeToXY();
    double *hr = fftw_hBuf;
    // write a dataset into it.
    hid_t file_id, dset_id;
    hid_t file_space, memory_space;
    hid_t plist_id; //property list id
    hid_t dims_timestep[1] = {1};
    hid_t rank_timestep = 1;
    plist_id = H5Pcreate(H5P_FILE_ACCESS); // access property list
    H5Pset_fapl_mpio(plist_id, mpi_cube_comm, info);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, info);
    /* open file to read and write */
    file_id = H5Fopen(filename,H5F_ACC_RDWR,plist_id);
    H5Pclose(plist_id);
    /* create file space and memory space. Memory space holds each process' data, file space going to hold whole data */
    file_space = H5Screate_simple(hdf_rank, dataspace_dims_r, NULL);
    memory_space = H5Screate_simple(hdf_rank, chunk_dims_r, NULL);

    /*
     * creating h dataset
     */
    plist_id = H5Pcreate(H5P_DATASET_CREATE); //creating chunked dataset
    H5Pset_chunk(plist_id, hdf_rank, chunk_dims_c);
    dset_id = H5Dcreate(file_id,"hr", H5T_NATIVE_DOUBLE, file_space, H5P_DEFAULT, plist_id, H5P_DEFAULT);
    H5Pclose(plist_id);
    H5Sclose(file_space);
    /* writing dataset collectively */
    file_space = H5Dget_space(dset_id);
    status = H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset, stride, count, chunk_dims_r);
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memory_space, file_space,plist_id, hr);
    /*closing dataset identifier, memory space and file space*/
    H5Dclose(dset_id);
    H5Sclose(file_space);
    H5Sclose(memory_space);
    H5Pclose(plist_id);
    H5Fclose(file_id);

    // transposing data to original shape
    fftw_transposeToYX();
}

void hdf_saveDistrib(COMPLEX* h, int timestep){
    char distrib_name[FILENAME_ID_LEN];
    sprintf(distrib_name, "%s%d.h5", DISTRIBUTION_FILENAME, timestep);
    hdf_dumpCheckpoint(h,  timestep, distrib_name);
    hdf_dumpCheckpointReal(h,  timestep, distrib_name);
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

    /* variables needed to write real fields to the file  */
    hid_t dims_field_r[4] = {0,array_global_size.nx, array_global_size.ny, array_global_size.nz + 2};
    hid_t chunk_field_r[4] = {1,array_global_size.nx, array_global_size.ny, array_global_size.nz + 2};
    hid_t max_dims_field_r[4] = {H5S_UNLIMITED,array_global_size.nx, array_global_size.ny, array_global_size.nz + 2};

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
    /*creating time dataset*/
    plist_id   = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist_id, 1, &chunk_field[0]);
    dspace_id = H5Screate_simple(1,&dims_field[0],&max_dims_field[0]);
    dset_id = H5Dcreate2(file_id, "time", H5T_NATIVE_INT, dspace_id, H5P_DEFAULT, plist_id,
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

        /* create phi real dataset */
        plist_id   = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 4, chunk_field_r);
        dspace_id = H5Screate_simple(4,dims_field_r,max_dims_field_r);
        dset_id = H5Dcreate2(file_id,
                             "phi_r",
                             H5T_NATIVE_DOUBLE,
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

        /* create phi real dataset */
        plist_id   = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 4, chunk_field_r);
        dspace_id = H5Screate_simple(4,dims_field_r,max_dims_field_r);
        dset_id = H5Dcreate2(file_id,
                             "phi_r",
                             H5T_NATIVE_DOUBLE,
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

        /* create A real dataset */
        plist_id   = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 4, chunk_field_r);
        dspace_id = H5Screate_simple(4,dims_field_r,max_dims_field_r);
        dset_id = H5Dcreate2(file_id,
                             "A_r",
                             H5T_NATIVE_DOUBLE,
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

        /* create B real dataset */
        plist_id   = H5Pcreate(H5P_DATASET_CREATE);
        H5Pset_chunk(plist_id, 4, chunk_field_r);
        dspace_id = H5Screate_simple(4,dims_field_r,max_dims_field_r);
        dset_id = H5Dcreate2(file_id,
                             "B_r",
                             H5T_NATIVE_DOUBLE,
                             dspace_id,
                             H5P_DEFAULT,
                             plist_id,
                             H5P_DEFAULT);

        H5Sclose(dspace_id);
        H5Dclose(dset_id);
        H5Pclose(plist_id);
    }
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

    hid_t dims_ext_full_r[4] = {1,array_global_size.nkx, array_global_size.nky, array_global_size.nz + 2};
    hid_t dims_ext_local_r[4] = {1,array_local_size.nkx, array_local_size.nky, array_local_size.nz + 2};
    hid_t max_dims_local_r[4] = {H5S_UNLIMITED,array_local_size.nkx, array_local_size.nky, array_local_size.nz + 2};
    hid_t max_dims_full_r[4] = {H5S_UNLIMITED,array_local_size.nkx, array_local_size.nky, array_local_size.nz + 2};

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
        /*write phi field to the file*/
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

        /*writing phi dataset*/
        /*opening a group and a dataset*/
        /*opening a group*/
        dset_id = H5Dopen2(file_id, "phi_r", H5P_DEFAULT);
        /*open a dataset*/
        dspace_id = H5Dget_space(dset_id);
        /*get dataset's dimensions */
        ndims = H5Sget_simple_extent_ndims(dspace_id);
        H5Sget_simple_extent_dims(dspace_id, dims, NULL);
        H5Sclose(dspace_id);
        /*extend dataset size*/
        size[0] = dims[0] + dims_ext_full_r[0];
        size[1] = dims[1];
        size[2] = dims[2];
        size[3] = dims[3];
        offset[0] = dims[0];
        offset[1] =  mpi_my_row_rank * parameters.nkx / mpi_dims[1];
        offset[2] =  0;
        offset[3] =  0;
        H5Dset_extent(dset_id, size);

        /*transform phi to real space*/
        fftw_copyFieldBuf_c(fftw_field, fields_fields.phi);
        fftw_c2r_field();
        fftw_transposeToXY_field();

        /*write phi real field to the file*/
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        dspace_id = H5Dget_space(dset_id);
        H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, stride,count, dims_ext_local_r);
        memspace = H5Screate_simple(3,&dims_ext_local_r[1],&max_dims_local_r[1]);
        if (mpi_my_col_rank != 0)
        {
            H5Sselect_none(dspace_id);
            H5Sselect_none(memspace);
        }
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, dspace_id, plist_id, fftw_field);
        H5Sclose(memspace);
        H5Sclose(dspace_id);
        H5Pclose(plist_id);
        H5Dclose(dset_id);

        // transposing data to original shape
        fftw_transposeToYX_field();
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

        /*writing phi real dataset*/
        /*opening a group and a dataset*/
        /*opening a group*/
        dset_id = H5Dopen2(file_id, "phi_r", H5P_DEFAULT);
        /*open a dataset*/
        dspace_id = H5Dget_space(dset_id);
        /*get dataset's dimensions */
        ndims = H5Sget_simple_extent_ndims(dspace_id);
        H5Sget_simple_extent_dims(dspace_id, dims, NULL);
        H5Sclose(dspace_id);
        /*extend dataset size*/
        size[0] = dims[0] + dims_ext_full_r[0];
        size[1] = dims[1];
        size[2] = dims[2];
        size[3] = dims[3];
        offset[0] = dims[0];
        offset[1] =  mpi_my_row_rank * parameters.nkx / mpi_dims[1];
        offset[2] =  0;
        offset[3] =  0;
        H5Dset_extent(dset_id, size);
        /*transform phi to real space*/
        fftw_copyFieldBuf_c(fftw_field, fields_fields.phi);
        fftw_c2r_field();
        fftw_transposeToXY_field();
        double *f_r = fftw_field;
        /*write phi real field to the file*/
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        dspace_id = H5Dget_space(dset_id);
        H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, stride,count, dims_ext_local_r);
        memspace = H5Screate_simple(3,&dims_ext_local_r[1],&max_dims_local_r[1]);
        if (mpi_my_col_rank != 0)
        {
            H5Sselect_none(dspace_id);
            H5Sselect_none(memspace);
        }
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, dspace_id, plist_id, f_r);
        H5Sclose(memspace);
        H5Sclose(dspace_id);
        H5Pclose(plist_id);
        H5Dclose(dset_id);
        // transposing data to original shape
        fftw_transposeToYX_field();

        /* writing B dataset */
        dset_id = H5Dopen2(file_id, "B", H5P_DEFAULT);
        dspace_id = H5Dget_space(dset_id);
        ndims = H5Sget_simple_extent_ndims(dspace_id);
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

        /*writing B real dataset*/
        /*opening a group and a dataset*/
        /*opening a group*/
        dset_id = H5Dopen2(file_id, "B_r", H5P_DEFAULT);
        /*open a dataset*/
        dspace_id = H5Dget_space(dset_id);
        /*get dataset's dimensions */
        ndims = H5Sget_simple_extent_ndims(dspace_id);
        H5Sget_simple_extent_dims(dspace_id, dims, NULL);
        H5Sclose(dspace_id);
        /*extend dataset size*/
        size[0] = dims[0] + dims_ext_full_r[0];
        size[1] = dims[1];
        size[2] = dims[2];
        size[3] = dims[3];
        offset[0] = dims[0];
        offset[1] =  mpi_my_row_rank * parameters.nkx / mpi_dims[1];
        offset[2] =  0;
        offset[3] =  0;
        H5Dset_extent(dset_id, size);
        /*transform B to real space*/
        fftw_copyFieldBuf_c(fftw_field, fields_fields.B);
        fftw_c2r_field();
        fftw_transposeToXY_field();
        f_r = fftw_field;
        /*write B real field to the file*/
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        dspace_id = H5Dget_space(dset_id);
        H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, stride,count, dims_ext_local_r);
        memspace = H5Screate_simple(3,&dims_ext_local_r[1],&max_dims_local_r[1]);
        if (mpi_my_col_rank != 0)
        {
            H5Sselect_none(dspace_id);
            H5Sselect_none(memspace);
        }
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, dspace_id, plist_id, f_r);
        H5Sclose(memspace);
        H5Sclose(dspace_id);
        H5Pclose(plist_id);
        H5Dclose(dset_id);
        // transposing data to original shape
        fftw_transposeToYX_field();

        /* writing A dataset */
        dset_id = H5Dopen2(file_id, "A", H5P_DEFAULT);
        dspace_id = H5Dget_space(dset_id);
        ndims = H5Sget_simple_extent_ndims(dspace_id);
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

        /*writing A real dataset*/
        /*opening a group and a dataset*/
        /*opening a group*/
        dset_id = H5Dopen2(file_id, "A_r", H5P_DEFAULT);
        /*open a dataset*/
        dspace_id = H5Dget_space(dset_id);
        /*get dataset's dimensions */
        ndims = H5Sget_simple_extent_ndims(dspace_id);
        H5Sget_simple_extent_dims(dspace_id, dims, NULL);
        H5Sclose(dspace_id);
        /*extend dataset size*/
        size[0] = dims[0] + dims_ext_full_r[0];
        size[1] = dims[1];
        size[2] = dims[2];
        size[3] = dims[3];
        offset[0] = dims[0];
        offset[1] =  mpi_my_row_rank * parameters.nkx / mpi_dims[1];
        offset[2] =  0;
        offset[3] =  0;
        H5Dset_extent(dset_id, size);
        /*transform A to real space*/
        fftw_copyFieldBuf_c(fftw_field, fields_fields.A);
        fftw_c2r_field();
        fftw_transposeToXY_field();
        f_r = fftw_field;
        /*write A real field to the file*/
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        dspace_id = H5Dget_space(dset_id);
        H5Sselect_hyperslab(dspace_id, H5S_SELECT_SET, offset, stride,count, dims_ext_local_r);
        memspace = H5Screate_simple(3,&dims_ext_local_r[1],&max_dims_local_r[1]);
        if (mpi_my_col_rank != 0)
        {
            H5Sselect_none(dspace_id);
            H5Sselect_none(memspace);
        }
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, dspace_id, plist_id, f_r);
        H5Sclose(memspace);
        H5Sclose(dspace_id);
        H5Pclose(plist_id);
        H5Dclose(dset_id);
        // transposing data to original shape
        fftw_transposeToYX_field();

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

    /*write a time dataset*/
    /*opening a group*/
    dset_id = H5Dopen2(file_id, "time", H5P_DEFAULT);
    /*open a dataset...*/
    dspace_id = H5Dget_space(dset_id);
    /*get dataset's dimensions */
    ndims = H5Sget_simple_extent_ndims(dspace_id);
    dims = malloc(ndims * sizeof(*dims));
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
        H5Dwrite(dset_id, H5T_NATIVE_DOUBLE,memspace, dspace_id,plist_id,&solver.curTime);
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
    for (size_t ii = 0; ii < array_local_size.total_comp; ii++){
        h[ii] = 0.;
    }
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
    int ndims = H5Sget_simple_extent_ndims(file_space);
    hsize_t *dims = malloc(ndims * sizeof(*dims));
    H5Sget_simple_extent_dims(file_space, dims, NULL);
    if (mpi_my_rank == 0){
        printf("DIMS OF DATASET = %zu,%zu,%zu,%zu,%zu,%zu\n",dims[0],dims[1],dims[2],dims[3],dims[4],dims[5]);
    }
    if(parameters.allow_rescale == 0){
        /*checking if dataset dimensions are the same as the simulation box resolution
         * if not, abort computation and exit*/
        if (array_global_size.nkx != dims[0]){
            if (mpi_my_rank == 0) printf("WARNING! nkx simulation resolution (%zu) not equal to file resolution (%zu)! Rescaling is not allowed! EXITING...\n", array_global_size.nkx, dims[0]);
            exit(1);
        }
        if (array_global_size.nky != dims[1]){
            if (mpi_my_rank == 0) printf("WARNING! nky simulation resolution (%zu) not equal to file resolution (%zu)! Rescaling is not allowed! EXITING...\n", array_global_size.nky, dims[1]);
            exit(1);
        }
        if (array_global_size.nkz != dims[2]){
            if (mpi_my_rank == 0) printf("WARNING! nkz simulation resolution (%zu) not equal to file resolution (%zu)! Rescaling is not allowed! EXITING...\n", array_global_size.nkz, dims[2]);
            exit(1);
        }
        if (array_global_size.nm != dims[3]){
            if (mpi_my_rank == 0) printf("WARNING! nm simulation resolution (%zu) not equal to file resolution (%zu)! Rescaling is not allowed! EXITING...\n", array_global_size.nm, dims[3]);
            exit(1);
        }
        /*read data if everything is fine with resolution*/
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset, stride, count, chunk_dims_c);
        memory_space = H5Screate_simple(hdf_rank, chunk_dims_c, NULL);

        /*
         * reading h dataset
         */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        H5Dread(dset_id,complex_id, memory_space, file_space, plist_id, h);
        /*close memory space*/
        H5Sclose(memory_space);
        H5Pclose(plist_id);
    }
    else{
        /*checking if dataset dimensions are the same as the simulation box resolution
         * if not, abort computation and exit*/
        if (array_global_size.nkx != dims[0]) {
            if (mpi_my_rank == 0)
                printf("WARNING! nkx simulation resolution (%zu) not equal to file resolution (%zu)! Rescaling will be performed!\n",
                       array_global_size.nkx, dims[0]);
        }
        if (array_global_size.nky != dims[1]) {
            if (mpi_my_rank == 0)
                printf("WARNING! nky simulation resolution (%zu) not equal to file resolution (%zu)! Rescaling will be performed!\n",
                       array_global_size.nky, dims[1]);
        }
        if (array_global_size.nkz != dims[2]) {
            if (mpi_my_rank == 0)
                printf("WARNING! nkz simulation resolution (%zu) not equal to file resolution (%zu)! Rescaling will be performed!\n",
                       array_global_size.nkz, dims[2]);
        }
        if (array_global_size.nm != dims[3]) {
            if (mpi_my_rank == 0)
                printf("WARNING! nm simulation resolution (%zu) not equal to file resolution (%zu)! Rescaling will be performed!\n",
                       array_global_size.nm, dims[3]);
        }
        /* we need take into account the normalization as well */
        COMPLEX old_norm = 1./(dims[0]*dims[1]*(dims[2] - 1) * 2);
        COMPLEX norm = 1./(fftw_norm/old_norm);

        /* next part is to read the correct hyperslab from the file.
         * If data in the file is bigger than simulation resolution, data will be truncated.
         * If it is smaller, than the data will be padded with zeros in the centre.
         * Algorithm is the following: first, data is treated as two separate arrays, one with positive kx and the other is with negative kx.
         * Array of positive kx is loaded to the left, and negative kx is loaded to the processors to the right.
         * Then they are being copied from buffer to a bigger/smaller array of data.*/
        hid_t offset_file[6];
        hid_t stride_file[6] = {1,1,1,1,1,1};
        hid_t count_file[6] = {1,1,1,1,1,1};
        hid_t chunk_dims_c_file[6];
        size_t yMax;
        size_t zMax;
        int procId_mMax;
        int procId_kMaxPos;
        int procId_kMaxNeg;
        int maxM;
        int maxXPositive;
        int maxXNegative;
        int maxXNegativeData;

        /*define maximum size of data in z and y*/
        yMax = (dims[1] > array_local_size.nky) ? array_local_size.nky : dims[1];
        zMax = (dims[2] > array_local_size.nkz) ? array_local_size.nkz : dims[2];

        /* first we want to find which m moment is the largest for the system or file, and at which processor it is located*/
        if(dims[3] > array_global_size.nm){
            maxM = array_global_size.nm - 1;
        }
        else{
            maxM = dims[3] - 1;
        }
        procId_mMax = mpi_whereIsM[2 * maxM];
        /* next we want to find out which is the processor where maximum positive and negative values of kx are located*/
        if(dims[0] > array_global_size.nkx){
            maxXPositive = array_global_size.nkx/2;
            maxXNegative = array_global_size.nkx/2 + 1;
            maxXNegativeData = dims[0] - (array_global_size.nkx/2) + 1;
        }
        else{
            maxXPositive = dims[0]/2;
            maxXNegative = array_global_size.nkx - dims[0]/2 + 1;
            maxXNegativeData = dims[0]/2 + 1;
        }
        procId_kMaxPos = mpi_whereIsX[2 * maxXPositive];
        procId_kMaxNeg = mpi_whereIsX[2 * maxXNegative];

        /* we now know processors id. Let's now set correct offsets to load data on processors (for now positive),
         * as well as array chunk dimensions.*/
        offset_file[0] = 0;
        offset_file[1] = 0;
        offset_file[2] = 0;
        offset_file[3] = 0;
        offset_file[4] = 0;
        offset_file[5] = 0;

        chunk_dims_c_file[0] = 0;
        chunk_dims_c_file[1] = dims[1];
        chunk_dims_c_file[2] = dims[2];
        chunk_dims_c_file[3] = 0;
        chunk_dims_c_file[4] = dims[4];
        chunk_dims_c_file[5] = dims[5];
        /* first we set offsets for m, as well as chunk sizes*/
        if(mpi_my_col_rank <= procId_mMax){
            offset_file[3] = mpi_my_col_rank * array_local_size.nm;
            chunk_dims_c_file[3] = array_local_size.nm;
            if (mpi_my_col_rank == procId_mMax) chunk_dims_c_file[3] = dims[3] - mpi_my_col_rank * array_local_size.nm;
        }
        else{
            offset_file[3] = 0;
            chunk_dims_c_file[3] = 0;
        }
        //printf("my_col_rank = %d, chunk_size for m = %zu, offset = %zu\n",mpi_my_col_rank, chunk_dims_c_file[3],offset_file[3]);

        /* setting offsets and chunk sizes for kx, positive array*/
        if(mpi_my_row_rank <= procId_kMaxPos){
            offset_file[0] = mpi_my_row_rank * array_local_size.nkx;
            chunk_dims_c_file[0] = array_local_size.nkx;
            if (mpi_my_row_rank == procId_kMaxPos) chunk_dims_c_file[0] = (maxXPositive + 1) - mpi_my_row_rank * array_local_size.nkx;
        }
        else{
            offset_file[0] = 0;
            chunk_dims_c_file[0] = 0;
        }

        /*now we are allocating buffers, in which the positive array will be stored*/
        size_t total_bufSize = chunk_dims_c_file[0] * chunk_dims_c_file[1]
                             * chunk_dims_c_file[2] * chunk_dims_c_file[3]
                             * chunk_dims_c_file[4] * chunk_dims_c_file[5];
        COMPLEX *buf = malloc(total_bufSize * sizeof(*buf));
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset_file, stride_file, count_file, chunk_dims_c_file);
        memory_space = H5Screate_simple(hdf_rank, chunk_dims_c_file, NULL);

        /*
         * reading h dataset
         */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        H5Dread(dset_id,complex_id, memory_space, file_space, plist_id, buf);

        /*fill positive part of array to data array*/
        for(size_t ix = 0; ix < chunk_dims_c_file[0]; ix++){
            //printf("ix = %zu\n",ix);
            for(size_t iy = 0; iy < yMax/2 + 1; iy++){
                for(size_t iz = 0; iz < zMax; iz++){
                    for(size_t im = 0; im < chunk_dims_c_file[3]; im++){
                        for(size_t il = 0; il < array_local_size.nl; il++){
                            for(size_t is = 0; is < array_local_size.ns; is++){
                                size_t indBufPosY = ix * chunk_dims_c_file[1] * chunk_dims_c_file[2] * chunk_dims_c_file[3] * chunk_dims_c_file[4] * chunk_dims_c_file[5]
                                                   +iy * chunk_dims_c_file[2] * chunk_dims_c_file[3] * chunk_dims_c_file[4] * chunk_dims_c_file[5]
                                                   +iz * chunk_dims_c_file[3] * chunk_dims_c_file[4] * chunk_dims_c_file[5]
                                                   +im * chunk_dims_c_file[4] * chunk_dims_c_file[5]
                                                   +il * chunk_dims_c_file[5]
                                                   +is;
                                size_t indArPosY = get_flat_c(is, il, im, ix, iy, iz);
                                h[indArPosY] = buf[indBufPosY];


                                if (iy != yMax/2){
                                    size_t indBufNegY = ix * chunk_dims_c_file[1] * chunk_dims_c_file[2] * chunk_dims_c_file[3] * chunk_dims_c_file[4] * chunk_dims_c_file[5]
                                                       +(chunk_dims_c_file[1] - iy - 1) * chunk_dims_c_file[2] * chunk_dims_c_file[3] * chunk_dims_c_file[4] * chunk_dims_c_file[5]
                                                       +iz * chunk_dims_c_file[3] * chunk_dims_c_file[4] * chunk_dims_c_file[5]
                                                       +im * chunk_dims_c_file[4] * chunk_dims_c_file[5]
                                                       +il * chunk_dims_c_file[5]
                                                       +is;
                                    size_t indArNegY = get_flat_c(is, il, im, ix, array_local_size.nky - iy - 1, iz);
                                    h[indArNegY] = buf[indBufNegY];
                                }
                            }
                        }
                    }
                }
            }
        }

        /*close memory space*/
        H5Sclose(memory_space);
        H5Pclose(plist_id);
        free(buf);

        /*Now we need to load negative k array. To do so, we first need to set correct negative offsets and chunk sizes. Only 0th dimension is required,
         * others are not changed*/
        offset_file[0] = 0;
        chunk_dims_c_file[0] = 0;
        int negativeOffsetStart = 0;
        if (mpi_my_row_rank == procId_kMaxNeg){
            offset_file[0] = maxXNegativeData;
            int local_kx = mpi_whereIsX[2 * (maxXNegative) + 1];
            chunk_dims_c_file[0] = array_local_size.nkx - local_kx;
            negativeOffsetStart = chunk_dims_c_file[0] + offset_file[0];
        }

        MPI_Bcast(&negativeOffsetStart,1,MPI_INT,procId_kMaxNeg,mpi_row_comm);
        if(mpi_my_row_rank > procId_kMaxNeg){
            offset_file[0] =  negativeOffsetStart + (mpi_my_row_rank - procId_kMaxNeg - 1) * array_local_size.nkx;
            chunk_dims_c_file[0] = array_local_size.nkx;
        }

        /* we have computed offsets and chunk dimensions, now we compute total size of buffer,
         * select the correct hyperslab and prepare memory space for it*/
        total_bufSize = chunk_dims_c_file[0] * chunk_dims_c_file[1]
                        * chunk_dims_c_file[2] * chunk_dims_c_file[3]
                        * chunk_dims_c_file[4] * chunk_dims_c_file[5];
        buf = malloc(total_bufSize * sizeof(*buf));
        H5Sselect_hyperslab(file_space, H5S_SELECT_SET, offset_file, stride_file, count_file, chunk_dims_c_file);
        memory_space = H5Screate_simple(hdf_rank, chunk_dims_c_file, NULL);

        /*
         * reading h dataset
         */
        plist_id = H5Pcreate(H5P_DATASET_XFER);
        H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
        H5Dread(dset_id,complex_id, memory_space, file_space, plist_id, buf);

        /*fill negative part of array to data array*/
        size_t kxNeg;
        for(size_t ix = 0; ix < chunk_dims_c_file[0]; ix++){
            //printf("ix = %zu\n",ix);
            if(mpi_my_row_rank == procId_kMaxNeg){
                kxNeg = mpi_whereIsX[2 * (maxXNegative + ix) + 1];
                //printf("kxNeg = %zu\n",kxNeg);
            }
            else{
                kxNeg = ix;
            }
            for(size_t iy = 0; iy < yMax/2 + 1; iy++){
                for(size_t iz = 0; iz < zMax; iz++){
                    for(size_t im = 0; im < chunk_dims_c_file[3]; im++){
                        for(size_t il = 0; il < array_local_size.nl; il++){
                            for(size_t is = 0; is < array_local_size.ns; is++){
                                size_t indBufPosY = ix * chunk_dims_c_file[1] * chunk_dims_c_file[2] * chunk_dims_c_file[3] * chunk_dims_c_file[4] * chunk_dims_c_file[5]
                                                    +iy * chunk_dims_c_file[2] * chunk_dims_c_file[3] * chunk_dims_c_file[4] * chunk_dims_c_file[5]
                                                    +iz * chunk_dims_c_file[3] * chunk_dims_c_file[4] * chunk_dims_c_file[5]
                                                    +im * chunk_dims_c_file[4] * chunk_dims_c_file[5]
                                                    +il * chunk_dims_c_file[5]
                                                    +is;
                                size_t indArPosY = get_flat_c(is, il, im, kxNeg, iy, iz);
                                h[indArPosY] = buf[indBufPosY];


                                if (iy != yMax/2){
                                    size_t indBufNegY = ix * chunk_dims_c_file[1] * chunk_dims_c_file[2] * chunk_dims_c_file[3] * chunk_dims_c_file[4] * chunk_dims_c_file[5]
                                                        +(chunk_dims_c_file[1] - iy - 1) * chunk_dims_c_file[2] * chunk_dims_c_file[3] * chunk_dims_c_file[4] * chunk_dims_c_file[5]
                                                        +iz * chunk_dims_c_file[3] * chunk_dims_c_file[4] * chunk_dims_c_file[5]
                                                        +im * chunk_dims_c_file[4] * chunk_dims_c_file[5]
                                                        +il * chunk_dims_c_file[5]
                                                        +is;
                                    size_t indArNegY = get_flat_c(is, il, im, kxNeg, array_local_size.nky - iy - 1, iz);
                                    h[indArNegY] = buf[indBufNegY];
                                }
                            }
                        }
                    }
                }
            }
        }

        /*close memory space*/
        H5Sclose(memory_space);
        H5Pclose(plist_id);
        free(buf);

        /*normalize data correctly*/
        for(size_t ii = 0; ii < array_local_size.total_comp;ii++){
            h[ii] *= norm;
        }
    }
    /*closing dataset identifier and file space*/
    H5Dclose(dset_id);
    H5Sclose(file_space);



};