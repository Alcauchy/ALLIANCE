//
// Created by alcauchy on 16/12/2021.
//

#include "hdf_utils.h"


int hdf_rank = 6;
hid_t complex_id;
hsize_t dataspace_dims_r[6];
hsize_t dataspace_dims_c[6];
hsize_t chunk_dims_r[6];
hsize_t chunk_dims_c[6];
hsize_t offset[6];
hsize_t count[6] = {1,1,1,1,1,1};
hsize_t stride[6] = {1,1,1,1,1,1};
herr_t	status;
MPI_Info info = MPI_INFO_NULL;
complex_t tmp;  //used only to compute offsets

void complex_t_init(){
    complex_id = H5Tcreate(H5T_COMPOUND, sizeof(tmp));
    H5Tinsert(complex_id, "r", HOFFSET(complex_t,re), H5T_NATIVE_DOUBLE);
    H5Tinsert(complex_id, "i", HOFFSET(complex_t,im), H5T_NATIVE_DOUBLE);
}

void hdf_init(){
    complex_t_init();
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

    printf("[MPI process %d] my coords are (%d, %d), offsets are (%d,%d,%d,%d,%d,%d)\n",mpi_my_rank,mpi_my_coords[0],mpi_my_coords[1],offset[0],offset[1],offset[2],offset[3],offset[4],offset[5]);
    printf("[MPI process %d] my coords are (%d, %d), chunk dims for real are (%d,%d,%d,%d,%d,%d)\n",mpi_my_rank,mpi_my_coords[0],mpi_my_coords[1],chunk_dims_r[0],chunk_dims_r[1],chunk_dims_r[2],chunk_dims_r[3],chunk_dims_r[4],chunk_dims_r[5]);
    printf("[MPI process %d] my coords are (%d, %d), chunk dims for complex are (%d,%d,%d,%d,%d,%d)\n",mpi_my_rank,mpi_my_coords[0],mpi_my_coords[1],chunk_dims_c[0],chunk_dims_c[1],chunk_dims_c[2],chunk_dims_c[3],chunk_dims_c[4],chunk_dims_c[5]);
};

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