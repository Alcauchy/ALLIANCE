//
// Created by alcauchy on 16/12/2021.
//

#ifndef ALLIANCE_ALPHA_1_0_HDF_UTILS_H
#define ALLIANCE_ALPHA_1_0_HDF_UTILS_H
#include <stdio.h> // standard C lib
#include <stdlib.h> //Standard C lib
#include <hdf5.h>
#include <hdf5_hl.h>
#include "mpi_utils.h"
#include "array.h"
#include "fields.h"
#include "init.h"

void hdf_init();                                            // inits hdf5 environment
void hdf_create_file_c(char *filename, COMPLEX *data);      // creates file with complex data
void hdf_create_file_r(char *filename, double *data);       // creates file with real data
void hdf_initField();
void hdf_saveFieldA(char *filename);
void hdf_saveFieldB(char *filename);
void hdf_saveFieldPhi(char *filename);
void hdf_saveData(COMPLEX *h, int timestep);
void hdf_createFiles();
void hdf_saveKSpec(int timestep);
void hdf_saveMSpec(int timestep);
void hdf_createFieldFile();
void hdf_saveFields(int timestep);
void hdf_createCheckpoint(COMPLEX *h, int timestep);
void hdf_initCheckpoints();
void hdf_dumpCheckpoint(COMPLEX *h, int timestep, char *filename);
void hdf_readData(char *filename, COMPLEX *h);
void hdf_saveDistrib(COMPLEX* h, int timestep);
typedef struct {                                            // structure needed to define the custom complex datatype. This datatype is then used to write complex data into the file
    double re;   //real part
    double im;   //imaginary part
} complex_t;






#endif //ALLIANCE_ALPHA_1_0_HDF_UTILS_H
