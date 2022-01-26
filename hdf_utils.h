//
// Created by alcauchy on 16/12/2021.
//

#ifndef ALLIANCE_ALPHA_1_0_HDF_UTILS_H
#define ALLIANCE_ALPHA_1_0_HDF_UTILS_H
#include "master.h"

void hdf_init();                                            // inits hdf5 environment
void hdf_create_file_c(char *filename, COMPLEX *data);      // creates file with complex data
void hdf_create_file_r(char *filename, double *data);       // creates file with real data

typedef struct {                                            // structure needed to define the custom complex datatype. This datatype is then used to write complex data into the file
    double re;   //real part
    double im;   //imaginary part
} complex_t;






#endif //ALLIANCE_ALPHA_1_0_HDF_UTILS_H
