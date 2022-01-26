//
// Created by alcauchy on 20/01/2022.
//

//
// Created by alcauchy on 20/01/2022.
//
#include "diagnostics.h"
#include <time.h>

void diag_multiply_ar_c(COMPLEX *ar1, COMPLEX *ar2, COMPLEX *ret){
    for(size_t ikx = 0; ikx < array_local_size.nkx; ikx++) {
        for (size_t iky = 0; iky < array_local_size.nky; iky++) {
            for (size_t ikz = 0; ikz < array_local_size.nkz; ikz++) {
                for (size_t im = 0; im < array_local_size.nm; im++) {
                    for (size_t il = 0; il < array_local_size.nl; il++) {
                        for (size_t is = 0; is < array_local_size.ns; is++) {
                            ret[get_flat_c(is,il,im,ikx,iky,ikz)] = ar1[get_flat_c(is,il,im,ikx,iky,ikz)]*ar2[get_flat_c(is,il,im,ikx,iky,ikz)];
                        }
                    }
                }
            }
        }
    }
}

void diag_multiply_ar_r(double *ar1, double *ar2, double *ret){
    for(size_t ikx = 0; ikx < array_local_size.nkx; ikx++) {
        for (size_t iky = 0; iky < array_local_size.nky; iky++) {
            for (size_t ikz = 0; ikz < array_local_size.nz+2; ikz++) {
                for (size_t im = 0; im < array_local_size.nm; im++) {
                    for (size_t il = 0; il < array_local_size.nl; il++) {
                        for (size_t is = 0; is < array_local_size.ns; is++) {
                            ret[get_flat_r(is,il,im,ikx,iky,ikz)] = ar1[get_flat_r(is,il,im,ikx,iky,ikz)]*ar2[get_flat_r(is,il,im,ikx,iky,ikz)];
                        }
                    }
                }
            }
        }
    }
}

void diag_fill_rand(COMPLEX *ar1){
    srand(time(NULL));
    for(size_t ikx = 0; ikx < array_local_size.nkx; ikx++) {
        for (size_t iky = 0; iky < array_local_size.nky; iky++) {
            for (size_t ikz = 0; ikz < array_local_size.nkz; ikz++) {
                for (size_t im = 0; im < array_local_size.nm; im++) {
                    for (size_t il = 0; il < array_local_size.nl; il++) {
                        for (size_t is = 0; is < array_local_size.ns; is++) {
                            ar1[get_flat_c(is,il,im,ikx,iky,ikz)] = (0.5 - (double)rand()/(double)(RAND_MAX))+(0.5 - (double)rand()/(double)(RAND_MAX))*1.j;
                        }
                    }
                }
            }
        }
    }
}